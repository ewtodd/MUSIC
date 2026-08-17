#include "StripSumScatter.hpp"
#include <fstream>
#include <queue>
#include <thread>

void StripSumScatter::FillScatters(const std::vector<Int_t> &runOrder,
                                   std::map<Int_t, TChain *> &chains) {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kXLo = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_HI;
  const Int_t kXBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.XBINS;
  const Double_t kXMin = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.XMIN;
  const Double_t kXMax = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.XMAX;
  const Int_t kYBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.YBINS;
  const Int_t kBeamReservoirCap =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRACES_PER_CLASS * 10;
  const Int_t nReac = kReacMax - kReacMin + 1;

  // Compute y-bounds once.
  YBounds(m_yLo, m_yHi);

  // Allocate scatter histograms.
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    Int_t ri = ReacIndex(reac);
    TH2F *h = new TH2F(
        Form("scatter_r%d", reac),
        Form(";norm. #DeltaE strips %d#rightarrow%d [a.u.];norm. #DeltaE "
             "strips %d#rightarrow%d [a.u.]",
             kXLo, kXHi, YLoOf(reac), YHiOf(reac)),
        kXBins, kXMin, kXMax, kYBins, m_yLo[ri], m_yHi[ri]);
    h->SetDirectory(nullptr);
    h->SetStats(0);
    m_scatter[reac] = h;
  }

  std::vector<GateSpec> activeGates = ActiveGates();

  // Pre-index chains so worker threads only touch per-run entries.
  const Int_t nRuns = Int_t(runOrder.size());
  std::vector<TChain *> chainVec(nRuns);
  for (Int_t i = 0; i < nRuns; i++)
    chainVec[i] = chains[runOrder[i]];

  Int_t n_workers =
      TMath::Min(Int_t(std::thread::hardware_concurrency()), nRuns);
  n_workers = TMath::Min(
      n_workers, Constants::cfg.STRIP_SUM_SCATTER_CONFIG.MAX_STRIP_SUM_WORKERS);
  if (n_workers < 1)
    n_workers = 1;
  std::cout << "strip-sum-scatter: " << nRuns << " runs on " << n_workers
            << " workers" << std::endl;
  std::cout
      << "  filters: noise(reject="
      << Int_t(Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REJECT_NOISE)
      << ", <= " << Constants::cfg.STRIP_SUM_SCATTER_CONFIG.NOISE_THRESHOLD
      << " on >= " << Constants::cfg.STRIP_SUM_SCATTER_CONFIG.NOISE_MIN_STRIPS
      << " strips) pileup(reject="
      << Int_t(Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REJECT_PILEUP)
      << ", >= " << Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PILEUP_THRESHOLD
      << " on >= " << Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PILEUP_MIN_STRIPS
      << " strips) offbeam(reject="
      << Int_t(Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REJECT_OFFBEAM)
      << ", dist=" << Constants::cfg.STRIP_SUM_SCATTER_CONFIG.OFFBEAM_DIST
      << " on >= " << Constants::cfg.STRIP_SUM_SCATTER_CONFIG.OFFBEAM_MIN_STRIPS
      << " strips) high-strip(reject="
      << Int_t(Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REJECT_HIGH_STRIP)
      << ", any strip > "
      << Constants::cfg.STRIP_SUM_SCATTER_CONFIG.HIGH_STRIP_THRESHOLD << ")"
      << std::endl;

  // --- Phase 1: beam classification ellipses + scatter filter gates ---
  // One worker task per run; within a run the series gates stay sequential
  // (each gate only sees events passing all prior gates). The per-gate
  // plots are serialized inside FindBeamGate by g_plot_mutex.
  std::vector<RunGateFit> fits(nRuns);
  {
    std::queue<Int_t> work;
    for (Int_t i = 0; i < nRuns; i++)
      work.push(i);
    std::mutex work_mutex;
    std::vector<std::thread> workers;
    for (Int_t w = 0; w < n_workers; w++) {
      workers.emplace_back([&]() {
        while (true) {
          Int_t i;
          {
            std::lock_guard<std::mutex> lk(work_mutex);
            if (work.empty())
              return;
            i = work.front();
            work.pop();
          }
          TChain *ch = chainVec[i];
          if (!ch || ch->GetEntries() == 0)
            continue;
          fits[i] = FitRunGates(runOrder[i], ch);
        }
      });
    }
    for (Int_t w = 0; w < Int_t(workers.size()); w++)
      workers[w].join();
  }

  std::map<Int_t, BeamEllipses> beamEllipses;
  std::map<Int_t, std::vector<BeamFit2D>> gates;
  std::map<Int_t, Bool_t> ok;
  for (Int_t i = 0; i < nRuns; i++) {
    Int_t run = runOrder[i];
    beamEllipses[run] = fits[i].beam;
    gates[run] = fits[i].gates;
    ok[run] = fits[i].ok;
  }

  // --- Phase 2: per-run event filling (parallel across runs) ---
  // Each worker fills private scatter histograms + a private reservoir
  // slice; both are merged below in run order so the combined result is
  // identical to a single-threaded fill.
  std::vector<FillRunResult> fills(nRuns);
  {
    std::queue<Int_t> work;
    for (Int_t i = 0; i < nRuns; i++)
      work.push(i);
    std::mutex work_mutex;
    std::vector<std::thread> workers;
    for (Int_t w = 0; w < n_workers; w++) {
      workers.emplace_back([&]() {
        while (true) {
          Int_t i;
          {
            std::lock_guard<std::mutex> lk(work_mutex);
            if (work.empty())
              return;
            i = work.front();
            work.pop();
          }
          Int_t run = runOrder[i];
          TChain *ch = chainVec[i];
          if (!ch || !ok[run])
            continue;
          fills[i] = FillRunScatters(run, ch, activeGates, gates[run],
                                     beamEllipses[run], m_yLo, m_yHi);
        }
      });
    }
    for (Int_t w = 0; w < Int_t(workers.size()); w++)
      workers[w].join();
  }

  // Merge per-run scatter histograms into m_scatter. Runs that were skipped
  // (no chain / no valid gates) never produced scatters -- their result was
  // left default-constructed with an empty vector, so skip them here.
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    Int_t ri = ReacIndex(reac);
    for (Int_t i = 0; i < nRuns; i++) {
      if (Int_t(fills[i].scatters.size()) != nReac)
        continue;
      TH2F *h = fills[i].scatters[ri];
      if (!h)
        continue;
      m_scatter[reac]->Add(h);
      delete h;
      fills[i].scatters[ri] = nullptr;
    }
  }

  // Merge reservoirs in run order. Every reaction event is kept; pure-beam
  // events are kept only while the global cap is not yet reached (the same
  // selection a single-threaded fill would make).
  Long64_t totalGated = 0, totalSeen = 0;
  Long64_t nRejGate = 0, nRejPileup = 0, nRejNoise = 0, nRejHighStrip = 0,
           nRejOffbeam = 0;
  Int_t nBeamKept = 0;
  m_reservoir.clear();
  for (Int_t i = 0; i < nRuns; i++) {
    totalSeen += fills[i].seen;
    totalGated += fills[i].gated;
    nRejGate += fills[i].n_rej_gate;
    nRejPileup += fills[i].n_rej_pileup;
    nRejNoise += fills[i].n_rej_noise;
    nRejHighStrip += fills[i].n_rej_highstrip;
    nRejOffbeam += fills[i].n_rej_offbeam;
    for (Int_t k = 0; k < Int_t(fills[i].reservoir.size()); k++) {
      TraceEvt e = fills[i].reservoir[k];
      if (e.reac_mask != 0)
        m_reservoir.push_back(e);
      else if (nBeamKept < kBeamReservoirCap) {
        m_reservoir.push_back(e);
        nBeamKept++;
      }
    }
    std::vector<TraceEvt>().swap(fills[i].reservoir);
  }

  std::cout << "Built scatters: " << totalGated
            << " reaction events across strips " << kReacMin << "-" << kReacMax
            << " (" << totalSeen << " seen), reservoir " << m_reservoir.size()
            << " events." << std::endl;
  std::cout << "  rejected: " << nRejGate << " gate, " << nRejPileup
            << " pileup, " << nRejNoise << " noise, " << nRejHighStrip
            << " high-strip, " << nRejOffbeam << " offbeam" << std::endl;
}

StripSumScatter::RunGateFit StripSumScatter::FitRunGates(Int_t run,
                                                         TChain *chain) {
  RunGateFit out;
  out.beam.ok = kFALSE;

  // --- Beam classification ellipses (no gating) ---
  std::vector<GateSpec> emptyPrior;
  std::vector<BeamFit2D> emptyGates;
  const TString tag = Form("run%d", run);
  const TString subdir = Form("strip_sum_scatter/run%d", run);
  Int_t ent_sx = 0, ent_sy = 1;
  const Char_t *ent_tag = "s0/s1";
  if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PURE_BEAM_GATE ==
      StripSumScatterConfig::PURE_BEAM_GATE_S1_S2) {
    ent_sx = 1;
    ent_sy = 2;
    ent_tag = "s1/s2";
  }
  BeamFit2D ent_ell =
      FindBeamGate(chain, ent_sx, ent_sy, emptyPrior, emptyGates, tag, subdir);
  if (ent_sx == 0)
    out.beam.s0_s1 = ent_ell;
  else
    out.beam.s1_s2 = ent_ell;
  if (ent_ell.ok)
    std::cout << "  run " << run << " beam ellipse " << ent_tag << ": mu=("
              << ent_ell.mu_x << "," << ent_ell.mu_y << ")" << std::endl;
  else {
    std::cerr << "  run " << run << " beam ellipse " << ent_tag
              << " failed; skipping run" << std::endl;
    return out;
  }
  if (Constants::cfg.IGNORE_STRIP_17) {
    out.beam.use_s15_s16 = kTRUE;
    out.beam.s15_s16 =
        FindBeamGate(chain, 15, 16, emptyPrior, emptyGates, tag, subdir);
    if (out.beam.s15_s16.ok)
      std::cout << "  run " << run << " beam ellipse s15/s16: mu=("
                << out.beam.s15_s16.mu_x << "," << out.beam.s15_s16.mu_y << ")"
                << std::endl;
    else {
      std::cerr << "  run " << run
                << " beam ellipse s15/s16 failed; skipping run" << std::endl;
      return out;
    }
  } else {
    out.beam.use_s15_s16 = kFALSE;
    out.beam.s16_s17 =
        FindBeamGate(chain, 16, 17, emptyPrior, emptyGates, tag, subdir);
    if (out.beam.s16_s17.ok)
      std::cout << "  run " << run << " beam ellipse s16/s17: mu=("
                << out.beam.s16_s17.mu_x << "," << out.beam.s16_s17.mu_y << ")"
                << std::endl;
    else {
      std::cerr << "  run " << run
                << " beam ellipse s16/s17 failed; skipping run" << std::endl;
      return out;
    }
  }
  out.beam.ok = kTRUE;

  // --- Scatter filter gates (series gating) ---
  // Each gate only sees events passing all prior gates; a failed gate is
  // still recorded so later gates keep their prior list consistent.
  std::vector<GateSpec> activeGates = ActiveGates();
  std::vector<BeamFit2D> runGates;
  std::vector<GateSpec> priorSpecs;
  Bool_t allOk = kTRUE;
  for (Int_t gi = 0; gi < activeGates.size(); gi++) {
    BeamFit2D g = FindBeamGate(chain, activeGates[gi].sx, activeGates[gi].sy,
                               priorSpecs, runGates, tag, subdir);
    if (g.ok)
      std::cout << "  run " << run << " beam gate s" << activeGates[gi].sx
                << "/s" << activeGates[gi].sy << ": mu=(" << g.mu_x << ","
                << g.mu_y << ")" << std::endl;
    else {
      std::cerr << "  run " << run << " beam gate s" << activeGates[gi].sx
                << "/s" << activeGates[gi].sy << " failed; skipping run"
                << std::endl;
      allOk = kFALSE;
    }
    runGates.push_back(g);
    priorSpecs.push_back(activeGates[gi]);
  }
  out.gates = runGates;
  out.ok = allOk;
  return out;
}

StripSumScatter::FillRunResult StripSumScatter::FillRunScatters(
    Int_t run, TChain *chain, const std::vector<GateSpec> &active_gates,
    const std::vector<BeamFit2D> &run_gates, const BeamEllipses &run_beam,
    const Double_t *y_lo, const Double_t *y_hi) {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kXLo = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_HI;
  const Int_t kXBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.XBINS;
  const Double_t kXMin = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.XMIN;
  const Double_t kXMax = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.XMAX;
  const Int_t kYBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.YBINS;
  const Int_t kBeamReservoirCap =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRACES_PER_CLASS * 10;

  FillRunResult out;
  const Int_t nReac = kReacMax - kReacMin + 1;
  out.scatters.resize(nReac);
  for (Int_t ri = 0; ri < nReac; ri++) {
    Int_t reac = kReacMin + ri;
    out.scatters[ri] = new TH2F(
        Form("scatter_r%d", reac),
        Form(";norm. #DeltaE strips %d#rightarrow%d [a.u.];norm. #DeltaE "
             "strips %d#rightarrow%d [a.u.]",
             kXLo, kXHi, YLoOf(reac), YHiOf(reac)),
        kXBins, kXMin, kXMax, kYBins, y_lo[ri], y_hi[ri]);
    out.scatters[ri]->SetDirectory(nullptr);
    out.scatters[ri]->SetStats(0);
  }

  EnergyView ev;
  ev.Attach(chain);
  EnableEventBranches(chain);
  Long64_t n = chain->GetEntries();
  std::cout << "Run " << run << ": filling " << nReac
            << " reaction-strip scatters over " << n << " events..."
            << std::endl;

  Int_t beam_kept = 0;
  for (Long64_t j = 0; j < n; j++) {
    chain->GetEntry(j);
    ev.Decode();
    out.seen++;

    Bool_t passesAll = kTRUE;
    for (Int_t gi = 0; gi < active_gates.size(); gi++)
      if (!PassesGate(run_gates[gi], ev, active_gates[gi].sx,
                      active_gates[gi].sy)) {
        passesAll = kFALSE;
        break;
      }
    if (!passesAll) {
      out.n_rej_gate++;
      continue;
    }
    if (IsPileup(ev)) {
      out.n_rej_pileup++;
      continue;
    }
    if (IsNoise(ev)) {
      out.n_rej_noise++;
      continue;
    }
    if (IsHighStrip(ev)) {
      out.n_rej_highstrip++;
      continue;
    }
    if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REJECT_OFFBEAM &&
        IsOffbeam(ev)) {
      out.n_rej_offbeam++;
      continue;
    }

    Double_t x = SumRange(ev.total, kXLo, kXHi);
    UInt_t mask = 0;
    for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
      if (!PassesReaction(ev, reac))
        continue;
      mask |= (1u << ReacIndex(reac));
      out.scatters[ReacIndex(reac)]->Fill(
          x, SumRange(ev.total, YLoOf(reac), YHiOf(reac)));
    }

    // Keep every reaction-passing event for traces; cap pure-beam events
    // (only ~TRACES_PER_CLASS are ever drawn). The two are mutually
    // exclusive -- a pure-beam event has no reaction jump.
    Bool_t beam = (mask == 0) && IsPureBeam(ev, run_beam);
    if (mask == 0 && !(beam && beam_kept < kBeamReservoirCap))
      continue;
    if (beam)
      beam_kept++;

    TraceEvt e;
    for (Int_t s = 0; s < 18; s++) {
      e.total[s] = Float_t(ev.total[s]);
      e.total_adc[s] =
          Float_t(ev.left_0_17_adc[s]) + Float_t(ev.rightdE_adc[s]);
    }
    // Mirror IGNORE_SHORT_STRIPS: the normed total keeps only the long side
    // of a split strip, so the raw trace must drop the same side to stay
    // comparable.
    if (Constants::cfg.IGNORE_SHORT_STRIPS)
      for (Int_t s = 1; s <= 16; s++)
        e.total_adc[s] = ((s % 2) != 0) ? Float_t(ev.left_0_17_adc[s])
                                        : Float_t(ev.rightdE_adc[s]);
    // Both-channel multiplicity: split strips (1-16) where both ends
    // FIRED. Read off the RAW ADC, not the calibrated ends -- the
    // short-end gains are 0 (uncalibrated, no sim anchor), so the
    // calibrated short ends are always zero; the raw ADC still carries
    // whether the channel fired.
    Int_t both = 0;
    for (Int_t s = 1; s <= 16; s++)
      if (ev.left_0_17_adc[s] > 0.0 && ev.rightdE_adc[s] > 0.0)
        both++;
    e.both_mult = both;
    e.reac_mask = mask;
    e.beam_flat = beam;
    out.reservoir.push_back(e);
    if (mask != 0)
      out.gated++;
  }
  return out;
}

void StripSumScatter::PlotScatters() {
  const Int_t kXLo = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_HI;
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;

  std::lock_guard<std::mutex> lock(g_plot_mutex);
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(m_scatter[reac], c);
    m_scatter[reac]->GetYaxis()->SetTitleOffset(1.3);
    c->SetLeftMargin(0.18);
    PlottingUtils::SaveFigure(c,
                              Form("normsumE_reac%d_s%d_%d_vs_s%d_%d", reac,
                                   YLoOf(reac), YHiOf(reac), kXLo, kXHi),
                              "strip_sum_scatter", PlotSaveOptions::kLINEAR);
    delete c;
  }
}

void StripSumScatter::InteractiveOverlay(Int_t reac) {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kXLo = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_HI;
  const Int_t kTracesPerRegion =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRACES_PER_CLASS;

  if (reac < kReacMin || reac > kReacMax) {
    std::cerr << "strip-sum-scatter: candidate reaction strip " << reac
              << " outside [" << kReacMin << "," << kReacMax
              << "]; skipping interactive overlay." << std::endl;
    return;
  }
  if (!gSystem->Getenv("DISPLAY")) {
    std::cerr << "strip-sum-scatter: no DISPLAY; skipping interactive "
                 "region-trace overlay (scatters already saved)."
              << std::endl;
    return;
  }

  Int_t app_argc = 1;
  char app_arg0[] = "strip-sum-scatter";
  char *app_argv[] = {app_arg0};
  TApplication app("strip-sum-scatter", &app_argc, app_argv);
  gROOT->SetBatch(kFALSE);

  // Swallow recoverable X protocol errors for the interactive cut session so
  // the cut-canvas teardown's BadWindow/BadDrawable doesn't trip ROOT's
  // crashing default handler at the next canvas paint (same guard the fit
  // editors use).
  const AUXErrorHandlerSave xerr_save = AUInstallTolerantXErrorHandler();

  TCanvas *cutCanvas = new TCanvas("c_strip_sum_regions",
                                   "Draw (a,n) then (a,a') regions", 900, 700);
  cutCanvas->SetLogz(kTRUE); // match the saved scatter's z-scale
  m_scatter[reac]->Draw("COLZ");
  cutCanvas->Update();
  TCutG *cutAn = PromptCut(cutCanvas, "region_an", "(a,n)");
  TCutG *cutAa = PromptCut(cutCanvas, "region_aa", "(a,a')");

  cutCanvas->GetListOfPrimitives()->Remove(cutAn);
  cutCanvas->GetListOfPrimitives()->Remove(cutAa);
  gROOT->SetEditorMode();
  delete cutCanvas;
  gSystem->ProcessEvents();
  // SetBatch only flips the flag; the real X11 backend (TGX11) stays active,
  // and batch-mode canvas creation against it hits a null window context in
  // TGX11::DrawBoxW (fWindows[-1] is a null map entry). Re-point gVirtualX at
  // the no-op batch backend, the same switch TApplication::MakeBatch() makes
  // (MakeBatch itself is protected). The interactive TGX11 is deliberately
  // leaked rather than deleted to avoid tearing down the X display while the
  // TApplication still lives.
  gROOT->SetBatch(kTRUE);
  gVirtualX = gGXBatch;

  ClusterVarHists(reac, cutAa, cutAn, "strip_sum_scatter");

  std::vector<TGraph *> tr_an, tr_aa, tr_beam;
  // Same selected events, raw (un-normalized) ADC -- one entry per normed
  // trace, kept in lock-step so the two overlays show the identical events.
  std::vector<TGraph *> tr_an_adc, tr_aa_adc, tr_beam_adc;
  // Same selected events again, Savitzky-Golay smoothed (normed a.u. space),
  // for the with-smoothing overlay -- also in lock-step with the raw normed
  // traces, so the two a.u. overlays show the identical events.
  const Bool_t kSkipSg =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.SKIP_SAVGOL_PLOTS;
  std::vector<TGraph *> tr_an_sg, tr_aa_sg, tr_beam_sg;
  UInt_t bit = (1u << ReacIndex(reac));

  for (Int_t k = 0; k < m_reservoir.size(); k++) {
    if (Int_t(tr_an.size()) >= kTracesPerRegion &&
        Int_t(tr_aa.size()) >= kTracesPerRegion &&
        Int_t(tr_beam.size()) >= kTracesPerRegion)
      break;
    const TraceEvt &e = m_reservoir[k];
    if (e.beam_flat && Int_t(tr_beam.size()) < kTracesPerRegion) {
      tr_beam.push_back(TraceFromTotal(e.total));
      tr_beam_adc.push_back(TraceFromTotal(e.total_adc));
      if (!kSkipSg)
        tr_beam_sg.push_back(SmoothedTraceFromTotal(e.total));
      continue;
    }
    if (!(e.reac_mask & bit))
      continue;
    Double_t td[18];
    for (Int_t s = 0; s < 18; s++)
      td[s] = Double_t(e.total[s]);
    Double_t x = SumRange(td, kXLo, kXHi);
    Double_t y = SumRange(td, YLoOf(reac), YHiOf(reac));
    if (cutAn && Int_t(tr_an.size()) < kTracesPerRegion &&
        cutAn->IsInside(x, y)) {
      tr_an.push_back(TraceFromTotal(e.total));
      tr_an_adc.push_back(TraceFromTotal(e.total_adc));
      if (!kSkipSg)
        tr_an_sg.push_back(SmoothedTraceFromTotal(e.total));
    } else if (cutAa && Int_t(tr_aa.size()) < kTracesPerRegion &&
               cutAa->IsInside(x, y)) {
      tr_aa.push_back(TraceFromTotal(e.total));
      tr_aa_adc.push_back(TraceFromTotal(e.total_adc));
      if (!kSkipSg)
        tr_aa_sg.push_back(SmoothedTraceFromTotal(e.total));
    }
  }

  std::cout << "Sampled traces: beam=" << tr_beam.size()
            << " (a,a')=" << tr_aa.size() << " (a,n)=" << tr_an.size()
            << std::endl;

  // Dump the exact sampled trace values (plus per-trace count of strips
  // over the pileup threshold) next to the plot, so filter behavior can be
  // checked against the raw numbers instead of eyeballing the PNG. The
  // graph point k holds strip (s_lo + k), where s_lo follows the
  // IGNORE_STRIP_0/17 config like the plots themselves.
  {
    TString dump_path = Paths::ResultsDir() + "/plots/strip_sum_scatter/" +
                        Form("region_traces_reac%d.txt", reac);
    const Int_t s_lo = Constants::cfg.IGNORE_STRIP_0 ? 1 : 0;
    const Int_t s_hi = Constants::cfg.IGNORE_STRIP_17 ? 16 : 17;
    std::ofstream os(dump_path.Data());
    os << "class trace strips_ge_1.3 total[" << s_lo << ".." << s_hi << "]"
       << std::endl;
    const Char_t *klass[3] = {"an", "aa", "beam"};
    const std::vector<TGraph *> *sets[3] = {&tr_an, &tr_aa, &tr_beam};
    for (Int_t ic = 0; ic < 3; ic++) {
      for (Int_t i = 0; i < Int_t(sets[ic]->size()); i++) {
        Int_t n_hi = 0;
        os << klass[ic] << " " << i;
        for (Int_t k = 0; k < (*sets[ic])[i]->GetN(); k++) {
          Double_t v = (*sets[ic])[i]->GetY()[k];
          Int_t strip = s_lo + k;
          os << " " << v;
          if (strip >= 1 && strip <= 16 && v >= 1.3)
            n_hi++;
        }
        os << " " << n_hi << std::endl;
      }
    }
  }

  DrawRegionTraces(Form("region_traces_reac%d", reac), "strip_sum_scatter",
                   tr_beam, tr_aa, tr_an, 0.6, 1.6, "#DeltaE [a.u.]");
  DrawRegionMeanTraces(Form("region_mean_traces_reac%d", reac),
                       "strip_sum_scatter", tr_beam, tr_aa, tr_an, 0.6, 1.6,
                       "#DeltaE [a.u.]");
  Double_t adc_y_lo = 0.0, adc_y_hi = 0.0;
  TraceYRange(tr_beam_adc, tr_aa_adc, tr_an_adc, adc_y_lo, adc_y_hi);
  DrawRegionTraces(Form("region_traces_reac%d_adc", reac), "strip_sum_scatter",
                   tr_beam_adc, tr_aa_adc, tr_an_adc, adc_y_lo, adc_y_hi,
                   "#DeltaE [ADC]");
  DrawRegionMeanTraces(Form("region_mean_traces_reac%d_adc", reac),
                       "strip_sum_scatter", tr_beam_adc, tr_aa_adc, tr_an_adc,
                       adc_y_lo, adc_y_hi, "#DeltaE [ADC]");
  if (!kSkipSg) {
    DrawRegionTraces(Form("region_traces_reac%d_sg", reac), "strip_sum_scatter",
                     tr_beam_sg, tr_aa_sg, tr_an_sg, 0.6, 1.6,
                     "#DeltaE [a.u.]");
    DrawRegionMeanTraces(Form("region_mean_traces_reac%d_sg", reac),
                         "strip_sum_scatter", tr_beam_sg, tr_aa_sg, tr_an_sg,
                         0.7, 1.3, "#DeltaE [a.u.]");
  }

  for (Int_t i = 0; i < tr_an.size(); i++)
    delete tr_an[i];
  for (Int_t i = 0; i < tr_aa.size(); i++)
    delete tr_aa[i];
  for (Int_t i = 0; i < tr_beam.size(); i++)
    delete tr_beam[i];
  for (Int_t i = 0; i < tr_an_adc.size(); i++)
    delete tr_an_adc[i];
  for (Int_t i = 0; i < tr_aa_adc.size(); i++)
    delete tr_aa_adc[i];
  for (Int_t i = 0; i < tr_beam_adc.size(); i++)
    delete tr_beam_adc[i];
  for (Int_t i = 0; i < tr_an_sg.size(); i++)
    delete tr_an_sg[i];
  for (Int_t i = 0; i < tr_aa_sg.size(); i++)
    delete tr_aa_sg[i];
  for (Int_t i = 0; i < tr_beam_sg.size(); i++)
    delete tr_beam_sg[i];

  delete cutAn;
  delete cutAa;

  AURestoreXErrorHandler(xerr_save);
}
