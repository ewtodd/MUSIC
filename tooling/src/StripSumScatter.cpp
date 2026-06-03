#include "StripSumScatter.hpp"
#include <TNamed.h>
#include <limits>

namespace {

// (a,n) reaction selection on the long-anode strips, mirroring the upstream
// dEE.C cuts. All energies are calibrated MeV normalized so the beam sits at
// Constants::NORM_MUSIC_MEV per strip.
const Int_t kSmoothHiStrip =
    12; // smoothness checked through this (absolute) strip (upstream k<13)
const Double_t kBeamFlatTol =
    1.0; // |E - NORM| tol for the pre-reaction strips (0..reac-1)
const Double_t kReacJumpMin =
    0.4; // reaction dE jump / excess over NORM, lower bound
const Double_t kReacJumpMax =
    5.0; // reaction dE jump / excess over NORM, upper bound
const Double_t kSmoothMaxStep =
    1.2; // max |dE| step between adjacent post-reaction strips
const Double_t kStrip17Max = 11.5; // Strip17 upper bound (upstream S17R[1])

// Candidate reaction strips: every scatter is computed in one pass over the
// data, one TH2F per reaction strip in [MIN, MAX].
const Int_t kReacMin = Constants::REACTION_STRIP_MIN;
const Int_t kReacMax = Constants::REACTION_STRIP_MAX;
const Int_t kNReac = kReacMax - kReacMin + 1;

// Max sampled per-strip traces drawn per region (beam / (a,a') / (a,n)).
const Int_t kTracesPerRegion = 40;
// Cap on beam-flat events kept in the trace reservoir (only ~40 are drawn).
const Int_t kBeamReservoirCap = 400;

// Sum windows: x = all long anodes (1-16); y = the downstream reaction window
// reac+1 .. reac+6, clamped to the last strip (17). The window shortens for
// reaction strips near the downstream end.
const Int_t kXLo = 1;
const Int_t kXHi = 16;

// Beam gate (kept): Strip0 vs Strip1 inside the calibration moments ellipse,
// axis-aligned at (kGateNSigmaX, kGateNSigmaY) exactly like the calibration.
const Int_t kGateStripX = 0; // Strip0 (gate x)
const Int_t kGateStripY = 1; // Strip1 (gate y)
const Double_t kGateNSigmaX = 2.0;
const Double_t kGateNSigmaY = 2.0;

const Double_t kGateMin = 0.0;
const Double_t kGateMax = 3.0 * Constants::NORM_MUSIC_MEV;
const Int_t kGateBins = 240;
const Int_t kSeedHalfBins = 40;
const Double_t kSeedFrac = 0.30;

const Double_t kXMin = Constants::STRIP_SUM_XMIN;
const Double_t kXMax = Constants::STRIP_SUM_XMAX;
const Int_t kXBins = Constants::STRIP_SUM_XBINS;
const Int_t kYBins = Constants::STRIP_SUM_YBINS;

// Sampled-pass budget for the beam-gate fit and the y-bound auto-range. Both
// only need to characterize a population, not visit every event.
const Long64_t kSampleMaxPoints = 2000000;

// Cache file (under root_files): all per-strip scatters + the trace reservoir,
// stamped with a fingerprint of the cut constants and input entry counts.
const char *kCacheName = "StripSumScatter_cache.root";

Int_t ReacIndex(Int_t reac) { return reac - kReacMin; }
Int_t YLoOf(Int_t reac) { return reac + 1; }
Int_t YHiOf(Int_t reac) { return TMath::Min(reac + 6, 17); }

// One sampled trace per kept event: the 18 strip totals (enough to redraw the
// per-strip trace and recompute cut membership), the bitmask of which reaction
// strips it passes, and whether it is clean flat beam.
struct TraceEvt {
  Float_t total[18];
  UInt_t reac_mask;
  Bool_t beam_flat;
};

void EnableEventBranches(TChain *chain) {
  // Skip decompressing the branches EnergyView does not read
  // (Hits/Grid/FlagsOR) on these full/sampled scans.
  chain->SetBranchStatus("*", 0);
  chain->SetBranchStatus("Left_0_17_dE", 1);
  chain->SetBranchStatus("RightdE", 1);
  chain->SetBranchStatus("Cathode", 1);
}

Bool_t AllStripsFired(const EnergyView &ev) {
  if (!(ev.total[0] > 0.0 && ev.total[17] > 0.0))
    return kFALSE;
  for (Int_t s = 1; s <= 16; s++)
    if (!(ev.total[s] > 0.0))
      return kFALSE;
  return kTRUE;
}

// Full (a,n) ladder (upstream dEE.C) for a given reaction strip, assuming the
// beam gate already passed. The smoothness check runs from reac+1 up to the
// absolute strip kSmoothHiStrip, so it is vacuous once reac >= kSmoothHiStrip
// (the high reaction strips therefore carry no smoothness constraint -- a
// physics knob to revisit if those strips matter).
Bool_t PassesReaction(const EnergyView &ev, Int_t reac) {
  const Double_t kNorm = Constants::NORM_MUSIC_MEV;
  if (!AllStripsFired(ev))
    return kFALSE;
  for (Int_t s = 0; s < reac; s++)
    if (TMath::Abs(ev.total[s] - kNorm) > kBeamFlatTol)
      return kFALSE;
  Double_t reac_jump = ev.total[reac] - ev.total[reac - 1];
  if (!(reac_jump > kReacJumpMin && reac_jump < kReacJumpMax))
    return kFALSE;
  if (!(ev.total[reac] > kNorm + kReacJumpMin &&
        ev.total[reac] < kNorm + kReacJumpMax))
    return kFALSE;
  for (Int_t s = reac + 1; s <= kSmoothHiStrip; s++)
    if (TMath::Abs(ev.total[s] - ev.total[s - 1]) > kSmoothMaxStep)
      return kFALSE;
  return ev.total[17] < kStrip17Max;
}

// Clean un-reacted beam: all long anodes flat near NORM.
Bool_t IsBeamFlat(const EnergyView &ev) {
  const Double_t kNorm = Constants::NORM_MUSIC_MEV;
  if (!AllStripsFired(ev))
    return kFALSE;
  for (Int_t s = 1; s <= 16; s++)
    if (TMath::Abs(ev.total[s] - kNorm) > kBeamFlatTol)
      return kFALSE;
  return kTRUE;
}

Double_t SumRange(const Double_t *total, Int_t lo, Int_t hi) {
  Double_t sum = 0.0;
  for (Int_t s = lo; s <= hi; s++)
    sum += total[s];
  return sum;
}

Bool_t PassesGate(const BeamFit2D &gate, const EnergyView &ev) {
  Double_t g0 = ev.total[kGateStripX];
  Double_t g1 = ev.total[kGateStripY];
  if (!(g0 > 0.0 && g1 > 0.0))
    return kFALSE;
  return BeamFitUtils::InEllipseXY(gate, g0, g1, kGateNSigmaX, kGateNSigmaY);
}

// Fit one run's Strip0-vs-Strip1 beam gate from a stride-sampled subset (the
// 2D peak is characterized fine without every event), saving the gate plot.
BeamFit2D FindBeamGate(TChain *chain, const TString &tag,
                       const TString &subdir) {
  BeamFit2D out;
  EnergyView ev;
  ev.Attach(chain);
  EnableEventBranches(chain);
  TH2F *h = new TH2F(
      Form("h2_beamgate_s%d_s%d_%s", kGateStripX, kGateStripY, tag.Data()),
      Form(";#DeltaE strip %d [MeV];#DeltaE strip %d [MeV]", kGateStripX,
           kGateStripY),
      kGateBins, kGateMin, kGateMax, kGateBins, kGateMin, kGateMax);
  h->SetDirectory(nullptr);
  Long64_t n = chain->GetEntries();
  Long64_t stride = FileSet::SampleStride(n, kSampleMaxPoints);
  for (Long64_t j = 0; j < n; j += stride) {
    chain->GetEntry(j);
    ev.Decode();
    Double_t x = ev.total[kGateStripX];
    Double_t y = ev.total[kGateStripY];
    if (x > 0.0 && y > 0.0)
      h->Fill(x, y);
  }
  if (h->GetEntries() < 100) {
    delete h;
    return out;
  }
  Double_t bw_x = h->GetXaxis()->GetBinWidth(1);
  Double_t bw_y = h->GetYaxis()->GetBinWidth(1);
  Int_t bx = 0, by = 0, bz = 0;
  h->GetMaximumBin(bx, by, bz);
  Double_t peak_val = h->GetBinContent(bx, by);
  Int_t lo_bx = std::max(1, bx - kSeedHalfBins);
  Int_t hi_bx = std::min(h->GetNbinsX(), bx + kSeedHalfBins);
  Int_t lo_by = std::max(1, by - kSeedHalfBins);
  Int_t hi_by = std::min(h->GetNbinsY(), by + kSeedHalfBins);
  Moments2D m = BeamFitUtils::ComputeMoments(h, lo_bx, hi_bx, lo_by, hi_by,
                                             kSeedFrac * peak_val, bw_x, bw_y);
  if (m.weight <= 0) {
    delete h;
    return out;
  }
  out.amp = peak_val;
  out.mu_x = m.mu_x;
  out.mu_y = m.mu_y;
  out.sigma_x = m.sigma_x;
  out.sigma_y = m.sigma_y;
  out.rho = m.rho;
  out.ok = kTRUE;

  if (Constants::SAVE_PLOTS) {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h, c);
    TEllipse *e = new TEllipse(out.mu_x, out.mu_y, kGateNSigmaX * out.sigma_x,
                               kGateNSigmaY * out.sigma_y);
    e->SetFillStyle(0);
    e->SetLineColor(kRed + 1);
    e->SetLineWidth(2);
    e->Draw();
    PlottingUtils::SaveFigure(
        c, Form("beam_gate_s%d_s%d", kGateStripX, kGateStripY), subdir,
        PlotSaveOptions::kLINEAR);
    delete c;
  }
  delete h;
  return out;
}

void DrawTraceSet(const std::vector<TGraph *> &traces, Int_t color) {
  for (std::size_t i = 0; i < traces.size(); i++) {
    traces[i]->SetLineColorAlpha(color, 0.45);
    traces[i]->SetLineWidth(1);
    traces[i]->Draw("L SAME");
  }
}

TGraph *TraceFromTotal(const Float_t *total) {
  Double_t td[18];
  for (Int_t s = 0; s < 18; s++)
    td[s] = Double_t(total[s]);
  return TraceCreator::BuildTraceFromTotals(td);
}

// Overlay the sampled per-strip traces, one colour per region, on one canvas.
void DrawRegionTraces(Int_t reac, const std::vector<TGraph *> &beam,
                      const std::vector<TGraph *> &aa,
                      const std::vector<TGraph *> &an) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);
  TH2F *frame =
      new TH2F("h_region_trace_frame", ";Strip;#DeltaE [a.u.]", 18, -0.5, 17.5,
               100, Constants::STRIP_E_MIN_MEV, Constants::STRIP_E_MAX_MEV);
  frame->SetStats(0);
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  frame->Draw();
  DrawTraceSet(beam, kGray + 2);
  DrawTraceSet(aa, kAzure + 2);
  DrawTraceSet(an, kRed + 1);

  TGraph *p_beam = new TGraph(1);
  TGraph *p_aa = new TGraph(1);
  TGraph *p_an = new TGraph(1);
  TGraph *proxies[3] = {p_beam, p_aa, p_an};
  Int_t pcol[3] = {kGray + 2, kAzure + 2, kRed + 1};
  for (Int_t i = 0; i < 3; i++) {
    proxies[i]->SetPoint(0, -1e9, -1e9);
    proxies[i]->SetLineColor(pcol[i]);
    proxies[i]->SetLineWidth(3);
  }
  TLegend *leg = PlottingUtils::AddLegend(0.725, 0.875, 0.70, 0.86);
  leg->AddEntry(p_beam, "Beam", "l");
  leg->AddEntry(p_aa, "(#alpha,#alpha')", "l");
  leg->AddEntry(p_an, "(#alpha,n)", "l");
  leg->Draw();

  PlottingUtils::SaveFigure(c, Form("region_traces_reac%d", reac),
                            "strip_sum_scatter", PlotSaveOptions::kLINEAR);
  delete leg;
  delete p_beam;
  delete p_aa;
  delete p_an;
  delete c;
  delete frame;
}

// Prompt the user to draw one graphical region on the already-drawn scatter
// canvas; blocks until the polygon is closed. Returns the cut (renamed) or
// null.
TCutG *PromptCut(TCanvas *c, const char *name, const char *label) {
  std::cout << "  >>> draw the " << label
            << " region: left-click vertices, double-click to close"
            << std::endl;
  c->cd();
  TCutG *cut = dynamic_cast<TCutG *>(c->WaitPrimitive("CUTG", "CutG"));
  if (!cut) {
    std::cerr << "  no " << label << " cut drawn" << std::endl;
    return nullptr;
  }
  cut->SetName(name);
  cut->SetLineColor(kBlack);
  cut->SetLineWidth(2);
  return cut;
}

TString BuildFingerprint(const std::vector<Int_t> &run_order,
                         std::map<Int_t, TChain *> &chains) {
  TString s = Form(
      "v1 reac[%d,%d] norm=%.4f flatTol=%.3f jump[%.3f,%.3f] smoothHi=%d "
      "step=%.3f s17=%.3f gate[s%d,s%d,%.2f,%.2f,%d,%.3f,%.3f] x[%.3f,%.3f,%d] "
      "ybins=%d",
      kReacMin, kReacMax, Constants::NORM_MUSIC_MEV, kBeamFlatTol, kReacJumpMin,
      kReacJumpMax, kSmoothHiStrip, kSmoothMaxStep, kStrip17Max, kGateStripX,
      kGateStripY, kGateNSigmaX, kGateNSigmaY, kGateBins, kGateMin, kGateMax,
      kXMin, kXMax, kXBins, kYBins);
  for (std::size_t i = 0; i < run_order.size(); i++) {
    Int_t run = run_order[i];
    s += Form(" r%d:%lld", run, chains[run]->GetEntries());
  }
  return s;
}

// Per-reaction-strip y-axis bounds from a stride-sampled pass over all runs
// (x is strip-independent and stays fixed). Falls back to the configured
// STRIP_SUM_Y range for strips with too few sampled events.
void FindYBounds(const std::vector<Int_t> &run_order,
                 std::map<Int_t, TChain *> &chains,
                 std::map<Int_t, BeamFit2D> &gates, Double_t *y_lo,
                 Double_t *y_hi) {
  for (Int_t r = 0; r < kNReac; r++) {
    y_lo[r] = std::numeric_limits<Double_t>::max();
    y_hi[r] = -std::numeric_limits<Double_t>::max();
  }
  for (std::size_t i = 0; i < run_order.size(); i++) {
    Int_t run = run_order[i];
    TChain *chain = chains[run];
    if (!chain || !gates[run].ok)
      continue;
    EnergyView ev;
    ev.Attach(chain);
    EnableEventBranches(chain);
    Long64_t n = chain->GetEntries();
    Long64_t stride = FileSet::SampleStride(n, kSampleMaxPoints);
    for (Long64_t j = 0; j < n; j += stride) {
      chain->GetEntry(j);
      ev.Decode();
      if (!PassesGate(gates[run], ev))
        continue;
      for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
        if (!PassesReaction(ev, reac))
          continue;
        Double_t y = SumRange(ev.total, YLoOf(reac), YHiOf(reac));
        Int_t ri = ReacIndex(reac);
        if (y < y_lo[ri])
          y_lo[ri] = y;
        if (y > y_hi[ri])
          y_hi[ri] = y;
      }
    }
  }
  for (Int_t r = 0; r < kNReac; r++) {
    if (y_lo[r] > y_hi[r]) { // no events sampled for this strip
      y_lo[r] = Constants::STRIP_SUM_YMIN;
      y_hi[r] = Constants::STRIP_SUM_YMAX;
      continue;
    }
    Double_t pad = 0.05 * (y_hi[r] - y_lo[r]);
    if (pad < 1.0)
      pad = 1.0;
    y_lo[r] -= pad;
    y_hi[r] += pad;
    if (y_lo[r] < 0.0)
      y_lo[r] = 0.0;
  }
}

} // namespace

void StripSumScatter::Run() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                Paths::ResultsDir() + "/plots",
                                Paths::ResultsDir() + "/root_files");
  gROOT->SetBatch(kTRUE);

  std::vector<Int_t> run_order;
  std::map<Int_t, TChain *> chain_by_run = FileSet::GroupEventsByRun(run_order);
  if (run_order.empty()) {
    std::cerr << "strip-sum-scatter: no runs found" << std::endl;
    return;
  }

  TString fingerprint = BuildFingerprint(run_order, chain_by_run);

  std::map<Int_t, TH2F *> scatter; // keyed by reaction strip
  std::vector<TraceEvt> reservoir;

  // ---- try the cache ----
  Bool_t loaded = kFALSE;
  TString cache_full = IO::GetRootFilesBaseDir() + TString("/") + kCacheName;
  if (!gSystem->AccessPathName(cache_full)) {
    TFile *cf = IO::OpenForReading(kCacheName);
    if (cf && !cf->IsZombie()) {
      TNamed *fp = dynamic_cast<TNamed *>(cf->Get("fingerprint"));
      if (fp && fingerprint == fp->GetTitle()) {
        Bool_t ok = kTRUE;
        for (Int_t reac = kReacMin; reac <= kReacMax && ok; reac++) {
          TH2F *h = dynamic_cast<TH2F *>(cf->Get(Form("scatter_r%d", reac)));
          if (!h) {
            ok = kFALSE;
            break;
          }
          TH2F *hc = static_cast<TH2F *>(h->Clone());
          hc->SetDirectory(nullptr);
          scatter[reac] = hc;
        }
        TTree *tt = dynamic_cast<TTree *>(cf->Get("traces"));
        if (ok && tt) {
          TraceEvt e;
          tt->SetBranchAddress("total", e.total);
          tt->SetBranchAddress("reac_mask", &e.reac_mask);
          tt->SetBranchAddress("beam_flat", &e.beam_flat);
          Long64_t nt = tt->GetEntries();
          reservoir.reserve(nt);
          for (Long64_t j = 0; j < nt; j++) {
            tt->GetEntry(j);
            reservoir.push_back(e);
          }
          loaded = ok;
        }
      }
      cf->Close();
      delete cf;
    }
    if (loaded)
      std::cout << "strip-sum-scatter: loaded cached scatters + "
                << reservoir.size() << " reservoir events (fingerprint match)."
                << std::endl;
    else
      std::cout << "strip-sum-scatter: cache present but stale; rebuilding."
                << std::endl;
  }

  // ---- build (one sampled gate pass + one sampled bounds pass + one full
  // fill pass over all runs, producing every reaction-strip scatter) ----
  if (!loaded) {
    std::map<Int_t, BeamFit2D> gates;
    for (std::size_t i = 0; i < run_order.size(); i++) {
      Int_t run = run_order[i];
      if (!chain_by_run[run] || chain_by_run[run]->GetEntries() == 0)
        continue;
      gates[run] = FindBeamGate(chain_by_run[run], Form("run%d", run),
                                Form("strip_sum_scatter/run%d", run));
      if (gates[run].ok)
        std::cout << "  run " << run << " beam gate: mu=(" << gates[run].mu_x
                  << "," << gates[run].mu_y << ")" << std::endl;
      else
        std::cerr << "  run " << run << " beam gate failed; skipping"
                  << std::endl;
    }

    Double_t y_lo[64], y_hi[64];
    FindYBounds(run_order, chain_by_run, gates, y_lo, y_hi);

    for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
      Int_t ri = ReacIndex(reac);
      TH2F *h = new TH2F(
          Form("scatter_r%d", reac),
          Form(";norm. #DeltaE strips %d#rightarrow%d [MeV];norm. #DeltaE "
               "strips %d#rightarrow%d [MeV]",
               kXLo, kXHi, YLoOf(reac), YHiOf(reac)),
          kXBins, kXMin, kXMax, kYBins, y_lo[ri], y_hi[ri]);
      h->SetDirectory(nullptr);
      h->SetStats(0);
      scatter[reac] = h;
    }

    Long64_t total_gated = 0, total_seen = 0;
    Int_t n_beam_kept = 0;
    for (std::size_t i = 0; i < run_order.size(); i++) {
      Int_t run = run_order[i];
      TChain *chain = chain_by_run[run];
      if (!chain || !gates[run].ok)
        continue;
      EnergyView ev;
      ev.Attach(chain);
      EnableEventBranches(chain);
      Long64_t n = chain->GetEntries();
      std::cout << "Run " << run << ": filling " << kNReac
                << " reaction-strip scatters over " << n << " events..."
                << std::endl;
      for (Long64_t j = 0; j < n; j++) {
        chain->GetEntry(j);
        ev.Decode();
        total_seen++;
        if (!PassesGate(gates[run], ev))
          continue;
        Double_t x = SumRange(ev.total, kXLo, kXHi);
        UInt_t mask = 0;
        for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
          if (!PassesReaction(ev, reac))
            continue;
          mask |= (1u << ReacIndex(reac));
          scatter[reac]->Fill(x, SumRange(ev.total, YLoOf(reac), YHiOf(reac)));
        }
        // Keep every reaction-passing event for traces; cap clean flat-beam
        // events (only ~kTracesPerRegion are ever drawn). The two are mutually
        // exclusive -- a flat-beam event has no reaction jump.
        Bool_t beam = (mask == 0) && IsBeamFlat(ev);
        if (mask == 0 && !(beam && n_beam_kept < kBeamReservoirCap))
          continue;
        if (beam)
          n_beam_kept++;
        TraceEvt e;
        for (Int_t s = 0; s < 18; s++)
          e.total[s] = Float_t(ev.total[s]);
        e.reac_mask = mask;
        e.beam_flat = beam;
        reservoir.push_back(e);
        if (mask != 0)
          total_gated++;
      }
    }
    std::cout << "Built scatters: " << total_gated
              << " reaction events across strips " << kReacMin << "-"
              << kReacMax << " (" << total_seen << " seen), reservoir "
              << reservoir.size() << " events." << std::endl;

    // ---- write cache ----
    TFile *out = IO::OpenForWriting(kCacheName, "RECREATE");
    if (out && !out->IsZombie()) {
      out->cd();
      TNamed fp("fingerprint", fingerprint.Data());
      fp.Write();
      for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
        scatter[reac]->Write(Form("scatter_r%d", reac));
      TTree *tt = new TTree("traces", "strip-sum trace reservoir");
      TraceEvt e;
      tt->Branch("total", e.total, "total[18]/F");
      tt->Branch("reac_mask", &e.reac_mask, "reac_mask/i");
      tt->Branch("beam_flat", &e.beam_flat, "beam_flat/O");
      for (std::size_t k = 0; k < reservoir.size(); k++) {
        e = reservoir[k];
        tt->Fill();
      }
      tt->Write();
      out->Close();
      delete out;
      std::cout << "strip-sum-scatter: wrote cache " << kCacheName << std::endl;
    }
  }

  // ---- save every reaction-strip scatter PNG ----
  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureAndDraw2DHistogram(scatter[reac], c);
      scatter[reac]->GetYaxis()->SetTitleOffset(1.3);
      c->SetLeftMargin(0.18);
      PlottingUtils::SaveFigure(c,
                                Form("normsumE_reac%d_s%d_%d_vs_s%d_%d", reac,
                                     YLoOf(reac), YHiOf(reac), kXLo, kXHi),
                                "strip_sum_scatter", PlotSaveOptions::kLINEAR);
      delete c;
    }
  }

  // ---- interactive region-trace overlay on the configured reaction strip ----
  // Needs a display; samples traces from the cached reservoir (no rescan).
  Int_t reac = Constants::STRIP_SUM_CANDIDATE_REACTION_STRIP;
  if (reac < kReacMin || reac > kReacMax) {
    std::cerr << "strip-sum-scatter: candidate reaction strip " << reac
              << " outside [" << kReacMin << "," << kReacMax
              << "]; skipping interactive overlay." << std::endl;
  } else if (!gSystem->Getenv("DISPLAY")) {
    std::cerr << "strip-sum-scatter: no DISPLAY; skipping interactive "
                 "region-trace overlay (scatters already saved)."
              << std::endl;
  } else {
    Int_t app_argc = 1;
    char app_arg0[] = "strip-sum-scatter";
    char *app_argv[] = {app_arg0};
    TApplication app("strip-sum-scatter", &app_argc, app_argv);
    gROOT->SetBatch(kFALSE);

    TCanvas *cut_canvas = new TCanvas(
        "c_strip_sum_regions", "Draw (a,n) then (a,a') regions", 900, 700);
    cut_canvas->SetLogz(kTRUE); // match the saved scatter's z-scale
    scatter[reac]->Draw("COLZ");
    cut_canvas->Update();
    TCutG *cut_an = PromptCut(cut_canvas, "region_an", "(a,n)");
    TCutG *cut_aa = PromptCut(cut_canvas, "region_aa", "(a,a')");

    std::vector<TGraph *> tr_an, tr_aa, tr_beam;
    UInt_t bit = (1u << ReacIndex(reac));
    for (std::size_t k = 0; k < reservoir.size(); k++) {
      if (tr_an.size() >= std::size_t(kTracesPerRegion) &&
          tr_aa.size() >= std::size_t(kTracesPerRegion) &&
          tr_beam.size() >= std::size_t(kTracesPerRegion))
        break;
      const TraceEvt &e = reservoir[k];
      if (e.beam_flat && tr_beam.size() < std::size_t(kTracesPerRegion)) {
        tr_beam.push_back(TraceFromTotal(e.total));
        continue;
      }
      if (!(e.reac_mask & bit))
        continue;
      Double_t td[18];
      for (Int_t s = 0; s < 18; s++)
        td[s] = Double_t(e.total[s]);
      Double_t x = SumRange(td, kXLo, kXHi);
      Double_t y = SumRange(td, YLoOf(reac), YHiOf(reac));
      if (cut_an && tr_an.size() < std::size_t(kTracesPerRegion) &&
          cut_an->IsInside(x, y))
        tr_an.push_back(TraceFromTotal(e.total));
      else if (cut_aa && tr_aa.size() < std::size_t(kTracesPerRegion) &&
               cut_aa->IsInside(x, y))
        tr_aa.push_back(TraceFromTotal(e.total));
    }

    std::cout << "Sampled traces: beam=" << tr_beam.size()
              << " (a,a')=" << tr_aa.size() << " (a,n)=" << tr_an.size()
              << std::endl;
    DrawRegionTraces(reac, tr_beam, tr_aa, tr_an);

    for (std::size_t i = 0; i < tr_an.size(); i++)
      delete tr_an[i];
    for (std::size_t i = 0; i < tr_aa.size(); i++)
      delete tr_aa[i];
    for (std::size_t i = 0; i < tr_beam.size(); i++)
      delete tr_beam[i];
  }

  for (Int_t reac2 = kReacMin; reac2 <= kReacMax; reac2++)
    delete scatter[reac2];
  for (std::size_t i = 0; i < run_order.size(); i++)
    delete chain_by_run[run_order[i]];
}
