#include "StripSumScatter.hpp"

namespace {

// (a,n) reaction selection on the long-anode strips, mirroring the upstream
// dEE.C cuts. All energies are calibrated MeV normalized so the beam sits at
// Constants::NORM_MUSIC_MEV per strip.
const Int_t kReacStrip =
    Constants::STRIP_SUM_CANDIDATE_REACTION_STRIP; // upstream anode_choice
const Int_t kSmoothHiStrip =
    12; // smoothness checked through this strip (upstream k<13)
const Double_t kBeamFlatTol =
    1.0; // |E - NORM| tol for the pre-reaction strips (0..reac-1)
const Double_t kReacJumpMin =
    0.4; // reaction dE jump / excess over NORM, lower bound
const Double_t kReacJumpMax =
    5.0; // reaction dE jump / excess over NORM, upper bound
const Double_t kSmoothMaxStep =
    1.2; // max |dE| step between adjacent post-reaction strips
const Double_t kStrip17Max = 11.5; // Strip17 upper bound (upstream S17R[1])

// Max sampled per-strip traces drawn per region (beam / (a,a') / (a,n)).
const Int_t kTracesPerRegion = 40;

// Sum windows: x = all long anodes (1-16), y = the downstream reaction window
// reac+1 .. reac+6 (strips 4-9 for reac=3, matching upstream length_avDe=6;
// the published 4-10 looks like an off-by-one).
const Int_t kXLo = 1;
const Int_t kXHi = 16;
const Int_t kYLo = kReacStrip + 1;
const Int_t kYHi = kReacStrip + 6;

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
const Double_t kYMin = Constants::STRIP_SUM_YMIN;
const Double_t kYMax = Constants::STRIP_SUM_YMAX;
const Int_t kYBins = Constants::STRIP_SUM_YBINS;

BeamFit2D FindBeamGate(TChain *chain, const TString &tag,
                       const TString &subdir) {
  BeamFit2D out;
  EnergyView ev;
  ev.Attach(chain);
  TH2F *h = new TH2F(
      Form("h2_beamgate_s%d_s%d_%s", kGateStripX, kGateStripY, tag.Data()),
      Form(";#DeltaE strip %d [MeV];#DeltaE strip %d [MeV]", kGateStripX,
           kGateStripY),
      kGateBins, kGateMin, kGateMax, kGateBins, kGateMin, kGateMax);
  h->SetDirectory(nullptr);
  Long64_t n = chain->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
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

  // Save the gate plot: Strip0-vs-Strip1 with the axis-aligned selection
  // ellipse (kGateNSigmaX in x, kGateNSigmaY in y) drawn on it.
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

Bool_t PassesGate(const BeamFit2D &gate, const EnergyView &ev) {
  Double_t g0 = ev.total[kGateStripX];
  Double_t g1 = ev.total[kGateStripY];
  if (!(g0 > 0.0 && g1 > 0.0))
    return kFALSE;
  return BeamFitUtils::InEllipseXY(gate, g0, g1, kGateNSigmaX, kGateNSigmaY);
}

Bool_t AllStripsFired(const EnergyView &ev) {
  if (!(ev.total[0] > 0.0 && ev.total[17] > 0.0))
    return kFALSE;
  for (Int_t s = 1; s <= 16; s++)
    if (!(ev.total[s] > 0.0))
      return kFALSE;
  return kTRUE;
}

// Full (a,n) ladder (upstream dEE.C), assuming the beam gate already passed.
Bool_t PassesReaction(const EnergyView &ev) {
  const Double_t kNorm = Constants::NORM_MUSIC_MEV;
  if (!AllStripsFired(ev))
    return kFALSE;
  // (2) Beam still un-reacted up to the reaction strip.
  for (Int_t s = 0; s < kReacStrip; s++)
    if (TMath::Abs(ev.total[s] - kNorm) > kBeamFlatTol)
      return kFALSE;
  // (3) Energy jump at the reaction strip (relative + absolute).
  Double_t reac_jump = ev.total[kReacStrip] - ev.total[kReacStrip - 1];
  if (!(reac_jump > kReacJumpMin && reac_jump < kReacJumpMax))
    return kFALSE;
  if (!(ev.total[kReacStrip] > kNorm + kReacJumpMin &&
        ev.total[kReacStrip] < kNorm + kReacJumpMax))
    return kFALSE;
  // (4) Smooth energy loss after the reaction strip.
  for (Int_t s = kReacStrip + 1; s <= kSmoothHiStrip; s++)
    if (TMath::Abs(ev.total[s] - ev.total[s - 1]) > kSmoothMaxStep)
      return kFALSE;
  // (5) Strip17 below its upper bound.
  return ev.total[17] < kStrip17Max;
}

// Clean un-reacted beam: all long anodes flat near NORM. Beam is filtered off
// the (a,n) scatter, so its sample traces come from this separate population.
Bool_t IsBeamFlat(const EnergyView &ev) {
  const Double_t kNorm = Constants::NORM_MUSIC_MEV;
  if (!AllStripsFired(ev))
    return kFALSE;
  for (Int_t s = 1; s <= 16; s++)
    if (TMath::Abs(ev.total[s] - kNorm) > kBeamFlatTol)
      return kFALSE;
  return kTRUE;
}

Double_t SumStrips(const EnergyView &ev, Int_t lo, Int_t hi) {
  Double_t sum = 0.0;
  for (Int_t s = lo; s <= hi; s++)
    sum += ev.total[s];
  return sum;
}

// Prompt the user to draw one graphical region on the already-drawn scatter
// canvas; blocks until the polygon is closed. Returns the cut (renamed) or
// null.
TCutG *PromptCut(TCanvas *c, const char *name, const char *label) {
  std::cout << "  >>> draw the " << label
            << " region: left-click vertices, double-click to close"
            << std::endl;
  c->cd();
  // emode is case-sensitive: "CutG" puts the pad in graphical-cut mode so
  // left-clicks place vertices; the drawn TCutG is named "CUTG".
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

// Re-scan one run and collect up to kTracesPerRegion sampled per-strip traces
// per region: beam (flat, un-reacted) plus the (a,n)/(a,a') TCutGs on the
// (Sigma_1-16, Sigma_4-9) plane. Accumulates across runs into the shared lists.
void SampleRunTraces(TChain *chain, const BeamFit2D &gate, TCutG *cut_an,
                     TCutG *cut_aa, std::vector<TGraph *> &an,
                     std::vector<TGraph *> &aa, std::vector<TGraph *> &beam) {
  if (!chain || !gate.ok)
    return;
  EnergyView ev;
  ev.Attach(chain);
  Long64_t n = chain->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
    if (an.size() >= std::size_t(kTracesPerRegion) &&
        aa.size() >= std::size_t(kTracesPerRegion) &&
        beam.size() >= std::size_t(kTracesPerRegion))
      return;
    chain->GetEntry(j);
    ev.Decode();
    if (!PassesGate(gate, ev))
      continue;
    if (beam.size() < std::size_t(kTracesPerRegion) && IsBeamFlat(ev)) {
      beam.push_back(TraceCreator::BuildEventTrace(ev));
      continue;
    }
    if (!PassesReaction(ev))
      continue;
    Double_t x = SumStrips(ev, kXLo, kXHi);
    Double_t y = SumStrips(ev, kYLo, kYHi);
    if (cut_an && an.size() < std::size_t(kTracesPerRegion) &&
        cut_an->IsInside(x, y))
      an.push_back(TraceCreator::BuildEventTrace(ev));
    else if (cut_aa && aa.size() < std::size_t(kTracesPerRegion) &&
             cut_aa->IsInside(x, y))
      aa.push_back(TraceCreator::BuildEventTrace(ev));
  }
}

void DrawTraceSet(const std::vector<TGraph *> &traces, Int_t color) {
  for (std::size_t i = 0; i < traces.size(); i++) {
    traces[i]->SetLineColorAlpha(color, 0.45);
    traces[i]->SetLineWidth(1);
    traces[i]->Draw("L SAME");
  }
}

// Overlay the sampled per-strip traces, one colour per region, on one canvas.
void DrawRegionTraces(const std::vector<TGraph *> &beam,
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

  PlottingUtils::SaveFigure(c, Form("region_traces_reac%d", kReacStrip),
                            "strip_sum_scatter", PlotSaveOptions::kLINEAR);
  delete leg;
  delete p_beam;
  delete p_aa;
  delete p_an;
  delete c;
  delete frame;
}

} // namespace

BeamFit2D StripSumScatter::FillRun(Int_t run, TChain *chain,
                                   const TString &gate_subdir, TH2F *h) {
  BeamFit2D gate;
  if (!chain || chain->GetEntries() == 0) {
    std::cerr << "Run " << run << ": empty chain; skipping" << std::endl;
    return gate;
  }

  gate = FindBeamGate(chain, Form("run%d", run), gate_subdir);
  if (!gate.ok) {
    std::cerr << "Run " << run << ": beam gate (strips " << kGateStripX << ","
              << kGateStripY << ") failed; skipping" << std::endl;
    return gate;
  }
  std::cout << "  beam gate (strips " << kGateStripX << "," << kGateStripY
            << "): mu=(" << gate.mu_x << "," << gate.mu_y << ") sigma=("
            << gate.sigma_x << "," << gate.sigma_y << ")" << std::endl;

  EnergyView ev;
  ev.Attach(chain);

  Long64_t n = chain->GetEntries();
  std::cout << "Run " << run << ": (a,n)-gating " << n
            << " events (beam ellipse on strips " << kGateStripX << ","
            << kGateStripY << ", then reaction-jump + smoothness on strip "
            << kReacStrip << ")..." << std::endl;

  Long64_t n_gated = 0;
  for (Long64_t j = 0; j < n; j++) {
    chain->GetEntry(j);
    ev.Decode();
    if (!PassesGate(gate, ev))
      continue;
    if (!PassesReaction(ev))
      continue;
    n_gated++;
    h->Fill(SumStrips(ev, kXLo, kXHi), SumStrips(ev, kYLo, kYHi));
  }
  std::cout << "  " << n_gated << " events passed the PID gate" << std::endl;
  return gate;
}

void StripSumScatter::Run() {
  const TString project_root = Paths::DatasetDir();
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");
  gROOT->SetBatch(kTRUE);

  // Group cal sidecars by run; each run becomes one TChain over events_cal.
  std::vector<Int_t> run_order;
  std::map<Int_t, TChain *> chain_by_run =
      FileSet::GroupCalSidecarsByRun(run_order);

  // One aggregate histogram over all runs. Each run is gated with its own beam
  // ellipse (the per-run calibration leaves residual epoch-to-epoch shifts),
  // but the selected events all accumulate here for a single combined scatter.
  TH2F *h =
      new TH2F(Form("h2_normsumE_s%d_%d_vs_s%d_%d", kYLo, kYHi, kXLo, kXHi),
               Form(";norm. #DeltaE strips %d#rightarrow%d [MeV];"
                    "norm. #DeltaE strips %d#rightarrow%d [MeV]",
                    kXLo, kXHi, kYLo, kYHi),
               kXBins, kXMin, kXMax, kYBins, kYMin, kYMax);
  h->SetDirectory(nullptr);
  h->SetStats(0);

  std::map<Int_t, BeamFit2D> gates;
  for (std::size_t i = 0; i < run_order.size(); i++) {
    Int_t run = run_order[i];
    gates[run] = FillRun(run, chain_by_run[run],
                         Form("strip_sum_scatter/run%d", run), h);
  }

  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h, c);
    h->GetYaxis()->SetTitleOffset(1.3);
    c->SetLeftMargin(0.18);
    PlottingUtils::SaveFigure(c,
                              Form("normsumE_reac%d_s%d_%d_vs_s%d_%d",
                                   kReacStrip, kYLo, kYHi, kXLo, kXHi),
                              "strip_sum_scatter", PlotSaveOptions::kLINEAR);
    delete c;
  }

  // Interactive region-trace overlay: draw the aggregate, let the user draw the
  // (a,n) and (a,a') regions, then sample per-strip traces from each region
  // (plus a separate flat-beam population) and overlay them on one plot. Needs
  // a display; headless runs just keep the aggregate scatter already saved
  // above.
  if (!gSystem->Getenv("DISPLAY")) {
    std::cerr << "strip-sum-scatter: no DISPLAY; skipping interactive "
                 "region-trace overlay (aggregate scatter already saved)."
              << std::endl;
  } else {
    Int_t app_argc = 1;
    char app_arg0[] = "strip-sum-scatter";
    char *app_argv[] = {app_arg0};
    TApplication app("strip-sum-scatter", &app_argc, app_argv);
    gROOT->SetBatch(kFALSE);

    TCanvas *cut_canvas = new TCanvas(
        "c_strip_sum_regions", "Draw (a,n) then (a,a') regions", 900, 700);
    h->Draw("COLZ");
    cut_canvas->Update();
    TCutG *cut_an = PromptCut(cut_canvas, "region_an", "(a,n)");
    TCutG *cut_aa = PromptCut(cut_canvas, "region_aa", "(a,a')");

    std::vector<TGraph *> tr_an, tr_aa, tr_beam;
    for (std::size_t i = 0; i < run_order.size(); i++)
      SampleRunTraces(chain_by_run[run_order[i]], gates[run_order[i]], cut_an,
                      cut_aa, tr_an, tr_aa, tr_beam);

    std::cout << "Sampled traces: beam=" << tr_beam.size()
              << " (a,a')=" << tr_aa.size() << " (a,n)=" << tr_an.size()
              << std::endl;
    DrawRegionTraces(tr_beam, tr_aa, tr_an);

    for (std::size_t i = 0; i < tr_an.size(); i++)
      delete tr_an[i];
    for (std::size_t i = 0; i < tr_aa.size(); i++)
      delete tr_aa[i];
    for (std::size_t i = 0; i < tr_beam.size(); i++)
      delete tr_beam[i];
  }

  delete h;

  for (std::size_t i = 0; i < run_order.size(); i++)
    delete chain_by_run[run_order[i]];
}
