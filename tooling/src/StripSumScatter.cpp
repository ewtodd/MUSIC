#include "StripSumScatter.hpp"

namespace {

// "4-10 vs 1-16" sum. The standard ("normMUSIC") calibration already scales
// each channel's beam to the common reference NORM_MUSIC_MEV, so the calibrated
// per-strip energies are already normalized -- we just sum them directly. A
// clean beam event sits at ~NORM_MUSIC_MEV per strip, so Sum 1-16 ~ 16*ref and
// Sum 4-10 ~ 7*ref for beam; reactions deviate from that ridge. Long side only
// (short sides are 0 in MeV; requires IGNORE_SHORT_STRIPS). y-axis = sum over
// [kYLo, kYHi]; x-axis = sum over [kXLo, kXHi].
const Int_t kXLo = 1;
const Int_t kXHi = 16;
const Int_t kYLo = 4;
const Int_t kYHi = 10;

// PID gate (hardcoded): beam selected by the calibration gate -- the
// Strip0-vs-Strip1 beam ellipse (cleaner than strips 1&2, which include the
// noisy R-even side), axis-aligned at (kGateNSigmaX, kGateNSigmaY) exactly like
// the calibration -- AND more dE in strip 3 than in strip 2.
const Int_t kGateStripX = 0;       // Strip0 (gate x)
const Int_t kGateStripY = 1;       // Strip1 (gate y)
const Int_t kReacStrip = 3;        // reaction strip
const Int_t kPrevStrip = 2;        // require dE(kReacStrip) > dE(kPrevStrip)
const Double_t kGateNSigmaX = 1.0; // matches calibration kEllipseNSigmaX
const Double_t kGateNSigmaY = 3.0; // matches calibration kEllipseNSigmaY
// Gate histogram spans 0..3x the beam reference (covers beam + single pileup).
const Double_t kGateMin = 0.0;
const Double_t kGateMax = 3.0 * Constants::NORM_MUSIC_MEV;
const Int_t kGateBins = 240;
const Int_t kSeedHalfBins = 40;
const Double_t kSeedFrac = 0.30;

const Double_t kXMin = 159.5;
const Double_t kXMax = 264.5;
const Int_t kXBins = 300;
const Double_t kYMin = 63.8;
const Double_t kYMax = 115.7;
const Int_t kYBins = 300;

// The calibration beam gate: the Strip0-vs-Strip1 dE peak, located via local
// moments around the largest bin (identical to CalibrateBeam's gate). Events
// are later gated on this ellipse with InEllipseXY at (kGateNSigmaX,Y).
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

} // namespace

void StripSumScatter::ProcessRun(Int_t run, TChain *chain,
                                 const TString &subdir) {
  if (!chain || chain->GetEntries() == 0) {
    std::cerr << "Run " << run << ": empty chain; skipping" << std::endl;
    return;
  }

  BeamFit2D gate = FindBeamGate(chain, Form("run%d", run), subdir);
  if (!gate.ok) {
    std::cerr << "Run " << run << ": beam gate (strips " << kGateStripX << ","
              << kGateStripY << ") failed; skipping" << std::endl;
    return;
  }
  std::cout << "  beam gate (strips " << kGateStripX << "," << kGateStripY
            << "): mu=(" << gate.mu_x << "," << gate.mu_y << ") sigma=("
            << gate.sigma_x << "," << gate.sigma_y << ")" << std::endl;

  // EnergyView reads MeV/ADC transparently and collapses each strip to its
  // long side when IGNORE_SHORT_STRIPS is set (in MeV the short side is
  // already 0 since it is uncalibrated).
  EnergyView ev;
  ev.Attach(chain);

  TH2F *h = new TH2F(
      Form("h2_normsumE_s%d_%d_vs_s%d_%d_run%d", kYLo, kYHi, kXLo, kXHi, run),
      Form(";norm. #DeltaE strips %d#rightarrow%d [MeV];"
           "norm. #DeltaE strips %d#rightarrow%d [MeV]",
           kXLo, kXHi, kYLo, kYHi),
      kXBins, kXMin, kXMax, kYBins, kYMin, kYMax);
  h->SetDirectory(nullptr);
  h->SetStats(0);

  Long64_t n = chain->GetEntries();
  std::cout << "Run " << run << ": PID-gating " << n
            << " events (calibration beam gate on strips " << kGateStripX << ","
            << kGateStripY << ", then strip " << kReacStrip << " > strip "
            << kPrevStrip << ")..." << std::endl;
  Long64_t n_gated = 0;
  for (Long64_t j = 0; j < n; j++) {
    chain->GetEntry(j);
    ev.Decode();

    // Condition 1: clean 87Rb beam via the calibration gate -- Strip0 vs Strip1
    // inside the moments ellipse, axis-aligned at (kGateNSigmaX, kGateNSigmaY).
    Double_t g0 = ev.total[kGateStripX];
    Double_t g1 = ev.total[kGateStripY];
    if (!(g0 > 0.0 && g1 > 0.0))
      continue;
    if (!BeamFitUtils::InEllipseXY(gate, g0, g1, kGateNSigmaX, kGateNSigmaY))
      continue;

    // Condition 2: more energy deposited in the reaction strip than the one
    // before it (calibrated dE, directly comparable).
    if (!(ev.total[kReacStrip] > ev.total[kPrevStrip]))
      continue;
    n_gated++;

    // Calibrated energies are already normalized to NORM_MUSIC_MEV per channel,
    // so just sum them (long side; short sides are 0 in MeV).
    Double_t xsum = 0.0;
    for (Int_t s = kXLo; s <= kXHi; s++)
      xsum += ev.total[s];
    Double_t ysum = 0.0;
    for (Int_t s = kYLo; s <= kYHi; s++)
      ysum += ev.total[s];
    h->Fill(xsum, ysum);
  }
  std::cout << "  " << n_gated << " events passed the PID gate" << std::endl;

  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h, c);
    h->GetYaxis()->SetTitleOffset(1.3);
    c->SetLeftMargin(0.18);
    PlottingUtils::SaveFigure(
        c, Form("normsumE_s%d_%d_vs_s%d_%d", kYLo, kYHi, kXLo, kXHi), subdir,
        PlotSaveOptions::kLINEAR);
    delete c;
  }
  delete h;
}

void StripSumScatter::Run() {
  // We sum long-side per-strip energies only, so the short sides must be
  // excluded from the per-strip totals (collapsed to the long side in MeV).
  if (!Constants::IGNORE_SHORT_STRIPS)
    throw std::runtime_error(
        "StripSumScatter: Constants::IGNORE_SHORT_STRIPS must be kTRUE "
        "(long-side per-strip sum; short sides must be excluded).");

  const TString project_root = Paths::DatasetDir();
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");

  // Group cal sidecars by run; each run becomes one TChain over events_cal.
  std::vector<FileSpec> all_specs = FileSet::BuildFileSpecs();
  std::map<Int_t, TChain *> chain_by_run;
  std::vector<Int_t> run_order;
  for (std::size_t i = 0; i < all_specs.size(); i++) {
    const FileSpec &s = all_specs[i];
    TString full = IO::GetRootFilesBaseDir() + "/" + FileSet::CalSidecarName(s);
    if (gSystem->AccessPathName(full)) {
      std::cerr << "Missing cal sidecar: " << full << std::endl;
      continue;
    }
    if (chain_by_run.find(s.run) == chain_by_run.end()) {
      chain_by_run[s.run] = new TChain("events_cal");
      run_order.push_back(s.run);
    }
    chain_by_run[s.run]->Add(full);
  }

  for (std::size_t i = 0; i < run_order.size(); i++) {
    Int_t run = run_order[i];
    TString subdir = Form("strip_sum_scatter/run%d", run);
    ProcessRun(run, chain_by_run[run], subdir);
  }

  for (std::size_t i = 0; i < run_order.size(); i++)
    delete chain_by_run[run_order[i]];
}
