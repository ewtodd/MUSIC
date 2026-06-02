#include "StripScatter.hpp"

namespace {

// Gate: 2D bigaus fit on Strip1 long-side (L, LeftdEMeV[1]) vs Strip2
// (R, RightdEMeV[2]). Mirrors the k=1 ellipse gate used in DeltaEScatter.
const Double_t kGateMin = 0.0;
const Double_t kGateMax = 12.0;
const Int_t kGateBins = 240;
const Double_t kGateNSigma = 2.0;
const Int_t kSeedHalfBins = 40;
const Double_t kSeedFrac = 0.30;

TH2F *BuildGateHist(TChain *chain, const TString &name) {
  Float_t leftdE[18] = {0};
  Float_t rightdE[18] = {0};
  chain->SetBranchStatus("*", 0);
  chain->SetBranchStatus("LeftdEMeV", 1);
  chain->SetBranchStatus("RightdEMeV", 1);
  chain->SetBranchAddress("LeftdEMeV", leftdE);
  chain->SetBranchAddress("RightdEMeV", rightdE);

  TH2F *h = new TH2F(name, ";#DeltaE S0 [MeV];#DeltaE L1 [MeV]", kGateBins,
                     kGateMin, kGateMax, kGateBins, kGateMin, kGateMax);
  h->SetDirectory(nullptr);

  Long64_t n = chain->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
    chain->GetEntry(j);
    Double_t x = Double_t(leftdE[1]);
    Double_t y = Double_t(rightdE[2]);
    if (x > 0 && y > 0)
      h->Fill(x, y);
  }
  return h;
}

// Moments-seeded bigaus fit on the single largest peak.
BeamFit2D FitBigPeak(TH2F *h, const TString &tag) {
  BeamFit2D out;
  if (!h || h->GetEntries() < 100)
    return out;

  Double_t bw_x = h->GetXaxis()->GetBinWidth(1);
  Double_t bw_y = h->GetYaxis()->GetBinWidth(1);

  Int_t bx = 0, by = 0, bz = 0;
  h->GetMaximumBin(bx, by, bz);
  Double_t peak_val = h->GetBinContent(bx, by);
  if (peak_val <= 0)
    return out;

  Int_t lo_bx = std::max(1, bx - kSeedHalfBins);
  Int_t hi_bx = std::min(h->GetNbinsX(), bx + kSeedHalfBins);
  Int_t lo_by = std::max(1, by - kSeedHalfBins);
  Int_t hi_by = std::min(h->GetNbinsY(), by + kSeedHalfBins);

  Moments2D m = BeamFitUtils::ComputeMoments(h, lo_bx, hi_bx, lo_by, hi_by,
                                             kSeedFrac * peak_val, bw_x, bw_y);
  if (m.weight <= 0) {
    std::cerr << "  FitBigPeak " << tag << ": no bins above seed threshold"
              << std::endl;
    return out;
  }

  Double_t x_lo = std::max(h->GetXaxis()->GetXmin(), m.mu_x - 3.0 * m.sigma_x);
  Double_t x_hi = std::min(h->GetXaxis()->GetXmax(), m.mu_x + 3.0 * m.sigma_x);
  Double_t y_lo = std::max(h->GetYaxis()->GetXmin(), m.mu_y - 3.0 * m.sigma_y);
  Double_t y_hi = std::min(h->GetYaxis()->GetXmax(), m.mu_y + 3.0 * m.sigma_y);

  out = BeamFitUtils::FitBigausInWindow(h, x_lo, x_hi, y_lo, y_hi, m, peak_val,
                                        Form("f_gate_%s", tag.Data()),
                                        TString("FitBigPeak ") + tag);

  std::cout << "  gate fit " << tag << ": mu=(" << out.mu_x << "," << out.mu_y
            << ") sigma=(" << out.sigma_x << "," << out.sigma_y
            << ") rho=" << out.rho << std::endl;
  return out;
}

// x = strip #DeltaE, y = full-event total #DeltaE.
TH2F *NewStripHist(Int_t strip, Int_t run, const char *suffix) {
  TH2F *h =
      new TH2F(Form("h2_totalE_vs_stripE_s%d_run%d%s", strip, run, suffix),
               Form(";Strip %d #DeltaE [MeV];Total #DeltaE [MeV]", strip), 200,
               Constants::STRIP_E_MIN_MEV, Constants::STRIP_E_MAX_MEV, 400,
               Constants::TOTAL_E_MIN_MEV, Constants::TOTAL_E_MAX_MEV);
  h->SetDirectory(nullptr);
  h->SetStats(0);
  return h;
}

// x = #DeltaE summed over {R-1, R, R+1}, y = #DeltaE summed from R to the
// downstream end of the detector.
TH2F *NewWindowHist(Int_t r, Int_t s_hi, Int_t run, const char *suffix) {
  TH2F *h = new TH2F(
      Form("h2_window_vs_downstream_R%d_run%d%s", r, run, suffix),
      Form(
          ";#DeltaE strips %d+%d+%d [MeV];#DeltaE strips %d#rightarrow%d [MeV]",
          r - 1, r, r + 1, r, s_hi),
      200, 3.0 * Constants::STRIP_E_MIN_MEV, 3.0 * Constants::STRIP_E_MAX_MEV,
      400, Constants::TOTAL_E_MIN_MEV, Constants::TOTAL_E_MAX_MEV);
  h->SetDirectory(nullptr);
  h->SetStats(0);
  return h;
}

} // namespace

void StripScatter::SaveH2(TH2F *h, const TString &save_name,
                          const TString &subdir) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDraw2DHistogram(h, c);
  h->GetYaxis()->SetTitleOffset(1.3);
  c->SetLeftMargin(0.18);
  PlottingUtils::SaveFigure(c, save_name, subdir, PlotSaveOptions::kLINEAR);
  delete c;
}

void StripScatter::SaveGatePlot(TH2F *gate_hist, const BeamFit2D &gate,
                                const TString &subdir) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDraw2DHistogram(gate_hist, c);
  if (gate.ok) {
    Double_t major, minor, theta_deg;
    BeamFitUtils::DiagonalizeCov(gate, major, minor, theta_deg);
    TEllipse *e = new TEllipse(gate.mu_x, gate.mu_y, kGateNSigma * major,
                               kGateNSigma * minor, 0.0, 360.0, theta_deg);
    e->SetFillStyle(0);
    e->SetLineColor(kRed + 1);
    e->SetLineWidth(2);
    e->Draw();
  }
  PlottingUtils::SaveFigure(c, "gate_S0_vs_L1", subdir,
                            PlotSaveOptions::kLINEAR);
  delete c;
}

void StripScatter::ProcessRun(Int_t run, TChain *chain, const TString &subdir) {
  if (!chain || chain->GetEntries() == 0) {
    std::cerr << "Run " << run << ": empty chain; skipping" << std::endl;
    return;
  }

  TString tag = Form("run%d", run);
  TH2F *gate_hist = BuildGateHist(chain, Form("h2_gate_S0_L1_%s", tag.Data()));
  BeamFit2D gate = FitBigPeak(gate_hist, tag);
  SaveGatePlot(gate_hist, gate, subdir);
  if (!gate.ok) {
    std::cerr << "Run " << run << ": gate fit failed; skipping scatter"
              << std::endl;
    delete gate_hist;
    return;
  }

  Float_t leftdE[18] = {0};
  Float_t rightdE[18] = {0};
  Float_t totaldE[18] = {0};
  Float_t cathode = 0;
  chain->SetBranchStatus("*", 0);
  chain->SetBranchStatus("LeftdEMeV", 1);
  chain->SetBranchStatus("RightdEMeV", 1);
  chain->SetBranchStatus("TotaldEMeV", 1);
  chain->SetBranchStatus("CathodeMeV", 1);
  chain->SetBranchAddress("LeftdEMeV", leftdE);
  chain->SetBranchAddress("RightdEMeV", rightdE);
  chain->SetBranchAddress("TotaldEMeV", totaldE);
  chain->SetBranchAddress("CathodeMeV", &cathode);

  const Int_t s_lo = Constants::IGNORE_STRIP_0 ? 1 : 0;
  const Int_t s_hi = Constants::IGNORE_STRIP_17 ? 16 : 17;
  const Int_t r_lo = std::max(Constants::REACTION_STRIP_MIN, s_lo);
  const Int_t r_hi = std::min(Constants::REACTION_STRIP_MAX, s_hi);

  // Three parallel sets per candidate strip R: all gated events, gated events
  // with cathode == 0, and gated events with cathode != 0. The cathode-split
  // sets are only saved if the run actually has nonzero-cathode events.
  std::map<Int_t, TH2F *> h_strip, h_strip_c0, h_strip_cn;
  std::map<Int_t, TH2F *> h_window, h_window_c0, h_window_cn;
  for (Int_t r = r_lo; r <= r_hi; r++) {
    h_strip[r] = NewStripHist(r, run, "");
    h_strip_c0[r] = NewStripHist(r, run, "_c0");
    h_strip_cn[r] = NewStripHist(r, run, "_cn");
    h_window[r] = NewWindowHist(r, s_hi, run, "");
    h_window_c0[r] = NewWindowHist(r, s_hi, run, "_c0");
    h_window_cn[r] = NewWindowHist(r, s_hi, run, "_cn");
  }

  Long64_t n = chain->GetEntries();
  std::cout << "Run " << run << ": aggregating " << n
            << " events over candidate strips " << r_lo << "-" << r_hi
            << " (gated)..." << std::endl;
  Long64_t n_gated = 0;
  Bool_t has_nonzero_cathode = kFALSE;
  for (Long64_t j = 0; j < n; j++) {
    chain->GetEntry(j);
    Double_t xg = Double_t(leftdE[1]);
    Double_t yg = Double_t(rightdE[2]);
    if (!(xg > 0 && yg > 0))
      continue;
    if (!BeamFitUtils::InEllipse(gate, xg, yg, kGateNSigma))
      continue;
    n_gated++;
    Bool_t cathode_zero = (cathode == 0);
    if (!cathode_zero)
      has_nonzero_cathode = kTRUE;

    Double_t total = 0.0;
    for (Int_t s = s_lo; s <= s_hi; s++)
      total += Double_t(totaldE[s]);
    for (Int_t r = r_lo; r <= r_hi; r++) {
      Double_t window = 0.0;
      for (Int_t s = r - 1; s <= r + 1; s++)
        if (s >= s_lo && s <= s_hi)
          window += Double_t(totaldE[s]);
      Double_t downstream = 0.0;
      for (Int_t s = r; s <= s_hi; s++)
        downstream += Double_t(totaldE[s]);

      h_strip[r]->Fill(Double_t(totaldE[r]), total);
      h_window[r]->Fill(window, downstream);
      if (cathode_zero) {
        h_strip_c0[r]->Fill(Double_t(totaldE[r]), total);
        h_window_c0[r]->Fill(window, downstream);
      } else {
        h_strip_cn[r]->Fill(Double_t(totaldE[r]), total);
        h_window_cn[r]->Fill(window, downstream);
      }
    }
  }
  std::cout << "  " << n_gated << " events passed the gate"
            << (has_nonzero_cathode ? "" : " (cathode all zero; no split)")
            << std::endl;

  // Namespace per candidate strip so each strip's plots are self-contained.
  for (Int_t r = r_lo; r <= r_hi; r++) {
    TString cand_subdir = subdir + Form("/candidate_s%d", r);
    SaveH2(h_strip[r], Form("totalE_vs_stripE_s%d", r), cand_subdir);
    SaveH2(h_window[r], "window_vs_downstream", cand_subdir);
    if (has_nonzero_cathode) {
      TString c0_subdir = cand_subdir + "/cathode0";
      TString cn_subdir = cand_subdir + "/cathode_nonzero";
      SaveH2(h_strip_c0[r], Form("totalE_vs_stripE_s%d", r), c0_subdir);
      SaveH2(h_window_c0[r], "window_vs_downstream", c0_subdir);
      SaveH2(h_strip_cn[r], Form("totalE_vs_stripE_s%d", r), cn_subdir);
      SaveH2(h_window_cn[r], "window_vs_downstream", cn_subdir);
    }
    delete h_strip[r];
    delete h_window[r];
    delete h_strip_c0[r];
    delete h_window_c0[r];
    delete h_strip_cn[r];
    delete h_window_cn[r];
  }
  delete gate_hist;
}

void StripScatter::Run() {
  const TString project_root = Paths::DatasetDir();
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, project_root + "/plots",
                                project_root + "/root_files");

  // Group cal sidecars by run; each run becomes one TChain that drives the
  // per-run aggregated, gated 2D histograms.
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
    TString subdir = Form("strip_scatter/run%d", run);
    ProcessRun(run, chain_by_run[run], subdir);
  }

  for (std::size_t i = 0; i < run_order.size(); i++)
    delete chain_by_run[run_order[i]];
}
