#include "CalibrateBeamInternal.hpp"

BeamFit2D FindBeamGateStp2VsStp1(const FileSpec &spec, const TString &run_label,
                                 const TString &plot_subdir, Bool_t save_plot) {
  BeamFit2D out;

  TString sub = FileSet::EventsName(spec) + ".root";
  TFile *sf = IO::OpenForReading(sub);
  if (!sf || sf->IsZombie()) {
    if (sf)
      delete sf;
    return out;
  }
  TTree *tree = static_cast<TTree *>(sf->Get("events"));
  if (!tree) {
    sf->Close();
    delete sf;
    return out;
  }
  // Raw ADC, read before any calibration exists (this fit produces it), so read
  // the branches directly. strip1 total = L1 + R1; strip2 total = L2 + R2.
  UShort_t left_0_17_adc[18], rightdE_adc[18];
  tree->SetBranchAddress("Left_0_17_dE", left_0_17_adc);
  tree->SetBranchAddress("RightdE", rightdE_adc);

  // Use 1024 bins (previously 256) so each bin is narrower (16 ADC for 37Cl).
  // The sigma floor in ComputeMoments is 2 * bin_width; at the old 64 ADC/bin
  // it clamped sigma to 128 ADC, hiding the true beam width and washing out
  // correlation (rho ≈ 0.09 even though strip1/strip2 track the same beam).
  const Int_t kBeamGateNBins = 1024;
  TH2F *h = new TH2F(Form("h2_stp2_vs_stp1_%s", run_label.Data()),
                     ";Strip1 #DeltaE [ADC];Strip2 #DeltaE [ADC]",
                     kBeamGateNBins, 0.0, Constants::cfg.STRIP_E_MAX_ADC,
                     kBeamGateNBins, 0.0, Constants::cfg.STRIP_E_MAX_ADC);
  h->SetDirectory(nullptr);
  Long64_t n = tree->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
    tree->GetEntry(j);
    Int_t stp1 = Int_t(left_0_17_adc[1]) + Int_t(rightdE_adc[1]);
    Int_t stp2 = Int_t(left_0_17_adc[2]) + Int_t(rightdE_adc[2]);
    if (stp1 > 0 && stp2 > 0)
      h->Fill(Double_t(stp1), Double_t(stp2));
  }
  sf->Close();
  delete sf;

  if (h->GetEntries() < 100) {
    std::cerr << "  " << run_label
              << ": too few events for Strip2-vs-Strip1 beam gate" << std::endl;
    delete h;
    return out;
  }

  const Double_t kSeedFrac = 0.30;
  // With 1024 bins instead of 256, scale from 10→40 to keep the ADC seed
  // window roughly 640 ADC (10 bins × 16384/256 = 640; 40 bins × 16384/1024 =
  // 640).
  const Int_t kSeedHalfBins = 40;
  const Int_t kMomentRefineIters = 4;
  const Double_t kMomentRefineNSigma = 2.5;
  Double_t bw_x = h->GetXaxis()->GetBinWidth(1);
  Double_t bw_y = h->GetYaxis()->GetBinWidth(1);
  Int_t bx, by, bz;
  h->GetMaximumBin(bx, by, bz);
  Double_t peak_val = h->GetBinContent(bx, by);
  Int_t lo_bx = std::max(1, bx - kSeedHalfBins);
  Int_t hi_bx = std::min(h->GetNbinsX(), bx + kSeedHalfBins);
  Int_t lo_by = std::max(1, by - kSeedHalfBins);
  Int_t hi_by = std::min(h->GetNbinsY(), by + kSeedHalfBins);
  Moments2D m = BeamFitUtils::ComputeMoments(h, lo_bx, hi_bx, lo_by, hi_by,
                                             kSeedFrac * peak_val, bw_x, bw_y);
  if (m.weight <= 0) {
    std::cerr << "  " << run_label << ": no bins above beam seed threshold"
              << std::endl;
    delete h;
    return out;
  }
  // Iteratively re-center: recompute the moments inside a ±2.5σ window
  // around the current centroid. At high rate (run84: 45 kHz) the pileup
  // blob at ~2x the beam and the correlated diagonal band both pass the
  // seed threshold; a single wide-window pass then reports a huge, highly
  // correlated pseudo-blob (sigma ~650, rho 0.95) instead of the beam. The
  // shrinking window converges onto the dominant (beam) blob.
  for (Int_t iter = 0; iter < kMomentRefineIters; iter++) {
    Int_t wlo_bx = std::max(
        1, h->GetXaxis()->FindBin(m.mu_x - kMomentRefineNSigma * m.sigma_x));
    Int_t whi_bx = std::min(
        h->GetNbinsX(),
        h->GetXaxis()->FindBin(m.mu_x + kMomentRefineNSigma * m.sigma_x));
    Int_t wlo_by = std::max(
        1, h->GetYaxis()->FindBin(m.mu_y - kMomentRefineNSigma * m.sigma_y));
    Int_t whi_by = std::min(
        h->GetNbinsY(),
        h->GetYaxis()->FindBin(m.mu_y + kMomentRefineNSigma * m.sigma_y));
    Moments2D m_ref = BeamFitUtils::ComputeMoments(
        h, wlo_bx, whi_bx, wlo_by, whi_by, kSeedFrac * peak_val, bw_x, bw_y);
    if (m_ref.weight <= 0)
      break;
    m = m_ref;
  }
  out.amp = peak_val;
  out.mu_x = m.mu_x;
  out.mu_y = m.mu_y;
  out.sigma_x = m.sigma_x;
  out.sigma_y = m.sigma_y;
  out.rho = m.rho;
  out.ok = kTRUE;
  std::cout << "  beam gate (Strip2 vs Strip1): mu=(" << out.mu_x << ","
            << out.mu_y << ") sigma=(" << out.sigma_x << "," << out.sigma_y
            << ") rho=" << out.rho << std::endl;

  if (save_plot) {
    TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h, cv);
    // Draw the correlated 2D Gaussian ellipse: the TEllipse rotates
    // according to the covariance eigen-decomposition so the drawn contour
    // matches the actual InEllipseXY gate (χ² < n²).
    Double_t sxx = out.sigma_x * out.sigma_x;
    Double_t syy = out.sigma_y * out.sigma_y;
    Double_t sxy = out.rho * out.sigma_x * out.sigma_y;
    Double_t sum = sxx + syy;
    Double_t diff = sxx - syy;
    Double_t det = TMath::Sqrt(diff * diff + 4.0 * sxy * sxy);
    Double_t lambda1 = 0.5 * (sum + det);
    Double_t lambda2 = 0.5 * (sum - det);
    Double_t theta = 0.5 * TMath::ATan2(2.0 * sxy, diff) * 180.0 / TMath::Pi();
    Double_t n = 0.5 * (kEllipseNSigmaX + kEllipseNSigmaY);
    TEllipse *e = new TEllipse(out.mu_x, out.mu_y, n * TMath::Sqrt(lambda1),
                               n * TMath::Sqrt(lambda2), 0, 360, theta);
    e->SetFillStyle(0);
    e->SetLineColor(kViolet + 2);
    e->SetLineWidth(2);
    e->Draw();
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(cv, "beam_gate_stp2_vs_stp1", plot_subdir,
                                PlotSaveOptions::kLINEAR);
    delete cv;
  }
  delete h;
  return out;
}

BeamFit2D FindBeamGateStp0VsGrid(const FileSpec &spec, const TString &run_label,
                                 const TString &plot_subdir,
                                 const BeamFit2D &beam, Bool_t save_plot) {
  BeamFit2D out;

  TString sub = FileSet::EventsName(spec) + ".root";
  TFile *sf = IO::OpenForReading(sub);
  if (!sf || sf->IsZombie()) {
    if (sf)
      delete sf;
    return out;
  }
  TTree *tree = static_cast<TTree *>(sf->Get("events"));
  if (!tree) {
    sf->Close();
    delete sf;
    return out;
  }
  UShort_t left_0_17_adc[18], rightdE_adc[18];
  Short_t grid_adc = 0;
  tree->SetBranchAddress("Left_0_17_dE", left_0_17_adc);
  tree->SetBranchAddress("RightdE", rightdE_adc);
  tree->SetBranchAddress("Grid", &grid_adc);

  const Int_t kBeamGateNBins = 256;
  const Double_t kGridMaxADC = 4096.0;
  const Double_t kStrip0MaxADC = 1024.0;
  TH2F *h = new TH2F(Form("h2_stp0_vs_grid_%s", run_label.Data()),
                     ";Grid #DeltaE [ADC];Strip0 #DeltaE [ADC]", kBeamGateNBins,
                     0.0, kGridMaxADC, kBeamGateNBins, 0.0, kStrip0MaxADC);
  h->SetDirectory(nullptr);
  Long64_t n = tree->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
    tree->GetEntry(j);
    if (beam.ok) {
      Double_t x = Double_t(left_0_17_adc[1]) + Double_t(rightdE_adc[1]);
      Double_t y = Double_t(left_0_17_adc[2]) + Double_t(rightdE_adc[2]);
      if (x <= 0 || y <= 0)
        continue;
      if (!BeamFitUtils::InEllipseXY(beam, x, y, kEllipseNSigmaX,
                                     kEllipseNSigmaY))
        continue;
    }
    Int_t s0 = Int_t(left_0_17_adc[0]);
    Int_t g = Int_t(grid_adc);
    if (s0 > 0 && g > 0)
      h->Fill(Double_t(g), Double_t(s0));
  }
  sf->Close();
  delete sf;

  if (h->GetEntries() < 100) {
    std::cerr << "  " << run_label
              << ": too few events for Strip0-vs-Grid beam gate" << std::endl;
    delete h;
    return out;
  }

  const Double_t kSeedFrac = 0.30;
  const Int_t kSeedHalfBins = 40;
  const Int_t kMomentRefineIters = 4;
  const Double_t kMomentRefineNSigma = 2.5;
  Double_t bw_x = h->GetXaxis()->GetBinWidth(1);
  Double_t bw_y = h->GetYaxis()->GetBinWidth(1);
  Int_t bx, by, bz;
  h->GetMaximumBin(bx, by, bz);
  Double_t peak_val = h->GetBinContent(bx, by);
  Int_t lo_bx = std::max(1, bx - kSeedHalfBins);
  Int_t hi_bx = std::min(h->GetNbinsX(), bx + kSeedHalfBins);
  Int_t lo_by = std::max(1, by - kSeedHalfBins);
  Int_t hi_by = std::min(h->GetNbinsY(), by + kSeedHalfBins);
  Moments2D m = BeamFitUtils::ComputeMoments(h, lo_bx, hi_bx, lo_by, hi_by,
                                             kSeedFrac * peak_val, bw_x, bw_y);
  if (m.weight <= 0) {
    std::cerr << "  " << run_label
              << ": no bins above beam seed threshold in Strip0-vs-Grid"
              << std::endl;
    delete h;
    return out;
  }
  for (Int_t iter = 0; iter < kMomentRefineIters; iter++) {
    Int_t wlo_bx = std::max(
        1, h->GetXaxis()->FindBin(m.mu_x - kMomentRefineNSigma * m.sigma_x));
    Int_t whi_bx = std::min(
        h->GetNbinsX(),
        h->GetXaxis()->FindBin(m.mu_x + kMomentRefineNSigma * m.sigma_x));
    Int_t wlo_by = std::max(
        1, h->GetYaxis()->FindBin(m.mu_y - kMomentRefineNSigma * m.sigma_y));
    Int_t whi_by = std::min(
        h->GetNbinsY(),
        h->GetYaxis()->FindBin(m.mu_y + kMomentRefineNSigma * m.sigma_y));
    Moments2D m_ref = BeamFitUtils::ComputeMoments(
        h, wlo_bx, whi_bx, wlo_by, whi_by, kSeedFrac * peak_val, bw_x, bw_y);
    if (m_ref.weight <= 0)
      break;
    m = m_ref;
  }
  out.amp = peak_val;
  out.mu_x = m.mu_x;
  out.mu_y = m.mu_y;
  out.sigma_x = m.sigma_x;
  out.sigma_y = m.sigma_y;
  out.rho = m.rho;
  out.ok = kTRUE;
  std::cout << "  beam gate (Strip0 vs Grid): mu=(" << out.mu_x << ","
            << out.mu_y << ") sigma=(" << out.sigma_x << "," << out.sigma_y
            << ") rho=" << out.rho << std::endl;

  if (save_plot) {
    TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h, cv);
    Double_t sxx = out.sigma_x * out.sigma_x;
    Double_t syy = out.sigma_y * out.sigma_y;
    Double_t sxy = out.rho * out.sigma_x * out.sigma_y;
    Double_t sum = sxx + syy;
    Double_t diff = sxx - syy;
    Double_t det = TMath::Sqrt(diff * diff + 4.0 * sxy * sxy);
    Double_t lambda1 = 0.5 * (sum + det);
    Double_t lambda2 = 0.5 * (sum - det);
    Double_t theta = 0.5 * TMath::ATan2(2.0 * sxy, diff) * 180.0 / TMath::Pi();
    Double_t n = 0.5 * (kEllipseNSigmaX + kEllipseNSigmaY);
    TEllipse *e = new TEllipse(out.mu_x, out.mu_y, n * TMath::Sqrt(lambda1),
                               n * TMath::Sqrt(lambda2), 0, 360, theta);
    e->SetFillStyle(0);
    e->SetLineColor(kViolet + 2);
    e->SetLineWidth(2);
    e->Draw();
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(cv, "beam_gate_stp0_vs_grid", plot_subdir,
                                PlotSaveOptions::kLINEAR);
    delete cv;
  }
  delete h;
  return out;
}
