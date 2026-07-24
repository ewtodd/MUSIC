#include "CalibrateBeam.hpp"
#include <Rtypes.h>

const Int_t kMaxChannels = 35;
const Double_t kEllipseNSigmaX = 3;
const Double_t kEllipseNSigmaY = 3;
const Long64_t kMinSamples = 200;
const Long64_t kSampleCap = 20000;
// (L, R) pair cap for the L/R gain-matching passes. Larger than kSampleCap
// because the short-side "shoulder" anchor is found in a narrow slice of the
// pairs (long side ≈ 300-350 ADC) that only a small fraction of events
// populate.
const Long64_t kPairCap = 100000;

// Defined below; forward-declared so both beam-peak fitters and the beam_peak
// diagnostic plot share one histogram recipe regardless of definition order.
void RobustPeakSeed(const std::vector<Float_t> &v, Double_t &mode,
                    Double_t &sigma);
TH1F *MakeBeamPeakHist(const TString &name, const TString &title,
                       const std::vector<Float_t> &v, Double_t mode,
                       Double_t sigma);

std::vector<ChannelCal> CalibrateBeam::BuildChannels() {
  std::vector<ChannelCal> chans;
  ChannelCal c{};
  c.name = "Strip0";
  c.side = 'S';
  c.strip = 0;
  chans.push_back(c);
  c.name = "Strip17";
  c.side = 'S';
  c.strip = 17;
  chans.push_back(c);
  for (Int_t s = 1; s <= 16; s++) {
    c.name = Form("L%d", s);
    c.side = 'L';
    c.strip = s;
    chans.push_back(c);
    c.name = Form("R%d", s);
    c.side = 'R';
    c.strip = s;
    chans.push_back(c);
  }
  c.name = "Cathode";
  c.side = 'C';
  c.strip = -1;
  chans.push_back(c);
  return chans;
}

Int_t ChannelToEresIndex(const ChannelCal &c) {
  if (c.side == 'C')
    return 0;
  if (c.side == 'S' && c.strip == 0)
    return 1;
  if (c.side == 'S' && c.strip == 17)
    return 2;
  if (c.side == 'L' && c.strip >= 1 && c.strip <= 16)
    return 3 + (c.strip - 1);
  if (c.side == 'R' && c.strip >= 1 && c.strip <= 16)
    return 19 + (c.strip - 1);
  return -1;
}

Char_t LongSide(Int_t strip) { return (strip % 2 == 0) ? 'R' : 'L'; }

Bool_t IsBeamdEChannel(const ChannelCal &c) {
  if (c.side == 'S')
    return kTRUE;
  if (c.side != 'L' && c.side != 'R')
    return kFALSE;
  if (c.strip < 1 || c.strip > 16)
    return kFALSE;
  return c.side == LongSide(c.strip);
}

Bool_t IsCalibrated(const ChannelCal &c) { return c.fit_adc > 0; }

Double_t Gain(const ChannelCal &c) {
  return c.gain >= 0.0 ? c.gain : 1.0 / c.fit_adc;
}

Double_t ResolutionFWHMPercent(const ChannelCal &c) {
  if (c.fit_adc <= 0)
    return 0.0;
  const Double_t kFwhmPerSigma = 2.0 * TMath::Sqrt(2.0 * TMath::Log(2.0));
  return 100.0 * kFwhmPerSigma * c.fit_sigma_adc / c.fit_adc;
}

inline Double_t ApplyCal(const ChannelCal &c, Double_t adc) {
  return Gain(c) * adc;
}

BeamFit2D FindBeamGateStp2VsStp1(const FileSpec &spec, const TString &run_label,
                                 const TString &plot_subdir,
                                 Bool_t save_plot = kTRUE) {
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
                                 const BeamFit2D &beam,
                                 Bool_t save_plot = kTRUE) {
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

inline Double_t Median(std::vector<Float_t> &v) {
  if (v.empty())
    return 0.0;
  Int_t n = Int_t(v.size());
  std::nth_element(v.begin(), v.begin() + n / 2, v.end());
  Double_t med = Double_t(v[n / 2]);
  if (n % 2 == 0) {
    std::nth_element(v.begin(), v.begin() + n / 2 - 1, v.end());
    med = 0.5 * (med + Double_t(v[n / 2 - 1]));
  }
  return med;
}

// IQR = Q3 - Q1; divide by 1.349 outside this helper for the Gaussian-sigma
// approximation when needed.
inline Double_t InterquartileRange(std::vector<Float_t> &v) {
  if (v.size() < 4)
    return 0.0;
  Int_t n = Int_t(v.size());
  Int_t i1 = n / 4;
  Int_t i3 = (3 * n) / 4;
  std::nth_element(v.begin(), v.begin() + i1, v.end());
  Double_t q1 = Double_t(v[i1]);
  std::nth_element(v.begin(), v.begin() + i3, v.end());
  Double_t q3 = Double_t(v[i3]);
  return q3 - q1;
}

// Robust peak/width estimate used to seed the beam-peak Gaussian fit, and as
// the fallback anchor when the fit fails. Histograms only the
// 5th-95th-percentile core so outlier ADC values can't stretch the binning,
// takes the modal bin centre as the peak and IQR/1.349 as the width. Peak-like:
// unlike the raw sample mean it is not pulled up by the straggling beam-dE
// tail.
void RobustPeakSeed(const std::vector<Float_t> &v, Double_t &mode,
                    Double_t &sigma) {
  mode = 0.0;
  sigma = 0.0;
  Int_t n = Int_t(v.size());
  if (n < 4)
    return;
  std::vector<Float_t> s(v);
  std::sort(s.begin(), s.end());
  Int_t i_lo = Int_t(0.05 * n);
  Int_t i_hi = Int_t(0.95 * n);
  Int_t i_q1 = Int_t(0.25 * n);
  Int_t i_q3 = Int_t(0.75 * n);
  if (i_hi >= n)
    i_hi = n - 1;
  Double_t p_lo = Double_t(s[i_lo]);
  Double_t p_hi = Double_t(s[i_hi]);
  Double_t med = Double_t(s[n / 2]);
  sigma = (Double_t(s[i_q3]) - Double_t(s[i_q1])) / 1.349;
  if (!(sigma > 0.0))
    sigma = 0.05 * (med > 0.0 ? med : 1.0);
  if (!(p_hi > p_lo)) {
    mode = med;
    return;
  }
  const Int_t nbins = 100;
  TH1F h("robust_peak_seed_h", "", nbins, p_lo, p_hi);
  h.SetDirectory(nullptr);
  for (Int_t j = 0; j < n; j++)
    if (Double_t(s[j]) >= p_lo && Double_t(s[j]) <= p_hi)
      h.Fill(Double_t(s[j]));
  mode = h.GetBinCenter(h.GetMaximumBin());
  if (!(mode > 0.0))
    mode = med;
}

// Per-channel beam-peak ADC histogram, shared by the fitters and the beam_peak
// diagnostic plot so a fitted Gaussian's amplitude (counts per bin) always
// matches the histogram it is drawn over. Range and bin count come from the
// robust mode/width (mode +/- 4 sigma, ~6 bins/sigma, never finer than 1 ADC),
// so the binning is chosen consistently per channel instead of from a fixed
// count over the outlier-stretched min..max range. Caller owns the histogram.
TH1F *MakeBeamPeakHist(const TString &name, const TString &title,
                       const std::vector<Float_t> &v, Double_t mode,
                       Double_t sigma) {
  Double_t xlo = mode - 4.0 * sigma;
  Double_t xhi = mode + 4.0 * sigma;
  if (xlo < 0.0)
    xlo = 0.0;
  Double_t bin_target = TMath::Max(1.0, sigma / 6.0);
  Int_t nbins = Int_t((xhi - xlo) / bin_target);
  if (nbins < 20)
    nbins = 20;
  if (nbins > 200)
    nbins = 200;
  TH1F *h = new TH1F(name, title, nbins, xlo, xhi);
  h->SetDirectory(nullptr);
  for (Int_t j = 0; j < Int_t(v.size()); j++)
    h->Fill(Double_t(v[j]));
  return h;
}

// Beam-peak Gaussian fit, robustly seeded: estimate the peak/width from the
// percentile-clipped mode (RobustPeakSeed), then fit "gaus" over only the peak
// core (mode +/- 2 sigma). A second pass refits inside a NARROW window around
// the first centroid: when the spectrum is bimodal (e.g. Strip0 with a
// secondary bump below the beam peak), the IQR-based seed sigma is inflated
// by the contamination and the wide first-pass window lets the fit average
// both components — the narrow refit locks onto the dominant peak. This is
// the primary (and only) fit in ReduceToAnchors.
Bool_t FitBeamPeakGaussian(const std::vector<Float_t> &v, const TString &fname,
                           Double_t &peak_adc, Double_t &sigma_adc,
                           TF1 *&fit_out) {
  fit_out = nullptr;
  if (v.size() < kMinSamples)
    return kFALSE;
  Double_t mode = 0.0, rsigma = 0.0;
  RobustPeakSeed(v, mode, rsigma);
  if (!(mode > 0.0) || !(rsigma > 0.0))
    return kFALSE;

  TH1F *h = MakeBeamPeakHist(fname + "_h", "", v, mode, rsigma);
  Double_t xlo = h->GetXaxis()->GetXmin();
  Double_t xhi = h->GetXaxis()->GetXmax();
  Double_t bw = h->GetBinWidth(1);
  Double_t fit_lo = mode - 2.0 * rsigma;
  Double_t fit_hi = mode + 2.0 * rsigma;
  if (fit_lo < xlo)
    fit_lo = xlo;
  if (fit_hi > xhi)
    fit_hi = xhi;
  Double_t amp_seed = h->GetBinContent(h->FindBin(mode));

  TF1 *f = new TF1(fname, "gaus", fit_lo, fit_hi);
  f->SetNpx(1000);
  f->SetParameters(amp_seed, mode, rsigma);
  f->SetParLimits(1, fit_lo, fit_hi);
  f->SetParLimits(2, bw, fit_hi - fit_lo);
  TFitResultPtr r = h->Fit(f, "QSRNL");
  if (!r.Get() || !r->IsValid()) {
    delete f;
    delete h;
    return kFALSE;
  }
  peak_adc = f->GetParameter(1);
  sigma_adc = std::fabs(f->GetParameter(2));

  // Second pass: refit inside mu ± min(1.5*sigma_fit, 12% of mu). For a
  // clean single Gaussian this window still spans the core and reproduces
  // the first-pass result; for a contaminated spectrum it excludes the
  // secondary component and re-centres onto the dominant peak.
  Double_t half = 1.5 * sigma_adc;
  Double_t half_cap = 0.12 * peak_adc;
  if (half_cap < half)
    half = half_cap;
  if (half > 2.0 * bw) {
    Double_t lo2 = peak_adc - half;
    Double_t hi2 = peak_adc + half;
    if (lo2 < xlo)
      lo2 = xlo;
    if (hi2 > xhi)
      hi2 = xhi;
    TF1 *f2 = new TF1(fname + "_p2", "gaus", lo2, hi2);
    f2->SetNpx(1000);
    f2->SetParameters(f->GetParameter(0), peak_adc, sigma_adc);
    f2->SetParLimits(1, lo2, hi2);
    f2->SetParLimits(2, bw, hi2 - lo2);
    TFitResultPtr r2 = h->Fit(f2, "QSRNL");
    if (r2.Get() && r2->IsValid() && f2->GetParameter(1) > 0 &&
        std::fabs(f2->GetParameter(2)) > 0) {
      peak_adc = f2->GetParameter(1);
      sigma_adc = std::fabs(f2->GetParameter(2));
      delete f;
      f = f2;
    } else {
      delete f2;
    }
  }
  delete h;
  if (!(peak_adc > 0.0) || !(sigma_adc > 0.0)) {
    delete f;
    return kFALSE;
  }
  fit_out = f;
  return kTRUE;
}

// Fit a Gaussian to the bucket and return (mu, sigma). Falls back to
// (median, IQR/1.349) on fit failure. Sim per-channel deposits are
// well-approximated by a Gaussian, so a direct fit gives a cleaner
// (mu, sigma) than median/IQR estimators.
Bool_t FitGaussianMuSigma(const std::vector<Float_t> &v, const TString &fname,
                          Double_t &mu, Double_t &sigma) {
  if (v.size() < 50)
    return kFALSE;
  Float_t lo = v[0], hi = v[0];
  for (Int_t j = 1; j < Int_t(v.size()); j++) {
    if (v[j] < lo)
      lo = v[j];
    if (v[j] > hi)
      hi = v[j];
  }
  Double_t pad = 0.05 * (Double_t(hi) - Double_t(lo));
  if (pad < 1e-6)
    pad = 1e-6;
  const Int_t nbins = 75;
  TH1F h(fname + "_h", "", nbins, Double_t(lo) - pad, Double_t(hi) + pad);
  h.SetDirectory(nullptr);
  for (Int_t j = 0; j < Int_t(v.size()); j++)
    h.Fill(Double_t(v[j]));
  TF1 fg(fname, "gaus", Double_t(lo) - pad, Double_t(hi) + pad);
  Int_t pb = h.GetMaximumBin();
  fg.SetParameters(h.GetBinContent(pb), h.GetBinCenter(pb), h.GetRMS());
  TFitResultPtr r = h.Fit(&fg, "QSRN");
  if (!r.Get() || !r->IsValid())
    return kFALSE;
  mu = fg.GetParameter(1);
  sigma = std::fabs(fg.GetParameter(2));
  return mu > 0 && sigma > 0;
}

// Paired (L, R) raw-ADC samples for one strip, collected UNGATED from events
// where both ends fire, plus the uncapped "shoulder" slice: short-side values
// from events where the LONG side reads low (charge went mostly to the short
// end). The slice is collected separately because those events are rare — the
// capped pair vectors fill up with beam events long before enough slice
// events arrive.
struct StripPairSamples {
  std::vector<Float_t> l;
  std::vector<Float_t> r;
  std::vector<Float_t> shoulder;
};

// check_LR gain-match slice constants (37Cl_an_check_LR.ipynb, cell "Save
// 2-pass calibration"). Defined before CollectAnchorSamplesOneSubfile because
// the shoulder slice is selected during sample collection.
const Double_t kGmOffsetLongL = 350.0; // ADC slice centre when L is long
const Double_t kGmOffsetLongR = 300.0; // ADC slice centre when R is long
const Double_t kGmOffsetWin = 80.0;    // ADC slice half-width
const Double_t kGmPeakWin = 150.0;     // ADC near-peak gate half-width

void CollectAnchorSamplesOneSubfile(const FileSpec &spec,
                                    const std::vector<ChannelCal> &chans,
                                    const BeamFit2D &beam,
                                    const BeamFit2D &beam0vGrid,
                                    std::vector<std::vector<Float_t>> &samples,
                                    StripPairSamples pairs[18]) {
  Int_t n_chans = Int_t(chans.size());
  samples.assign(n_chans, std::vector<Float_t>());
  if (!beam.ok)
    return;

  TString sub = FileSet::EventsName(spec) + ".root";
  TFile *sf = IO::OpenForReading(sub);
  if (!sf || sf->IsZombie()) {
    if (sf)
      sf->Close();
    return;
  }
  TTree *tree = static_cast<TTree *>(sf->Get("events"));
  if (!tree) {
    sf->Close();
    delete sf;
    return;
  }
  // Raw ADC, pre-calibration. Guard strips (S) and left ends (L) live in
  // Left_0_17_dE; right ends in RightdE. Strip totals are L+R; the gate uses
  // the strip1 total (L1+R1) and the strip2 total (L2+R2).
  UShort_t left_0_17_adc[18], rightdE_adc[18];
  Short_t cathode_adc = 0, grid_adc = 0;
  tree->SetBranchAddress("Left_0_17_dE", left_0_17_adc);
  tree->SetBranchAddress("RightdE", rightdE_adc);
  tree->SetBranchAddress("Cathode", &cathode_adc);
  tree->SetBranchAddress("Grid", &grid_adc);

  Long64_t n = tree->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
    tree->GetEntry(j);
    // (L, R) pairs for strips 1–16 where both ends fired — UNGATED, matching
    // the check_LR notebook which runs on all events. The short-side anchor
    // ("shoulder") is found in a slice where the LONG side reads low
    // (~300-350 ADC), i.e. events where the charge went mostly to the short
    // end. Those are reaction/off-position events that a beam gate would
    // remove, so the pairs must not be beam-gated.
    if (pairs) {
      for (Int_t s = 1; s <= 16; s++) {
        Int_t lv = Int_t(left_0_17_adc[s]);
        Int_t rv = Int_t(rightdE_adc[s]);
        if (lv > 0 && rv > 0 && Long64_t(pairs[s].l.size()) < kPairCap) {
          pairs[s].l.push_back(Float_t(lv));
          pairs[s].r.push_back(Float_t(rv));
        }
        // Shoulder slice: short-side value when the long side reads ≈ the
        // offset. Collected independently of the capped pairs above so rare
        // slice events are never crowded out by beam statistics.
        Bool_t l_is_long = (LongSide(s) == 'L');
        Int_t long_v = l_is_long ? lv : rv;
        Int_t short_v = l_is_long ? rv : lv;
        Double_t offset = l_is_long ? kGmOffsetLongL : kGmOffsetLongR;
        if (short_v > 0 &&
            TMath::Abs(Double_t(long_v) - offset) < kGmOffsetWin &&
            Long64_t(pairs[s].shoulder.size()) < kPairCap)
          pairs[s].shoulder.push_back(Float_t(short_v));
      }
    }
    Double_t x = Double_t(left_0_17_adc[1]) + Double_t(rightdE_adc[1]);
    Double_t y = Double_t(left_0_17_adc[2]) + Double_t(rightdE_adc[2]);
    if (x <= 0 || y <= 0)
      continue;
    if (!BeamFitUtils::InEllipseXY(beam, x, y, kEllipseNSigmaX,
                                   kEllipseNSigmaY))
      continue;
    if (beam0vGrid.ok) {
      Double_t g0 = Double_t(grid_adc);
      Double_t s0 = Double_t(left_0_17_adc[0]);
      if (g0 <= 0 || s0 <= 0)
        continue;
      if (!BeamFitUtils::InEllipseXY(beam0vGrid, g0, s0, kEllipseNSigmaX,
                                     kEllipseNSigmaY))
        continue;
    }
    for (Int_t i = 0; i < n_chans; i++) {
      if (Long64_t(samples[i].size()) >= kSampleCap)
        continue;
      const ChannelCal &c = chans[i];
      Int_t v = 0;
      if (c.side == 'S' || c.side == 'L')
        v = Int_t(left_0_17_adc[c.strip]);
      else if (c.side == 'R')
        v = Int_t(rightdE_adc[c.strip]);
      else if (c.side == 'C')
        v = Int_t(cathode_adc);
      if (v > 0)
        samples[i].push_back(Float_t(v));
    }
  }
  sf->Close();
  delete sf;
}

// L/R gain matching, replicating the check_LR notebook's two-pass recipe
// (37Cl_an_check_LR.ipynb, "Save 2-pass calibration" cell) exactly, minus the
// eta/position correction:
//
// Pass 1 — per-side anchors:
//   * LONG side anchor  = beam peak of the long side (histogram mode). The
//     per-channel anchor from ReduceToAnchors is exactly this (the beam
//     dominates the spectrum), so it is reused.
//   * SHORT side anchor = "shoulder": the mode of the short side among events
//     where the LONG side reads low (≈350 ADC for L-long strips, ≈300 for
//     R-long, ±80). Those are events where the charge went mostly to the
//     short end, so the short side is measured at the same charge scale as
//     the long side's beam peak. Search window (500, 2000) ADC.
//   * gain = TARGET / anchor per side, with TARGET = 1.0 a.u. (the notebook
//     uses 1000 ADC; only the overall scale differs).
//
// Pass 2 — per-strip eSum alignment via the SHORT side only:
//   * eSum = gain_L·L + gain_R·R over the pairs; per-strip beam peak found
//     in (0.8, 2.5) a.u. (notebook: TARGET=1000, eSum ≈1200 →
//     (1100/1000, 2500/1000) = (1.1, 2.5); widened slightly here).
//   * reference = median of the per-strip peaks.
//   * strips with |peak − ref| > 0.05 get gain_short *= (ref−1)/(peak−1),
//     which moves the strip's eSum peak onto the reference without touching
//     the long side.
//
// Pairs are collected UNGATED (all events with both ends firing) because the
// shoulder slice needs reaction/off-position events that the beam gate
// removes. Strips whose slice is still too small (upstream strips see few
// reaction products, so slice statistics grow with strip number) fall back to
// the MEDIAN shoulder of the strips that measured one — the short-side
// full-charge response is the same electronics scale on every strip, and
// pass 2 then fine-tunes each strip individually using its own eSum peak
// (full beam statistics).
const Double_t kGmShortLo = 500.0; // ADC short-side shoulder search window
const Double_t kGmShortHi = 2000.0;
const Int_t kGmBins = 512;
const Double_t kGmSkipFrac = 0.05;
const Long64_t kGmMinSlice = 100; // min entries in the shoulder slice
const Double_t kGmEsumLo = 0.8;   // a.u. eSum peak search window (pass 2)
const Double_t kGmEsumHi = 2.5;
const Double_t kGmCorrThresh = 0.05; // a.u. (notebook: 30/1000 = 3% of TARGET)

// Histogram-mode peak finder, mirroring the notebook's find_peak(): histogram
// `v` over [lo, hi] with kGmBins bins, skip the first skip_frac of bins (to
// avoid the threshold pile), return the max-bin centre. Returns 0 when empty.
Double_t GmFindPeak(const std::vector<Float_t> &v, Double_t lo, Double_t hi,
                    Double_t skip_frac) {
  if (v.empty() || hi <= lo)
    return 0.0;
  std::vector<Long64_t> h(kGmBins, 0);
  Double_t bw = (hi - lo) / kGmBins;
  for (Int_t j = 0; j < Int_t(v.size()); j++) {
    Double_t x = Double_t(v[j]);
    if (x < lo || x >= hi)
      continue;
    Int_t b = Int_t((x - lo) / bw);
    if (b >= 0 && b < kGmBins)
      h[b]++;
  }
  Int_t skip = Int_t(kGmBins * skip_frac);
  Int_t bmax = -1;
  Long64_t vmax = 0;
  for (Int_t b = skip; b < kGmBins; b++) {
    if (h[b] > vmax) {
      vmax = h[b];
      bmax = b;
    }
  }
  if (bmax < 0 || vmax <= 0)
    return 0.0;
  return lo + (bmax + 0.5) * bw;
}

void ComputeLRGainMatch(std::vector<ChannelCal> &chans,
                        const StripPairSamples pairs[18]) {
  Int_t idx_l[18], idx_r[18];
  for (Int_t s = 0; s < 18; s++) {
    idx_l[s] = -1;
    idx_r[s] = -1;
  }
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    if (chans[i].strip >= 1 && chans[i].strip <= 16) {
      if (chans[i].side == 'L')
        idx_l[chans[i].strip] = i;
      else if (chans[i].side == 'R')
        idx_r[chans[i].strip] = i;
    }
  }

  // ── Pass 1: per-side peak anchors, matching the notebook cell 7 ──
  // Even strips (s % 2 == 0, R=LONG): shoulder method — find the short
  //   side (L) peak in events where R reads low (~300 ADC ± 80).
  // Odd strips (s % 2 == 1, L=LONG): 2D correlation ridge — build a 2D
  //   histogram of (R, L), find the most probable R value given the L
  //   beam peak, then extract peak_R.
  // Strips without measured anchors fall back to the median.
  const Int_t kNCorrBins = 256;
  const Double_t kCorrRange = 2000.0;
  const Int_t kCorrMinCounts = 20;

  Bool_t matched[18] = {kFALSE};
  Double_t short_anchor_adc[18] = {0};
  std::vector<Double_t> anchors_found;
  for (Int_t s = 1; s <= 16; s++) {
    if (idx_l[s] < 0 || idx_r[s] < 0)
      continue;
    Bool_t l_is_long = (LongSide(s) == 'L');
    ChannelCal &c_long = chans[l_is_long ? idx_l[s] : idx_r[s]];
    ChannelCal &c_short = chans[l_is_long ? idx_r[s] : idx_l[s]];
    if (!IsCalibrated(c_long)) {
      std::cerr << "  strip " << s
                << ": long side uncalibrated; skipping L/R gain match"
                << std::endl;
      continue;
    }
    Double_t peak_short = 0.0;

    if (l_is_long) {
      // Odd strip (L=LONG, R=SHORT): 2D correlation ridge for peak_R
      const StripPairSamples &p = pairs[s];
      if (Long64_t(p.l.size()) > 100) {
        TH2F *h2d = new TH2F(Form("h2d_lr_corr_s%d", s),
                             ";R (short) [ADC];L (long) [ADC]", kNCorrBins, 0.0,
                             kCorrRange, kNCorrBins, 0.0, kCorrRange);
        h2d->SetDirectory(nullptr);
        for (Int_t j = 0; j < Int_t(p.l.size()); j++)
          h2d->Fill(Double_t(p.r[j]), Double_t(p.l[j]));

        // Notebook: scan from high R downward, find first R bin whose
        // column (projection onto L) has any bin >= 20 counts.
        Double_t rough_R = 0.0;
        for (Int_t ix = kNCorrBins; ix >= 1; ix--) {
          Double_t col_max = 0.0;
          for (Int_t iy = 1; iy <= kNCorrBins; iy++) {
            Double_t v = h2d->GetBinContent(ix, iy);
            if (v > col_max)
              col_max = v;
          }
          if (col_max >= kCorrMinCounts) {
            rough_R = h2d->GetXaxis()->GetBinCenter(ix);
            break;
          }
        }

        // Find offset_r: R peak in events where L ≈ long peak
        Double_t long_peak = c_long.fit_adc;
        std::vector<Float_t> r_near_long;
        for (Int_t j = 0; j < Int_t(p.l.size()); j++) {
          if (TMath::Abs(Double_t(p.l[j]) - long_peak) < kGmPeakWin)
            r_near_long.push_back(p.r[j]);
        }
        Double_t offset_r = 0.0;
        if (r_near_long.size() > 50)
          offset_r = GmFindPeak(r_near_long, 50.0, 600.0, 0.05);

        if (rough_R > 0) {
          peak_short = GmFindPeak(p.r, rough_R * 0.85, rough_R * 1.15, 0.0);
        } else {
          Double_t lo = (offset_r > 0) ? offset_r + 100.0 : 500.0;
          peak_short = GmFindPeak(p.r, lo, 2000.0, 0.05);
        }
        delete h2d;
      }
      if (peak_short <= 0) {
        // Fall back to shoulder if 2D correlation fails
        const std::vector<Float_t> &slice = pairs[s].shoulder;
        if (Long64_t(slice.size()) >= kGmMinSlice)
          peak_short = GmFindPeak(slice, kGmShortLo, kGmShortHi, kGmSkipFrac);
      }
    } else {
      // Even strip (L=SHORT, R=LONG): shoulder method (notebook cell 7)
      const std::vector<Float_t> &slice = pairs[s].shoulder;
      if (Long64_t(slice.size()) >= kGmMinSlice)
        peak_short = GmFindPeak(slice, kGmShortLo, kGmShortHi, kGmSkipFrac);
    }

    if (peak_short <= 0)
      continue;
    short_anchor_adc[s] = peak_short;
    anchors_found.push_back(peak_short);
    std::cout << "  strip " << s << " short_anchor=" << Form("%.1f", peak_short)
              << " ADC (method=" << (l_is_long ? "2Dcorr" : "shoulder") << ")"
              << std::endl;
  }

  Double_t median_anchor = 0.0;
  if (!anchors_found.empty()) {
    std::sort(anchors_found.begin(), anchors_found.end());
    Int_t ms = Int_t(anchors_found.size());
    median_anchor =
        (ms % 2 == 1)
            ? anchors_found[ms / 2]
            : 0.5 * (anchors_found[ms / 2 - 1] + anchors_found[ms / 2]);
  }

  for (Int_t s = 1; s <= 16; s++) {
    if (idx_l[s] < 0 || idx_r[s] < 0)
      continue;
    Bool_t l_is_long = (LongSide(s) == 'L');
    ChannelCal &c_long = chans[l_is_long ? idx_l[s] : idx_r[s]];
    ChannelCal &c_short = chans[l_is_long ? idx_r[s] : idx_l[s]];
    if (!IsCalibrated(c_long))
      continue;
    if (short_anchor_adc[s] <= 0) {
      if (median_anchor <= 0) {
        std::cerr << "  strip " << s
                  << ": no anchor and no median fallback; keeping "
                     "independent gains"
                  << std::endl;
        continue;
      }
      short_anchor_adc[s] = median_anchor;
      std::cout << "  strip " << s << " short_anchor=median fallback "
                << Form("%.1f", median_anchor) << " ADC" << std::endl;
    }
    c_long.gain = 1.0 / c_long.fit_adc;
    c_short.gain = 1.0 / short_anchor_adc[s];
    matched[s] = kTRUE;
    std::cout << "  strip " << s
              << " L/R match: long peak=" << Form("%.1f", c_long.fit_adc)
              << " ADC  short anchor=" << Form("%.1f", short_anchor_adc[s])
              << " ADC" << std::endl;
  }

  // ── Pass 2: per-strip eSum peak, median reference, short-side correction ──
  Double_t esum_peak[18] = {0};
  std::vector<Double_t> peaks_for_median;
  for (Int_t s = 1; s <= 16; s++) {
    if (!matched[s])
      continue;
    const StripPairSamples &p = pairs[s];
    Double_t gl = Gain(chans[idx_l[s]]);
    Double_t gr = Gain(chans[idx_r[s]]);
    std::vector<Float_t> esum;
    esum.reserve(p.l.size());
    for (Int_t j = 0; j < Int_t(p.l.size()); j++)
      esum.push_back(Float_t(gl * Double_t(p.l[j]) + gr * Double_t(p.r[j])));
    esum_peak[s] = GmFindPeak(esum, kGmEsumLo, kGmEsumHi, 0.0);
    if (esum_peak[s] > 0)
      peaks_for_median.push_back(esum_peak[s]);
  }
  if (peaks_for_median.size() < 4) {
    std::cerr << "  L/R gain match pass 2: too few strips with eSum peaks ("
              << peaks_for_median.size() << "); skipping short-side correction"
              << std::endl;
  } else {
    std::sort(peaks_for_median.begin(), peaks_for_median.end());
    Int_t m = Int_t(peaks_for_median.size());
    Double_t esum_ref =
        (m % 2 == 1)
            ? peaks_for_median[m / 2]
            : 0.5 * (peaks_for_median[m / 2 - 1] + peaks_for_median[m / 2]);
    std::cout << "  L/R gain match: eSum reference (median) = "
              << Form("%.4f", esum_ref) << " a.u." << std::endl;

    for (Int_t s = 1; s <= 16; s++) {
      if (!matched[s] || esum_peak[s] <= 0)
        continue;
      Double_t delta = esum_peak[s] - esum_ref;
      Double_t short_now = esum_peak[s] - 1.0;
      Double_t short_want = esum_ref - 1.0;
      Bool_t l_is_long = (LongSide(s) == 'L');
      ChannelCal &c_short = chans[l_is_long ? idx_r[s] : idx_l[s]];
      if (TMath::Abs(delta) > kGmCorrThresh && short_now > 0 &&
          short_want > 0) {
        Double_t corr = short_want / short_now;
        c_short.gain *= corr;
        std::cout << "  strip " << s
                  << " eSum peak=" << Form("%.4f", esum_peak[s])
                  << " delta=" << Form("%+.4f", delta) << " -> short gain x"
                  << Form("%.4f", corr) << std::endl;
      } else {
        std::cout << "  strip " << s
                  << " eSum peak=" << Form("%.4f", esum_peak[s])
                  << " delta=" << Form("%+.4f", delta) << " OK" << std::endl;
      }
    }
  }
}

// Cathode uses median + IQR/1.349 (asymmetric tail not as clean and the user
// prefers to keep cathode on the existing approach). All other channels
// (S guard strips + L/R long anodes) use a robust mode-seeded Gaussian fit
// over the peak core (mode ± 2σ), anchored to the fitted centroid. Fall back
// to the robust mode itself on fit failure.
void ReduceToAnchors(std::vector<ChannelCal> &chans,
                     std::vector<std::vector<Float_t>> &samples,
                     std::vector<TF1 *> &fits_out, const TString &run_label,
                     const StripPairSamples pairs[18]) {
  Int_t n_chans = Int_t(chans.size());
  fits_out.assign(n_chans, nullptr);

  for (Int_t i = 0; i < n_chans; i++) {
    ChannelCal &c = chans[i];
    std::vector<Float_t> &v = samples[i];
    c.n_samples = Long64_t(v.size());

    if (Long64_t(v.size()) < kMinSamples) {
      c.fit_adc = 0;
      c.fit_sigma_adc = 0;
      continue;
    }
    if (c.side == 'C') {
      // Cathode: median + IQR (asymmetric tail, no clean peak).
      c.fit_adc = Median(v);
      c.fit_sigma_adc = InterquartileRange(v) / 1.349;
    } else {
      // Primary: robust mode-seeded Gaussian fit of the peak core. Last
      // resort: the robust mode itself (peak-like), never the tail-biased
      // sample mean.
      Double_t peak = 0, sig = 0;
      TF1 *fit = nullptr;
      TString fname =
          Form("f_peak_gaus_%s_%s", c.name.Data(), run_label.Data());
      if (FitBeamPeakGaussian(v, fname, peak, sig, fit)) {
        c.fit_adc = peak;
        c.fit_sigma_adc = sig;
        fits_out[i] = fit;
      } else {
        // Fit failed; anchor on the robust mode. Still "calibrated", but
        // no fit curve is drawn -- flag it, tagged long/short, since a
        // long-side fallback is a real miscalibration risk.
        Double_t mode = 0.0, rsigma = 0.0;
        RobustPeakSeed(v, mode, rsigma);
        c.fit_adc = mode;
        c.fit_sigma_adc = rsigma;
        TString kind = (c.side == 'S')                 ? "guard"
                       : (c.side == LongSide(c.strip)) ? "long"
                                                       : "short";
        std::cerr << "  [fit-fallback " << kind << "] " << c.name
                  << ": peak fit failed; using mode anchor "
                  << Form("%.1f", c.fit_adc) << " ADC (n=" << c.n_samples << ")"
                  << std::endl;
      }
    }
    std::cout << "  " << c.name << " anchor[ADC]=" << c.fit_adc
              << " sig=" << c.fit_sigma_adc << " (n=" << c.n_samples << ")"
              << std::endl;
  }

  // After all per-channel peaks are fitted, run the check_LR-style two-pass
  // L/R gain matching for strips 1–16: long-side beam peak + short-side
  // shoulder anchors, then a per-strip eSum median alignment applied to the
  // short side only. Sets the ChannelCal::gain overrides; strips where the
  // shoulder cannot be found keep the independent 1/fit_adc gains.
  if (pairs)
    ComputeLRGainMatch(chans, pairs);
}

void WriteEresTomlRaw(const TString &out_subpath,
                      const Double_t eres_vals[35]) {
  toml::table eres_tbl;
  eres_tbl.insert("Cathode", eres_vals[0]);
  eres_tbl.insert("S0", eres_vals[1]);
  eres_tbl.insert("S17", eres_vals[2]);
  for (Int_t s = 1; s <= 16; s++) {
    std::string key = "L" + std::to_string(s);
    eres_tbl.insert(key, eres_vals[3 + (s - 1)]);
  }
  for (Int_t s = 1; s <= 16; s++) {
    std::string key = "R" + std::to_string(s);
    eres_tbl.insert(key, eres_vals[19 + (s - 1)]);
  }
  toml::table detector_tbl;
  detector_tbl.insert("eres", eres_tbl);
  toml::table root_tbl;
  root_tbl.insert("detector", detector_tbl);

  // The eres calibration TOML is a small, version-controlled input (a control
  // file), not bulk output: write it into the repo's control/ dir alongside the
  // other Calibration_Run*_eres.toml, regardless of where root_files point.
  TString out_dir = Paths::DatasetDir() + "/sim_control";
  gSystem->mkdir(out_dir, kTRUE);
  TString out_full = out_dir + "/" + out_subpath;
  std::ofstream f(out_full.Data());
  if (!f) {
    std::cerr << "Cannot write eres TOML: " << out_full << std::endl;
    return;
  }
  f << root_tbl << std::endl;
  std::cout << "  wrote eres TOML: " << out_full << std::endl;
}

// Writes a one-row `calibration` tree into the open file `dst` (typically a
// per-subfile .cal.root). Layout matches AggregateEresTomlForRun's reader.
// `align` carries the beam-energy window and per-strip alignment factors;
// pass nullptr when neither has been computed yet (branches are written as
// zero).
void WriteCalibrationTree(TFile *dst, const std::vector<ChannelCal> &chans,
                          const StripAlignmentResult *align) {
  dst->cd();
  if (TObject *old = dst->Get("calibration"))
    old->Delete();
  TTree *cal = new TTree("calibration", "Per-channel normMUSIC calibration");
  Float_t gain[kMaxChannels] = {0};
  Float_t fit_adc[kMaxChannels] = {0}, fit_sigma[kMaxChannels] = {0};
  Long64_t fit_n[kMaxChannels] = {0};
  Bool_t ok[kMaxChannels] = {0};
  // Per-strip gains laid out to match the events tree exactly: GainLeft[s]
  // multiplies Left_0_17_dE[s] (s=0/17 are the single-ended guards, s=1..16 the
  // left ends), GainRight[s] multiplies RightdE[s] (0 at the guards). This is
  // what EnergyView reads to calibrate on the fly -- no per-event a.u. is
  // stored.
  Float_t gain_left[18] = {0}, gain_right[18] = {0};
  Float_t gain_cathode = 0.0f;
  Int_t n_actual = TMath::Min(Int_t(chans.size()), kMaxChannels);
  for (Int_t k = 0; k < n_actual; k++) {
    const ChannelCal &c = chans[k];
    ok[k] = IsCalibrated(c) || c.gain > 0;
    gain[k] = ok[k] ? Float_t(Gain(c)) : 0.0f;
    fit_adc[k] = Float_t(c.fit_adc);
    fit_sigma[k] = Float_t(c.fit_sigma_adc);
    fit_n[k] = c.n_samples;
    if (c.side == 'S' && c.strip >= 0 && c.strip <= 17)
      gain_left[c.strip] = gain[k];
    else if (c.side == 'L' && c.strip >= 1 && c.strip <= 16)
      gain_left[c.strip] = gain[k];
    else if (c.side == 'R' && c.strip >= 1 && c.strip <= 16)
      gain_right[c.strip] = gain[k];
    else if (c.side == 'C')
      gain_cathode = gain[k];
  }
  cal->Branch("Gain", gain, Form("Gain[%d]/F", kMaxChannels));
  cal->Branch("Ok", ok, Form("Ok[%d]/O", kMaxChannels));
  cal->Branch("FitADC", fit_adc, Form("FitADC[%d]/F", kMaxChannels));
  cal->Branch("FitSigmaADC", fit_sigma,
              Form("FitSigmaADC[%d]/F", kMaxChannels));
  cal->Branch("FitN", fit_n, Form("FitN[%d]/L", kMaxChannels));
  cal->Branch("GainLeft", gain_left, "GainLeft[18]/F");
  cal->Branch("GainRight", gain_right, "GainRight[18]/F");
  cal->Branch("GainCathode", &gain_cathode, "GainCathode/F");

  // Beam-energy window and per-strip multiplicative alignment factors
  // matching the notebook approach (pol3 reference / centroid).
  // EnergyView applies total_corrected = factor * total after the
  // per-channel gain. Default factor = 1.0 (identity).
  Float_t beam_e_min = 0.0f, beam_e_max = 0.0f;
  Float_t strip_factor[18];
  for (Int_t s = 0; s < 18; s++)
    strip_factor[s] = 1.0f;
  if (align) {
    beam_e_min = Float_t(align->beam_e_min);
    beam_e_max = Float_t(align->beam_e_max);
    if (align->ok) {
      for (Int_t s = 0; s < 18; s++)
        strip_factor[s] = Float_t(align->factors[s]);
    }
  }
  cal->Branch("BeamEMin", &beam_e_min, "BeamEMin/F");
  cal->Branch("BeamEMax", &beam_e_max, "BeamEMax/F");
  cal->Branch("StripFactor", strip_factor, "StripFactor[18]/F");
  cal->Fill();
  cal->Write("calibration", TObject::kOverwrite);
}

// Per-channel ADC histogram of the samples that fed each beam anchor. One file
// per channel under <plot_subdir>/beam_peak, named beam_peak_<channel>.
void SaveBeamPeakChannelHistograms(
    const std::vector<ChannelCal> &chans,
    const std::vector<std::vector<Float_t>> &samples,
    const std::vector<TF1 *> &fits, const TString &plot_subdir) {
  TString subdir = plot_subdir + "/beam_peak";
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    const ChannelCal &c = chans[i];
    const std::vector<Float_t> &v = samples[i];
    if (Long64_t(v.size()) < kMinSamples)
      continue;
    // Same binning recipe the fit used, so the overlaid Gaussian's amplitude
    // matches this histogram exactly.
    Double_t mode = 0.0, sigma = 0.0;
    RobustPeakSeed(v, mode, sigma);
    TH1F *h = MakeBeamPeakHist(Form("h_beam_peak_%s", c.name.Data()),
                               Form(";%s #DeltaE [ADC];Counts", c.name.Data()),
                               v, mode, sigma);
    TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDrawHistogram(h, kBlack);
    TF1 *fit = fits[i];
    if (fit) {
      fit->SetLineColor(kViolet + 2);
      fit->SetLineWidth(2);
      fit->Draw("L SAME");
    }
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(cv, Form("beam_peak_%s", c.name.Data()), subdir,
                                PlotSaveOptions::kLINEAR);
    delete cv;
    delete h;
  }
}

// Writes the per-channel gain table (tree "calibration") into the subfile's own
// events file. No per-event calibrated tree is produced: downstream readers
// recover a.u. on the fly via gain x raw ADC (EnergyView), so the raw events
// tree plus this one-row gain table fully determine every calibrated value.
void WriteCalibrationToEvents(const FileSpec &spec,
                              const std::vector<ChannelCal> &chans,
                              const StripAlignmentResult *align) {
  TString events_subpath = FileSet::EventsName(spec) + ".root";
  TFile *f = IO::OpenForWriting(events_subpath, "UPDATE");
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open " << events_subpath << " to write calibration"
              << std::endl;
    if (f)
      delete f;
    return;
  }
  WriteCalibrationTree(f, chans, align);
  std::cout << "  wrote calibration into " << events_subpath << std::endl;
  f->Close();
  delete f;
}

// Per-channel calibrated overlay for one subfile, via AttachCalSidecar (the
// same path downstream macros use). The sidecar must already be on disk.
void SaveDynamicRangeOverlay(const FileSpec &spec,
                             const std::vector<ChannelCal> &chans,
                             const TString &plot_subdir,
                             const TString &file_label) {
  const Int_t kNStrips = 18;
  const Int_t nbins = 300;
  const Double_t emin = Constants::cfg.STRIP_DE_MIN_NORMED;
  const Double_t emax = Constants::cfg.STRIP_DE_MAX_NORMED;
  TH1D *h[kNStrips];
  for (Int_t s = 0; s < kNStrips; s++) {
    h[s] = new TH1D(Form("h_dynrange_%s_S%d", file_label.Data(), s),
                    ";#DeltaE [a.u.];Counts", nbins, emin, emax);
    h[s]->SetDirectory(nullptr);
  }
  TString sub = FileSet::EventsName(spec) + ".root";
  TFile *sf = IO::OpenForReading(sub);
  if (!sf || sf->IsZombie()) {
    if (sf)
      sf->Close();
    for (Int_t s = 0; s < kNStrips; s++)
      delete h[s];
    return;
  }
  TTree *tree = static_cast<TTree *>(sf->Get("events"));
  if (!tree) {
    sf->Close();
    delete sf;
    for (Int_t s = 0; s < kNStrips; s++)
      delete h[s];
    return;
  }
  EnergyView ev;
  ev.Attach(tree);
  if (!ev.is_normed) {
    sf->Close();
    delete sf;
    for (Int_t s = 0; s < kNStrips; s++)
      delete h[s];
    return;
  }
  Long64_t n = tree->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
    tree->GetEntry(j);
    ev.Decode();
    for (Int_t s = 0; s < kNStrips; s++) {
      Double_t v = ev.total[s];
      if (v > 0)
        h[s]->Fill(v);
    }
  }
  sf->Close();
  delete sf;
  std::vector<Int_t> colors = PlottingUtils::GetDefaultColors();
  Double_t y_top = 0;
  for (Int_t s = 0; s < kNStrips; s++) {
    Double_t m = h[s]->GetMaximum();
    if (m > y_top)
      y_top = m;
  }
  TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
  cv->SetRightMargin(0.20);
  Bool_t first = kTRUE;
  for (Int_t s = 0; s < kNStrips; s++) {
    Int_t color = colors[s % Int_t(colors.size())];
    h[s]->SetLineColor(color);
    h[s]->SetLineWidth(2);
    h[s]->SetMaximum(1.15 * y_top);
    h[s]->Draw(first ? "HIST" : "HIST SAME");
    first = kFALSE;
  }
  TLegend *leg = PlottingUtils::AddLegend(0.81, 0.99, 0.10, 0.95);
  for (Int_t s = 0; s < kNStrips; s++)
    leg->AddEntry(h[s], Form("S%d", s), "l");
  leg->Draw();
  if (Constants::cfg.SAVE_PLOTS)
    PlottingUtils::SaveFigure(cv, "dynamic_range_check", plot_subdir,
                              PlotSaveOptions::kLOG);
  delete cv;
  delete leg;
  for (Int_t s = 0; s < kNStrips; s++)
    delete h[s];
}

// Overlay (one color per channel, log-y) of ONLY the events used for
// calibration: the beam anchor samples, converted to a.u. via each channel's
// gain. Same axes/style as SaveDynamicRangeOverlay but restricted to
// calibration events rather than the full spectrum.
void CalibrateBeam::SaveCalibSampleOverlay(
    const std::vector<ChannelCal> &chans,
    const std::vector<std::vector<Float_t>> &samples,
    const TString &plot_subdir, const TString &file_label) {
  const Int_t n_chans = Int_t(chans.size());
  const Int_t nbins = 300;
  const Double_t emin = Constants::cfg.STRIP_DE_MIN_NORMED;
  const Double_t emax = Constants::cfg.STRIP_DE_MAX_NORMED;
  std::vector<TH1D *> h(n_chans, nullptr);
  for (Int_t i = 0; i < n_chans; i++) {
    const ChannelCal &c = chans[i];
    if (!IsBeamdEChannel(c) || !IsCalibrated(c))
      continue;
    TString hname =
        Form("h_calibrange_%s_%s", file_label.Data(), c.name.Data());
    h[i] = new TH1D(hname, ";#DeltaE [a.u.];Counts", nbins, emin, emax);
    h[i]->SetDirectory(nullptr);
    const std::vector<Float_t> &v = samples[i];
    for (Int_t j = 0; j < Int_t(v.size()); j++) {
      Double_t val = ApplyCal(c, Double_t(v[j]));
      if (val > 0)
        h[i]->Fill(val);
    }
  }

  std::vector<Int_t> colors = PlottingUtils::GetDefaultColors();
  Double_t y_top = 0;
  for (Int_t i = 0; i < n_chans; i++) {
    if (!h[i])
      continue;
    Double_t m = h[i]->GetMaximum();
    if (m > y_top)
      y_top = m;
  }
  TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
  cv->SetRightMargin(0.20);
  Bool_t first = kTRUE;
  for (Int_t i = 0; i < n_chans; i++) {
    if (!h[i])
      continue;
    Int_t color = colors[i % Int_t(colors.size())];
    h[i]->SetLineColor(color);
    h[i]->SetLineWidth(2);
    h[i]->SetMaximum(1.15 * y_top);
    h[i]->Draw(first ? "HIST" : "HIST SAME");
    first = kFALSE;
  }
  TLegend *leg = PlottingUtils::AddLegend(0.81, 0.99, 0.10, 0.95);
  for (Int_t i = 0; i < n_chans; i++) {
    if (!h[i])
      continue;
    leg->AddEntry(h[i], chans[i].name.Data(), "l");
  }
  leg->Draw();
  if (Constants::cfg.SAVE_PLOTS)
    PlottingUtils::SaveFigure(cv, "dynamic_range_calib_events", plot_subdir,
                              PlotSaveOptions::kLOG);
  delete cv;
  delete leg;
  for (Int_t i = 0; i < n_chans; i++)
    delete h[i];
}

// Derives the beam-energy window from the Strip0 (entrance guard) beam-peak
// fit. The notebook (37Cl_an.ipynb cell 10) fits a Gaussian to raw stp0 ADC
// and takes mu ± 3*sigma. Here, Strip0's beam peak is already fitted in
// ReduceToAnchors (fit_adc / fit_sigma_adc), so we convert to a.u. via the
// channel's own gain. In a.u. the peak sits at 1.0 by construction, so the
// window is 1.0 ± 3*sigma/fit_adc.
void DeriveBeamEnergyWindow(const std::vector<ChannelCal> &chans,
                            StripAlignmentResult &align) {
  const Double_t kBeamNSigma = 3.0;
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    const ChannelCal &c = chans[i];
    if (c.side == 'S' && c.strip == 0 && IsCalibrated(c)) {
      Double_t g = Gain(c);
      align.beam_e_min = g * (c.fit_adc - kBeamNSigma * c.fit_sigma_adc);
      align.beam_e_max = g * (c.fit_adc + kBeamNSigma * c.fit_sigma_adc);
      std::cout << "  beam energy window (from Strip0): [" << align.beam_e_min
                << ", " << align.beam_e_max << "] a.u." << std::endl;
      return;
    }
  }
  std::cerr << "  beam energy window: Strip0 not calibrated, using [0, 0]"
            << std::endl;
}

// Per-strip multiplicative alignment, matching the notebook
// (37Cl_an.ipynb cell 5). Decodes events with the per-channel gains already on
// disk, finds each strip's beam-peak centroid from the eSum 2D histogram,
// fits a robust degree-3 polynomial reference trend through the centroids
// (strips 1-16 only — the single-ended anodes sit at a different scale), and
// derives a multiplicative factor = reference[s] / centroid[s] for 1-16.
// Strips 0/17 get factor = 1.0 / centroid (push beam peak to 1.0 a.u.).
//
// Uses ALL events (not beam-gated): the beam dominates every strip's
// histogram by a wide margin.
StripAlignmentResult FindStripCentroidAlignment(const FileSpec &spec,
                                                const TString &plot_subdir,
                                                const TString &file_label) {
  const Int_t kNStrips = 18;
  const Int_t kNHistBins = 200;
  const Double_t kHistMin = 0.0;
  const Double_t kHistMax = 10.0;
  const Int_t kSmoothTimes = 5;
  const Double_t kMisalignPct = 1.5;
  const Int_t kMaxIter = 4;
  const Int_t kPolyDeg = 3;
  const Double_t kGausFitHalfWidth = 0.15;

  StripAlignmentResult result;
  for (Int_t s = 0; s < kNStrips; s++) {
    result.factors[s] = 1.0;
    result.centroids[s] = 0.0;
  }

  TString sub = FileSet::EventsName(spec) + ".root";
  TFile *sf = IO::OpenForReading(sub);
  if (!sf || sf->IsZombie()) {
    if (sf)
      sf->Close();
    delete sf;
    return result;
  }
  TTree *tree = static_cast<TTree *>(sf->Get("events"));
  if (!tree) {
    sf->Close();
    delete sf;
    return result;
  }
  EnergyView ev;
  ev.Attach(tree);
  if (!ev.is_normed) {
    std::cerr << "  " << file_label
              << ": calibration tree not found -- skipping alignment"
              << std::endl;
    sf->Close();
    delete sf;
    return result;
  }

  TH2D *h2 = new TH2D(Form("h2_strip_align_%s", file_label.Data()),
                      ";Strip number;#DeltaE [a.u.]", kNStrips, -0.5,
                      kNStrips - 0.5, kNHistBins, kHistMin, kHistMax);
  h2->SetDirectory(nullptr);

  Long64_t n = tree->GetEntries();
  Long64_t n_used = 0;
  for (Long64_t j = 0; j < n; j++) {
    tree->GetEntry(j);
    ev.Decode();
    Bool_t any = kFALSE;
    for (Int_t s = 0; s < kNStrips; s++) {
      Double_t v = ev.total[s];
      if (v <= 0)
        continue;
      h2->Fill(Double_t(s), v);
      any = kTRUE;
    }
    if (any)
      n_used++;
  }
  sf->Close();
  delete sf;

  std::cout << "  strip alignment: " << n_used << " events decoded"
            << std::endl;

  Double_t beam_centroids[kNStrips] = {0};
  Bool_t beam_ok[kNStrips] = {kFALSE};

  for (Int_t s = 0; s < kNStrips; s++) {
    Int_t bin_ix = s + 1;
    TH1D *proj = h2->ProjectionY(
        Form("hproj_align_%s_s%d", file_label.Data(), s), bin_ix, bin_ix);
    proj->SetDirectory(nullptr);
    Long64_t n_entries = Long64_t(proj->GetEntries());
    if (n_entries < kMinSamples) {
      std::cerr << "  strip " << s << ": too few entries for alignment ("
                << n_entries << ")" << std::endl;
      delete proj;
      continue;
    }
    proj->Smooth(kSmoothTimes);

    // Strip 0/17: single-ended anode, beam at ~1.0, pileup at ~2.0.
    // Strips 1-16: total is bimodal when the short side doesn't fire
    // (long-only ≈ 0.5) — pick the peak nearest 1.0, which is the
    // full-strip beam peak (both sides contributing, eSum ≈ 1.0).
    Double_t peak_target = 1.0;
    Double_t search_lo, search_hi;
    if (s == 0 || s == 17) {
      search_lo = 0.3;
      search_hi = 1.6;
    } else {
      search_lo = 0.6;
      search_hi = 1.5;
    }
    Int_t b_lo = proj->FindBin(search_lo);
    Int_t b_hi = proj->FindBin(search_hi);
    Int_t b_max = -1;
    Double_t best_dist = 1e9;
    for (Int_t b = b_lo; b <= b_hi; b++) {
      Double_t v = proj->GetBinContent(b);
      if (v <= 0)
        continue;
      // Local maximum check: higher than neighbours
      if (b > b_lo && proj->GetBinContent(b - 1) >= v)
        continue;
      if (b < b_hi && proj->GetBinContent(b + 1) > v)
        continue;
      Double_t bc = proj->GetBinCenter(b);
      Double_t d = TMath::Abs(bc - peak_target);
      if (d < best_dist) {
        best_dist = d;
        b_max = b;
      }
    }
    // Fallback: global maximum if no local peaks found near 1.0
    if (b_max < 0) {
      Double_t val_max = 0;
      for (Int_t b = b_lo; b <= b_hi; b++) {
        Double_t v = proj->GetBinContent(b);
        if (v > val_max) {
          val_max = v;
          b_max = b;
        }
      }
    }
    if (b_max < 0) {
      delete proj;
      continue;
    }
    Double_t seed_peak = proj->GetBinCenter(b_max);
    Double_t seed_h = proj->GetBinContent(b_max);
    if (seed_h <= 0 || seed_peak <= 0) {
      delete proj;
      continue;
    }

    // Sub-bin precision: Gaussian fit around the smoothed max-bin seed
    Double_t fit_lo = seed_peak - kGausFitHalfWidth;
    Double_t fit_hi = seed_peak + kGausFitHalfWidth;
    if (fit_lo < kHistMin)
      fit_lo = kHistMin;
    if (fit_hi > kHistMax)
      fit_hi = kHistMax;
    TF1 *fg = new TF1("f_align_peak_refine", "gaus", fit_lo, fit_hi);
    fg->SetParameters(seed_h, seed_peak, 0.05);
    fg->SetParLimits(1, fit_lo, fit_hi);
    TFitResultPtr r = proj->Fit(fg, "QSRN");
    Double_t beam_peak = fg->GetParameter(1);
    if (!(r.Get() && r->IsValid() && beam_peak > 0))
      beam_peak = seed_peak;
    delete fg;

    beam_centroids[s] = beam_peak;
    beam_ok[s] = kTRUE;
    std::cout << "  strip " << s << " beam=" << Form("%.4f", beam_peak)
              << " a.u.  (n=" << n_entries << ")" << std::endl;
    delete proj;
  }

  // Robust pol3 reference trend through strips 1-16 centroids.
  // Iteratively drop the strip with the worst residual > kMisalignPct.
  TGraph *g_cent = new TGraph(kNStrips);
  Int_t np = 0;
  for (Int_t s = 1; s <= 16; s++) {
    if (beam_ok[s]) {
      g_cent->SetPoint(np, Double_t(s), beam_centroids[s]);
      np++;
    }
  }
  g_cent->Set(np);

  std::set<Int_t> outliers;
  for (Int_t iter = 0; iter < kMaxIter; iter++) {
    if (g_cent->GetN() <= kPolyDeg + 1)
      break;
    TF1 *fpol = new TF1(Form("fpol_align_%s_iter%d", file_label.Data(), iter),
                        "pol3", -0.5, kNStrips - 0.5);
    TFitResultPtr r = g_cent->Fit(fpol, "QSRN");
    if (!r.Get() || !r->IsValid()) {
      delete fpol;
      break;
    }
    Double_t worst_pct = 0;
    Int_t worst_idx = -1;
    for (Int_t i = 0; i < g_cent->GetN(); i++) {
      Double_t x, y;
      g_cent->GetPoint(i, x, y);
      Double_t pred = fpol->Eval(x);
      if (pred <= 0)
        continue;
      Double_t resid_pct = TMath::Abs(y - pred) / pred * 100.0;
      if (resid_pct > worst_pct) {
        worst_pct = resid_pct;
        worst_idx = i;
      }
    }
    delete fpol;
    if (worst_pct < kMisalignPct || worst_idx < 0)
      break;
    Double_t rx, ry;
    g_cent->GetPoint(worst_idx, rx, ry);
    std::cout << "  strip " << Int_t(rx) << " outlier "
              << Form("%.1f", worst_pct) << "% — dropped" << std::endl;
    outliers.insert(Int_t(rx));
    g_cent->RemovePoint(worst_idx);
  }

  TF1 *fbeam = nullptr;
  if (g_cent->GetN() > kPolyDeg + 1) {
    fbeam = new TF1(Form("fbeam_align_%s", file_label.Data()), "pol3", -0.5,
                    kNStrips - 0.5);
    g_cent->Fit(fbeam, "QSRN");
    std::cout << "  strip alignment reference: pol3 fitted through "
              << g_cent->GetN() << " strips" << std::endl;
  }

  // Notebook (37Cl_an.ipynb cell 5): factors = reference / centroid,
  // where reference is the pol3 trend through strips 1-16 centroids.
  // Strips 0/17 are single-ended anodes not in the pol3 fit; push to 1.0.
  Int_t valid_strips = 0;
  for (Int_t s = 0; s < kNStrips; s++) {
    if (!beam_ok[s])
      continue;
    Double_t centro = beam_centroids[s];
    if (centro <= 0)
      continue;
    result.centroids[s] = centro;
    if (s >= 1 && s <= 16 && fbeam) {
      Double_t ref = fbeam->Eval(Double_t(s));
      if (ref > 0)
        result.factors[s] = ref / centro;
      else
        result.factors[s] = 1.0 / centro;
    } else {
      result.factors[s] = 1.0 / centro;
    }
    valid_strips++;
    std::cout << "  strip " << s << " centroid=" << Form("%.4f", centro)
              << " factor=" << Form("%.4f", result.factors[s]) << std::endl;
  }
  result.ok = (valid_strips >= 4) ? kTRUE : kFALSE;

  // Diagnostic plot
  {
    TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::ConfigureAndDraw2DHistogram(h2, cv);
    TGraph *g_beam_plot = new TGraph(kNStrips);
    Int_t nb = 0;
    for (Int_t s = 0; s < kNStrips; s++) {
      if (beam_ok[s]) {
        g_beam_plot->SetPoint(nb, Double_t(s), beam_centroids[s]);
        nb++;
      }
    }
    g_beam_plot->Set(nb);
    if (nb > 0) {
      g_beam_plot->SetMarkerStyle(20);
      g_beam_plot->SetMarkerColor(kOrange);
      g_beam_plot->Draw("P SAME");
    }
    if (fbeam) {
      fbeam->SetLineColor(kViolet + 2);
      fbeam->SetLineWidth(2);
      fbeam->Draw("SAME");
    }
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(cv, "strip_alignment_check", plot_subdir,
                                PlotSaveOptions::kLINEAR);
    delete cv;
    delete g_beam_plot;
  }

  delete fbeam;
  delete g_cent;
  delete h2;

  return result;
}

void CalibrateBeam::CalibrateBeamOneSubfile(
    const FileSpec &spec, const std::vector<ChannelCal> &chans_template) {
  TString file_label = FileSet::FileLabel(spec);
  TString plot_subdir = "beam_calibration/" + file_label;
  std::cout << "Beam calibration: " << file_label << std::endl;

  BeamFit2D beam;
  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    beam = FindBeamGateStp2VsStp1(spec, file_label, plot_subdir);
  }
  if (!beam.ok) {
    std::cerr << "  " << file_label << ": Strip2-vs-Strip1 beam gate failed"
              << std::endl;
    return;
  }

  BeamFit2D beam0vGrid;
  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    beam0vGrid = FindBeamGateStp0VsGrid(spec, file_label, plot_subdir, beam);
  }
  if (!beam0vGrid.ok)
    std::cerr << "  " << file_label
              << ": Strip0-vs-Grid beam gate failed — continuing without it"
              << std::endl;

  std::vector<ChannelCal> chans = chans_template;
  std::vector<std::vector<Float_t>> samples;
  StripPairSamples pairs[18];
  CollectAnchorSamplesOneSubfile(spec, chans, beam, beam0vGrid, samples, pairs);
  std::vector<TF1 *> peak_fits;
  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    ReduceToAnchors(chans, samples, peak_fits, file_label, pairs);
    SaveBeamPeakChannelHistograms(chans, samples, peak_fits, plot_subdir);
  }
  for (Int_t i = 0; i < Int_t(peak_fits.size()); i++)
    delete peak_fits[i];
  peak_fits.clear();

  // Every channel that failed to calibrate is silently forced to gain 0 (reads
  // 0 a.u. and drops out of the strip total), so spell out which ones and why.
  // The long/short tag makes a long-side failure -- the dominant signal, which
  // should never starve -- easy to spot: grep "[uncalibrated long]".
  for (Int_t i = 0; i < Int_t(chans.size()); i++) {
    const ChannelCal &c = chans[i];
    if (IsCalibrated(c))
      continue;
    TString kind;
    if (c.side == 'C')
      kind = "cathode";
    else if (c.side == 'S')
      kind = "guard";
    else
      kind = (c.side == LongSide(c.strip)) ? "long" : "short";
    TString why;
    if (c.n_samples < kMinSamples)
      why =
          Form("too few beam samples (%lld < %lld)", c.n_samples, kMinSamples);
    else
      why = Form("bad exp anchor (fit_adc=%.1f)", c.fit_adc);
    std::cerr << "  [uncalibrated " << kind << "] " << c.name << " -> gain 0; "
              << why << " (beam n=" << c.n_samples << ")" << std::endl;
  }

  for (Int_t i = 0; i < Int_t(chans.size()); i++)
    if (IsCalibrated(chans[i]))
      std::cout << "  " << chans[i].name << " gain=" << Gain(chans[i])
                << " a.u./ADC  resolution="
                << Form("%.2f", ResolutionFWHMPercent(chans[i])) << "% FWHM"
                << std::endl;

  // Derive beam energy window from Strip0 (mu ± 3*sigma in a.u.).
  StripAlignmentResult align;
  DeriveBeamEnergyWindow(chans, align);

  // Write initial calibration tree: per-channel gains + beam window. The
  // alignment step needs this on disk so EnergyView can decode events in a.u.
  WriteCalibrationToEvents(spec, chans, &align);

  // Per-strip multiplicative alignment: decode events with the per-channel
  // gains just written, find each strip's beam-peak centroid, fit a pol3
  // reference trend, and derive factors = reference / centroid. Stored in
  // the calibration tree and applied by EnergyView as total *= factor.
  // gains just written, find each strip's beam-peak AND pileup-peak centroids,
  // and derive a linear correction (slope + intercept) that maps beam→1.0 and
  // pileup→2.0, flattening the Bragg curve. The polynomial trend is used only
  // Preserve beam energy window from DeriveBeamEnergyWindow across the
  // alignment call (which returns a fresh StripAlignmentResult).
  Double_t beam_e_min = align.beam_e_min;
  Double_t beam_e_max = align.beam_e_max;
  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    align = FindStripCentroidAlignment(spec, plot_subdir, file_label);
  }
  align.beam_e_min = beam_e_min;
  align.beam_e_max = beam_e_max;
  if (align.ok) {
    WriteCalibrationToEvents(spec, chans, &align);
  }

  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    SaveDynamicRangeOverlay(spec, chans, plot_subdir, file_label);
  }
  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    SaveCalibSampleOverlay(chans, samples, plot_subdir, file_label);
  }

  std::cout << "  " << file_label << " calibration complete." << std::endl;
}

void CalibrateBeam::AggregateEresTomlForRun(
    Int_t run, const std::vector<FileSpec> &specs) {
  const Int_t n_eres = 35;
  std::vector<std::vector<Double_t>> fwhm_per_chan(n_eres);

  std::vector<ChannelCal> tmpl = CalibrateBeam::BuildChannels();
  Int_t n_chans = Int_t(tmpl.size());

  for (Int_t s = 0; s < Int_t(specs.size()); s++) {
    // The calibration tree now lives inside each subfile's events file.
    TString cal_sub = FileSet::EventsName(specs[s]) + ".root";
    TFile *cf = IO::OpenForReading(cal_sub);
    if (!cf || cf->IsZombie()) {
      if (cf)
        delete cf;
      continue;
    }
    TTree *t = static_cast<TTree *>(cf->Get("calibration"));
    if (!t) {
      cf->Close();
      delete cf;
      continue;
    }
    Float_t fit_adc[kMaxChannels] = {0};
    Float_t fit_sigma[kMaxChannels] = {0};
    Bool_t ok[kMaxChannels] = {0};
    t->SetBranchAddress("FitADC", fit_adc);
    t->SetBranchAddress("Ok", ok);
    t->SetBranchAddress("FitSigmaADC", fit_sigma);
    if (t->GetEntries() < 1) {
      cf->Close();
      delete cf;
      continue;
    }
    t->GetEntry(0);

    for (Int_t i = 0; i < n_chans && i < kMaxChannels; i++) {
      if (!ok[i])
        continue;
      Double_t sig_adc = fit_sigma[i];
      if (sig_adc <= 0 || fit_adc[i] <= 0)
        continue;
      // Relative resolution in % FWHM, straight from the raw-ADC peak fit —
      // independent of the normMUSIC gain/normalization by construction.
      const Double_t kFwhmPerSigma = 2.0 * TMath::Sqrt(2.0 * TMath::Log(2.0));
      Double_t fwhm_pct =
          100.0 * kFwhmPerSigma * Double_t(sig_adc) / Double_t(fit_adc[i]);
      Int_t idx = ChannelToEresIndex(tmpl[i]);
      if (idx >= 0 && idx < n_eres)
        fwhm_per_chan[idx].push_back(fwhm_pct);
    }
    cf->Close();
    delete cf;
  }

  Double_t eres_vals[35];
  for (Int_t i = 0; i < n_eres; i++)
    eres_vals[i] = -1.0;
  for (Int_t i = 0; i < n_eres; i++) {
    std::vector<Double_t> &v = fwhm_per_chan[i];
    if (v.empty())
      continue;
    std::sort(v.begin(), v.end());
    Int_t m = Int_t(v.size());
    eres_vals[i] = (m % 2 == 1) ? v[m / 2] : 0.5 * (v[m / 2 - 1] + v[m / 2]);
  }
  std::cout << "Run " << run << ": writing per-channel %FWHM medians"
            << std::endl;
  WriteEresTomlRaw(Form("Calibration_Run%d_eres.toml", run), eres_vals);
}

void CalibrateBeam::Run(const TString &file_label) {
  const TString project_root = Paths::DatasetDir();
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                Paths::ResultsDir() + "/plots",
                                Paths::ResultsDir() + "/root_files");
  gROOT->SetBatch(kTRUE);

  std::vector<FileSpec> specs;
  if (file_label.IsNull()) {
    specs = FileSet::BuildProcessedFileSpecs();
    if (specs.empty()) {
      std::cerr << "No file specs from FileSet::BuildProcessedFileSpecs()"
                << std::endl;
      return;
    }
  } else {
    FileSpec s = FileSet::ResolveFileSpec(file_label);
    if (s.run < 0) {
      std::cerr << "Could not resolve file label '" << file_label << "'"
                << std::endl;
      return;
    }
    specs.push_back(s);
  }

  std::vector<ChannelCal> chans = CalibrateBeam::BuildChannels();

  Int_t n_specs = Int_t(specs.size());

  std::set<Int_t> runs;
  for (Int_t k = 0; k < n_specs; k++)
    runs.insert(specs[k].run);

  Int_t n_workers =
      TMath::Min(Int_t(std::thread::hardware_concurrency()), n_specs);
  n_workers = TMath::Min(n_workers, Constants::cfg.MAX_FUSED_WORKERS);
  if (n_workers < 1)
    n_workers = 1;
  std::cout << "calibrate-beam: " << n_specs << " subfiles on " << n_workers
            << " workers" << std::endl;

  std::queue<Int_t> work;
  for (Int_t k = 0; k < n_specs; k++)
    work.push(k);
  std::mutex work_mutex;

  std::vector<std::thread> workers;
  for (Int_t w = 0; w < n_workers; w++) {
    workers.emplace_back([&]() {
      while (true) {
        Int_t k;
        {
          std::lock_guard<std::mutex> lk(work_mutex);
          if (work.empty())
            return;
          k = work.front();
          work.pop();
        }
        CalibrateBeamOneSubfile(specs[k], chans);
      }
    });
  }
  for (Int_t w = 0; w < Int_t(workers.size()); w++)
    workers[w].join();

  for (std::set<Int_t>::const_iterator it = runs.begin(); it != runs.end();
       ++it) {
    std::vector<FileSpec> run_specs;
    for (Int_t k = 0; k < n_specs; k++)
      if (specs[k].run == *it)
        run_specs.push_back(specs[k]);
    AggregateEresTomlForRun(*it, run_specs);
  }
}
