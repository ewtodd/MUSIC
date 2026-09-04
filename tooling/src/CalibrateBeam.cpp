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
                     kBeamGateNBins, 0.0, Constants::ActiveStripEMaxAdc(),
                     kBeamGateNBins, 0.0, Constants::ActiveStripEMaxAdc());
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
  // Beam-gated, long-side-triggered pairs. Unlike `l`/`r` these keep events
  // where the short end did not fire (recorded as 0), so they describe the
  // population the decode actually sums. `l`/`r` require both ends and are
  // therefore biased towards anomalously large short-side signals.
  std::vector<Float_t> gated_long;
  std::vector<Float_t> gated_short;
};

// check_LR gain-match slice constants (37Cl_an_check_LR.ipynb, cell "Save
// 2-pass calibration"). Defined before CollectAnchorSamplesOneSubfile because
// the shoulder slice is selected during sample collection.

void CollectAnchorSamplesOneSubfile(const FileSpec &spec,
                                    const std::vector<ChannelCal> &chans,
                                    const BeamFit2D &beam,
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
  Short_t cathode_adc = 0;
  tree->SetBranchAddress("Left_0_17_dE", left_0_17_adc);
  tree->SetBranchAddress("RightdE", rightdE_adc);
  tree->SetBranchAddress("Cathode", &cathode_adc);

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
      }
    }
    Double_t x = Double_t(left_0_17_adc[1]) + Double_t(rightdE_adc[1]);
    Double_t y = Double_t(left_0_17_adc[2]) + Double_t(rightdE_adc[2]);
    if (x <= 0 || y <= 0)
      continue;
    if (!BeamFitUtils::InEllipseXY(beam, x, y, kEllipseNSigmaX,
                                   kEllipseNSigmaY))
      continue;
    if (pairs) {
      for (Int_t s = 1; s <= 16; s++) {
        if (Long64_t(pairs[s].gated_long.size()) >= kPairCap)
          continue;
        Bool_t l_is_long = (LongSide(s) == 'L');
        Int_t long_v =
            l_is_long ? Int_t(left_0_17_adc[s]) : Int_t(rightdE_adc[s]);
        Int_t short_v =
            l_is_long ? Int_t(rightdE_adc[s]) : Int_t(left_0_17_adc[s]);
        if (long_v > 0) {
          pairs[s].gated_long.push_back(Float_t(long_v));
          pairs[s].gated_short.push_back(Float_t(short_v > 0 ? short_v : 0));
        }
      }
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
//   * SHORT side anchor = C_long/|slope| of the charge-sharing ridge; see
//     RidgeShortAnchor. The short end is never observed collecting the full
//     deposit -- the beam is collimated onto the long end -- so the anchor is
//     reached by extrapolating the ridge to the short axis rather than by
//     finding a peak.
//   * gain = TARGET / anchor per side, with TARGET = 1.0 a.u. (the notebook
//     uses 1000 ADC; only the overall scale differs).
//
// Pass 2 — per-strip normalisation of the summed beam peak to 1.0 a.u.:
//   * pass 1 anchors the LONG side's own beam peak at 1.0, but at that peak
//     the long end only carries fraction (1-f) of the strip's charge, so the
//     summed total reads 1+f. f varies strip to strip (it grows downstream and
//     differs by parity), which is what shows up as a sawtooth in the summed
//     trace.
//   * the peak is measured on `gated_long`/`gated_short`: beam-gated events
//     where the LONG end fired, with a silent short end recorded as 0. That is
//     the population the decode sums. The `l`/`r` pairs must not be used here:
//     they require both ends to fire and so are biased towards large
//     short-side signals.
//   * both gains are scaled by 1/peak, which leaves the L/R ratio from pass 1
//     untouched, so a wrong short anchor misallocates charge within a strip
//     but no longer shifts the strip total.
//   * when IGNORE_SHORT_STRIPS is set the decode keeps only the long end, so
//     the long end alone is normalised to 1.0 and pass 1 already provides it.
//
// Pairs are collected UNGATED (all events with both ends firing) because the
// ridge needs the off-centre crossings that a beam gate removes. Strips whose
// ridge cannot be fitted fall back to the MEDIAN anchor of the strips that
// could; that fallback assumes a common electronics scale across strips, which
// the measured ridge slopes can be used to check rather than assume.
const Int_t kGmBins = 512;
// Ridge-slope short anchor. In the long-vs-short plane a fixed deposit divided
// between the two ends of a strip traces
//     ADC_long/C_long + ADC_short/C_short = 1,
// a line of slope -C_long/C_short. The slope is the same for every deposit
// energy, so the ridge direction alone gives the anchor ratio and C_short need
// never be observed directly -- which matters because the beam is collimated
// onto the long end and almost never deposits its full charge on the short one.
// The plane also contains pile-up bands at 2x, 3x the single-particle deposit;
// the 2-particle band falls into the single-particle window once
// short > ~0.55*C_short, so the fit is capped well below that.
const Double_t kRidgeShortMaxFrac =
    0.25;                           // cap on short, as a fraction of C_long
const Double_t kRidgeBandLo = 0.30; // single-particle band, x C_long
const Double_t kRidgeBandHi = 1.45;
const Int_t kRidgeSlices = 60;
const Long64_t kRidgeMinPerSlice = 200;
const Int_t kRidgeMinPts = 6;
// C_short/C_long is a preamp gain ratio, so it is order unity; outside this
// range the ridge fit has failed in a way the intercept test cannot see.
const Double_t kRidgeRatioLo = 0.30;
const Double_t kRidgeRatioHi = 3.00;
const Double_t kGmEsumLo = 0.8; // a.u. eSum peak search window (pass 2)
const Double_t kGmEsumHi = 2.5;

// Everything SaveRidgeFitPlots needs to redraw a fit: the slice medians it was
// fitted to, the window they were taken from, and the resulting line. Filled
// even when the fit is later rejected, so a bad fit can be looked at.
struct RidgeFit {
  std::vector<Double_t> x, y, ey; // slice centres, medians, median errors
  Double_t lo = 0.0, hi = 0.0;    // short-axis fit window
  Double_t slope = 0.0, intercept = 0.0;
  Double_t c_long = 0.0, c_short = 0.0;
  Bool_t fitted = kFALSE; // a line was fitted (it may still be rejected)
};

// Slope of the single-particle charge-sharing ridge, returning C_short via
// C_short = C_long/|slope|. Returns 0 when the ridge is not measurable.
Double_t RidgeShortAnchor(const std::vector<Float_t> &v_short,
                          const std::vector<Float_t> &v_long, Double_t c_long,
                          Double_t &slope_out, Double_t &intercept_out,
                          RidgeFit *dbg = nullptr) {
  slope_out = 0.0;
  intercept_out = 0.0;
  if (dbg)
    dbg->c_long = c_long;
  if (c_long <= 0 || v_short.size() != v_long.size() || v_short.size() < 500)
    return 0.0;
  const Double_t hi = kRidgeShortMaxFrac * c_long;
  const Double_t lo = 0.04 * c_long;
  std::vector<std::vector<Double_t>> slice(kRidgeSlices);
  for (Int_t j = 0; j < Int_t(v_short.size()); j++) {
    Double_t sh = Double_t(v_short[j]), lg = Double_t(v_long[j]);
    if (sh < lo || sh >= hi)
      continue;
    if (lg <= kRidgeBandLo * c_long || lg >= kRidgeBandHi * c_long)
      continue;
    Int_t b = Int_t((sh - lo) / (hi - lo) * kRidgeSlices);
    if (b >= 0 && b < kRidgeSlices)
      slice[b].push_back(lg);
  }
  std::vector<Double_t> x, y, ey;
  for (Int_t b = 0; b < kRidgeSlices; b++) {
    if (Long64_t(slice[b].size()) < kRidgeMinPerSlice)
      continue;
    std::sort(slice[b].begin(), slice[b].end());
    Double_t med = slice[b][slice[b].size() / 2];
    Double_t iqr =
        slice[b][slice[b].size() * 3 / 4] - slice[b][slice[b].size() / 4];
    x.push_back(lo + (b + 0.5) * (hi - lo) / kRidgeSlices);
    y.push_back(med);
    ey.push_back(1.253 * (iqr / 1.349) /
                 TMath::Sqrt(Double_t(slice[b].size())));
  }
  if (dbg) {
    dbg->x = x;
    dbg->y = y;
    dbg->ey = ey;
    dbg->lo = lo;
    dbg->hi = hi;
  }
  if (Int_t(x.size()) < kRidgeMinPts)
    return 0.0;
  TGraphErrors g(Int_t(x.size()), &x[0], &y[0], nullptr, &ey[0]);
  TF1 fit("f_ridge", "pol1", x.front(), x.back());
  if (g.Fit(&fit, "QN") != 0)
    return 0.0;
  Double_t slope = fit.GetParameter(1);
  Double_t inter = fit.GetParameter(0);
  slope_out = slope;
  intercept_out = inter;
  if (dbg) {
    dbg->slope = slope;
    dbg->intercept = inter;
    dbg->fitted = kTRUE;
    if (slope < 0)
      dbg->c_short = c_long / TMath::Abs(slope);
  }
  if (slope >= 0)
    return 0.0;
  // The line must pass through C_long on the long axis; a large departure means
  // the band selection did not isolate the single-particle ridge.
  if (inter < 0.80 * c_long || inter > 1.20 * c_long)
    return 0.0;
  return c_long / TMath::Abs(slope);
}

// One plot per strip under <plot_subdir>/ridge, named ridge_s<NN>: the
// beam-gated long-vs-short plane, the slice medians the fit was actually given
// (black), and the fitted line (violet). Drawn for every strip, including the
// ones whose fit was rejected, so a bad ridge can be seen rather than inferred
// from the slope in the log.
void SaveRidgeFitPlots(const StripPairSamples pairs[18], const RidgeFit dbg[18],
                       const TString &plot_subdir) {
  TString subdir = plot_subdir + "/ridge";
  for (Int_t s = 1; s <= 16; s++) {
    const RidgeFit &d = dbg[s];
    const StripPairSamples &p = pairs[s];
    if (d.c_long <= 0 || p.gated_short.empty())
      continue;
    Bool_t l_is_long = (LongSide(s) == 'L');
    // A little past the fit window, so the cap the fit stops at is visible.
    Double_t xhi = (d.hi > 0 ? d.hi : kRidgeShortMaxFrac * d.c_long) * 1.35;
    TH2F *h = new TH2F(Form("h_ridge_s%d", s),
                       Form(";Strip %d Short (%c) #DeltaE [ADC];Strip %d Long "
                            "(%c) #DeltaE [ADC]",
                            s, l_is_long ? 'R' : 'L', s, l_is_long ? 'L' : 'R'),
                       200, 0.0, xhi, 200, 0.0, 1.60 * d.c_long);
    for (Int_t j = 0; j < Int_t(p.gated_short.size()); j++)
      h->Fill(Double_t(p.gated_short[j]), Double_t(p.gated_long[j]));
    TCanvas *cv = PlottingUtils::GetConfiguredCanvas(kFALSE);
    PlottingUtils::Configure2DHistogram(h, cv);
    h->Draw("COLZ");
    TGraphErrors *g = nullptr;
    if (!d.x.empty()) {
      g = new TGraphErrors(Int_t(d.x.size()), &d.x[0], &d.y[0], nullptr,
                           &d.ey[0]);
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.8);
      g->SetMarkerColor(kBlack);
      g->SetLineColor(kBlack);
      g->Draw("P SAME");
    }
    // Extended to short = 0 on purpose: where the line lands there is the
    // anchor the calibration takes, and it should sit on C_long.
    TLine *lf = nullptr;
    if (d.fitted) {
      lf = new TLine(0.0, d.intercept, xhi, d.intercept + d.slope * xhi);
      lf->SetLineColor(kViolet + 2);
      lf->SetLineWidth(2);
      lf->Draw();
    }
    if (Constants::cfg.SAVE_PLOTS)
      PlottingUtils::SaveFigure(cv, Form("ridge_s%02d", s), subdir,
                                PlotSaveOptions::kLINEAR);
    delete cv;
    delete g;
    delete lf;
    delete h;
  }
}

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
                        const StripPairSamples pairs[18],
                        const TString &plot_subdir) {
  RidgeFit ridge_dbg[18];
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

  // ── Pass 1: per-side anchors ──
  // LONG side: its own beam peak, already fitted by ReduceToAnchors.
  // SHORT side: from the charge-sharing ridge slope, the same method on both
  //   parities. Strips whose ridge is not measurable fall back to the median
  //   of the strips that did measure one.

  Bool_t matched[18] = {kFALSE};
  Double_t short_anchor_adc[18] = {0};
  std::vector<Double_t> ratios_found[2];
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

    {
      const StripPairSamples &p = pairs[s];
      // Fit the ridge on beam-gated events. The gate cuts on the strip-1 and
      // strip-2 SUMS, so it removes pile-up and junk without touching where a
      // particle crossed on any other strip -- exactly the off-centre
      // crossings the ridge is made of. Fitting the ungated pairs instead
      // admits the 2- and 3-particle bands and a low-long background that
      // grows with short, which drags the slice medians down and steepens the
      // slope.
      const std::vector<Float_t> &v_short = p.gated_short;
      const std::vector<Float_t> &v_long = p.gated_long;
      Double_t slope = 0.0, inter = 0.0;
      peak_short = RidgeShortAnchor(v_short, v_long, c_long.fit_adc, slope,
                                    inter, &ridge_dbg[s]);
      if (peak_short > 0)
        std::cout << "  strip " << s << " ridge slope=" << Form("%.3f", slope)
                  << " intercept/C_long="
                  << Form("%.3f",
                          c_long.fit_adc > 0 ? inter / c_long.fit_adc : 0.0)
                  << "  short_anchor=" << Form("%.1f", peak_short) << " ADC"
                  << std::endl;
      else
        std::cerr << "  strip " << s
                  << ": ridge not measurable (slope=" << Form("%.3f", slope)
                  << ", intercept=" << Form("%.1f", inter) << ")" << std::endl;
    }

    if (peak_short <= 0)
      continue;
    // C_short/C_long is the ratio of the two preamp gains, so it is order
    // unity. A ridge too flat to measure still crosses the long axis near
    // C_long, so the intercept test above cannot catch it -- but it sends
    // C_long/|slope| to absurd values, which this does catch.
    Double_t ratio = peak_short / c_long.fit_adc;
    if (ratio < kRidgeRatioLo || ratio > kRidgeRatioHi) {
      std::cerr << "  strip " << s << ": ridge ratio " << Form("%.2f", ratio)
                << " outside [" << kRidgeRatioLo << ", " << kRidgeRatioHi
                << "]; rejecting anchor " << Form("%.1f", peak_short) << " ADC"
                << std::endl;
      continue;
    }
    short_anchor_adc[s] = peak_short;
    ratios_found[s % 2].push_back(ratio);
    chans[l_is_long ? idx_r[s] : idx_l[s]].ridge_ratio = ratio;
    std::cout << "  strip " << s << " short_anchor=" << Form("%.1f", peak_short)
              << " ADC (ratio " << Form("%.3f", ratio) << ")" << std::endl;
  }

  SaveRidgeFitPlots(pairs, ridge_dbg, plot_subdir);

  // Fall back on the median RATIO rather than the median anchor: the ratio is a
  // property of the two preamps, so it carries across strips, whereas an anchor
  // in ADC does not -- each strip has its own C_long. Taken per parity because
  // the L and R channels are on separate preamps and their ratios differ
  // systematically (odd ~1.2, even ~0.9 on 87Rb).
  Double_t median_ratio[2] = {0.0, 0.0};
  std::vector<Double_t> all_ratios;
  for (Int_t par = 0; par < 2; par++) {
    std::vector<Double_t> &v = ratios_found[par];
    all_ratios.insert(all_ratios.end(), v.begin(), v.end());
    if (v.size() < 2)
      continue;
    std::sort(v.begin(), v.end());
    Int_t m = Int_t(v.size());
    median_ratio[par] =
        (m % 2 == 1) ? v[m / 2] : 0.5 * (v[m / 2 - 1] + v[m / 2]);
  }
  Double_t global_ratio = 0.0;
  if (!all_ratios.empty()) {
    std::sort(all_ratios.begin(), all_ratios.end());
    Int_t m = Int_t(all_ratios.size());
    global_ratio = (m % 2 == 1)
                       ? all_ratios[m / 2]
                       : 0.5 * (all_ratios[m / 2 - 1] + all_ratios[m / 2]);
  }
  std::cout << "  ridge ratio medians: odd=" << Form("%.3f", median_ratio[1])
            << " even=" << Form("%.3f", median_ratio[0])
            << " all=" << Form("%.3f", global_ratio) << std::endl;

  for (Int_t s = 1; s <= 16; s++) {
    if (idx_l[s] < 0 || idx_r[s] < 0)
      continue;
    Bool_t l_is_long = (LongSide(s) == 'L');
    ChannelCal &c_long = chans[l_is_long ? idx_l[s] : idx_r[s]];
    ChannelCal &c_short = chans[l_is_long ? idx_r[s] : idx_l[s]];
    if (!IsCalibrated(c_long))
      continue;
    if (short_anchor_adc[s] <= 0) {
      Double_t r = median_ratio[s % 2] > 0 ? median_ratio[s % 2] : global_ratio;
      if (r <= 0) {
        std::cerr << "  strip " << s
                  << ": no ridge and no ratio fallback; keeping "
                     "independent gains"
                  << std::endl;
        continue;
      }
      short_anchor_adc[s] = r * c_long.fit_adc;
      std::cout << "  strip " << s << " short_anchor=ratio fallback "
                << Form("%.3f", r)
                << " x C_long = " << Form("%.1f", short_anchor_adc[s]) << " ADC"
                << (median_ratio[s % 2] > 0 ? "" : " (global, parity had none)")
                << std::endl;
    }
    c_long.gain = 1.0 / c_long.fit_adc;
    c_short.gain = 1.0 / short_anchor_adc[s];
    matched[s] = kTRUE;
    std::cout << "  strip " << s
              << " L/R match: long peak=" << Form("%.1f", c_long.fit_adc)
              << " ADC  short anchor=" << Form("%.1f", short_anchor_adc[s])
              << " ADC" << std::endl;
  }

  // ── Pass 2: check, do not correct ──
  // Both anchors are now measured (long from its beam peak, short from the
  // ridge slope), so gain_L*L + gain_R*R already peaks at 1.0 a.u. for a beam
  // event by construction. A departure means the ridge fit for that strip is
  // wrong, and is reported rather than absorbed into the gains -- rescaling
  // here would hide exactly the failure worth seeing.
  Double_t esum_peak[18] = {0};
  for (Int_t s = 1; s <= 16; s++) {
    if (!matched[s])
      continue;
    const StripPairSamples &p = pairs[s];
    Bool_t l_is_long = (LongSide(s) == 'L');
    Double_t g_long = Gain(chans[l_is_long ? idx_l[s] : idx_r[s]]);
    Double_t g_short = Gain(chans[l_is_long ? idx_r[s] : idx_l[s]]);
    std::vector<Float_t> esum;
    esum.reserve(p.gated_long.size());
    for (Int_t j = 0; j < Int_t(p.gated_long.size()); j++)
      esum.push_back(Float_t(g_long * Double_t(p.gated_long[j]) +
                             g_short * Double_t(p.gated_short[j])));
    esum_peak[s] = GmFindPeak(esum, kGmEsumLo, kGmEsumHi, 0.0);
  }
  for (Int_t s = 1; s <= 16; s++) {
    if (!matched[s])
      continue;
    if (esum_peak[s] <= 0) {
      std::cerr << "  strip " << s << ": no summed beam peak in (" << kGmEsumLo
                << ", " << kGmEsumHi << ") a.u." << std::endl;
      continue;
    }
    Double_t dev = esum_peak[s] - 1.0;
    std::cout << "  strip " << s
              << " summed beam peak=" << Form("%.4f", esum_peak[s]) << " a.u. ("
              << Form("%+.1f%%", 100.0 * dev) << ")"
              << (TMath::Abs(dev) > 0.05 ? "  <-- check ridge fit" : "")
              << std::endl;
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
                     const StripPairSamples pairs[18],
                     const TString &plot_subdir) {
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
    ComputeLRGainMatch(chans, pairs, plot_subdir);
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
  // Per-strip ridge ratio measured in THIS subfile, 0 where the ridge was not
  // measurable. AggregateRidgeRatiosForRun medians these across a run.
  Float_t ridge_ratio[18] = {0};
  Float_t long_anchor[18] = {0};
  Int_t n_actual = TMath::Min(Int_t(chans.size()), kMaxChannels);
  for (Int_t k = 0; k < n_actual; k++) {
    const ChannelCal &c = chans[k];
    ok[k] = IsCalibrated(c) || c.gain > 0;
    gain[k] = ok[k] ? Float_t(Gain(c)) : 0.0f;
    fit_adc[k] = Float_t(c.fit_adc);
    fit_sigma[k] = Float_t(c.fit_sigma_adc);
    fit_n[k] = c.n_samples;
    if (c.strip >= 1 && c.strip <= 16 && (c.side == 'L' || c.side == 'R')) {
      if (c.ridge_ratio > 0)
        ridge_ratio[c.strip] = Float_t(c.ridge_ratio);
      if (c.side == LongSide(c.strip))
        long_anchor[c.strip] = Float_t(c.fit_adc);
    }
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
  cal->Branch("RidgeRatio", ridge_ratio, "RidgeRatio[18]/F");
  cal->Branch("LongAnchor", long_anchor, "LongAnchor[18]/F");

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

  std::vector<ChannelCal> chans = chans_template;
  std::vector<std::vector<Float_t>> samples;
  StripPairSamples pairs[18];
  CollectAnchorSamplesOneSubfile(spec, chans, beam, samples, pairs);
  std::vector<TF1 *> peak_fits;
  {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    ReduceToAnchors(chans, samples, peak_fits, file_label, pairs, plot_subdir);
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

// Replace each subfile's short-side gain with one built from the run-level
// median ridge ratio.
//
// C_short/C_long is a ratio of preamp gains, so it is fixed per channel and
// does not vary subfile to subfile -- the measured scatter is ~1%, well below
// the ~7% spread between strips. But a single subfile often cannot fit the
// ridge on the low-occupancy short ends (on 87Rb the even strips fit in only
// ~45% of subfiles, and 37Cl's upstream strips almost never do), so per-subfile
// fitting leaves a large fraction of strips on a fallback.
//
// Aggregating fixes that without re-reading any event data: take the median
// ratio per strip over the subfiles that did measure it, then rebuild every
// subfile's short gain as 1/(ratio * C_long), using that subfile's own C_long
// so per-subfile gain drift is preserved. Strips that no subfile could fit keep
// whatever the per-subfile fallback gave them.
void CalibrateBeam::AggregateRidgeRatiosForRun(
    Int_t run, const std::vector<FileSpec> &specs) {
  std::vector<Double_t> per_strip[18];
  for (Int_t k = 0; k < Int_t(specs.size()); k++) {
    TString sub = FileSet::EventsName(specs[k]) + ".root";
    TFile *cf = IO::OpenForReading(sub);
    if (!cf || cf->IsZombie()) {
      if (cf)
        delete cf;
      continue;
    }
    TTree *t = static_cast<TTree *>(cf->Get("calibration"));
    if (!t || t->GetEntries() < 1 || !t->GetBranch("RidgeRatio")) {
      cf->Close();
      delete cf;
      continue;
    }
    Float_t rr[18] = {0};
    t->SetBranchAddress("RidgeRatio", rr);
    t->GetEntry(0);
    for (Int_t s = 1; s <= 16; s++)
      if (rr[s] > 0)
        per_strip[s].push_back(Double_t(rr[s]));
    cf->Close();
    delete cf;
  }

  Double_t med[18] = {0};
  std::cout << "Run " << run << " ridge-ratio aggregation:" << std::endl;
  for (Int_t s = 1; s <= 16; s++) {
    std::vector<Double_t> &v = per_strip[s];
    if (v.size() < 3) {
      std::cerr << "  strip " << s << ": only " << v.size()
                << " subfiles measured the ridge; leaving per-subfile gains"
                << std::endl;
      continue;
    }
    std::sort(v.begin(), v.end());
    Int_t m = Int_t(v.size());
    med[s] = (m % 2 == 1) ? v[m / 2] : 0.5 * (v[m / 2 - 1] + v[m / 2]);
    Double_t lo = v[m / 4], hi = v[(3 * m) / 4];
    std::cout << "  strip " << s << " ratio=" << Form("%.4f", med[s])
              << "  IQR " << Form("%.4f", lo) << "-" << Form("%.4f", hi)
              << "  from " << m << " subfiles" << std::endl;
  }

  Int_t n_rewritten = 0;
  for (Int_t k = 0; k < Int_t(specs.size()); k++) {
    TString sub = FileSet::EventsName(specs[k]) + ".root";
    TFile *cf = IO::OpenForWriting(sub, "UPDATE");
    if (!cf || cf->IsZombie()) {
      if (cf)
        delete cf;
      continue;
    }
    TTree *t = static_cast<TTree *>(cf->Get("calibration"));
    if (!t || t->GetEntries() < 1 || !t->GetBranch("LongAnchor")) {
      cf->Close();
      delete cf;
      continue;
    }
    Float_t gl[18] = {0}, gr[18] = {0}, la[18] = {0}, rr[18] = {0};
    t->SetBranchAddress("GainLeft", gl);
    t->SetBranchAddress("GainRight", gr);
    t->SetBranchAddress("LongAnchor", la);
    t->SetBranchAddress("RidgeRatio", rr);
    t->GetEntry(0);
    Bool_t changed = kFALSE;
    for (Int_t s = 1; s <= 16; s++) {
      if (med[s] <= 0 || la[s] <= 0)
        continue;
      Double_t anchor = med[s] * Double_t(la[s]);
      if (anchor <= 0)
        continue;
      if (LongSide(s) == 'L')
        gr[s] = Float_t(1.0 / anchor);
      else
        gl[s] = Float_t(1.0 / anchor);
      rr[s] = Float_t(med[s]);
      changed = kTRUE;
    }
    if (changed) {
      TTree *nt = t->CloneTree(0);
      nt->Fill();
      cf->cd();
      nt->Write("calibration", TObject::kOverwrite);
      n_rewritten++;
    }
    cf->Close();
    delete cf;
  }
  std::cout << "  rewrote short gains in " << n_rewritten << " subfiles"
            << std::endl;
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
