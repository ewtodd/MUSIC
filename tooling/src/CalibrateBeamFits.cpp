#include "CalibrateBeamInternal.hpp"

// nth_element permutes its input, so both helpers work on a copy and leave
// the caller's sample vector untouched.
Double_t Median(const std::vector<Float_t> &v_in) {
  if (v_in.empty())
    return 0.0;
  std::vector<Float_t> v(v_in);
  Int_t n = Int_t(v.size());
  std::nth_element(v.begin(), v.begin() + n / 2, v.end());
  Double_t med = Double_t(v[n / 2]);
  if (n % 2 == 0) {
    std::nth_element(v.begin(), v.begin() + n / 2 - 1, v.end());
    med = 0.5 * (med + Double_t(v[n / 2 - 1]));
  }
  return med;
}

Double_t InterquartileRange(const std::vector<Float_t> &v_in) {
  std::vector<Float_t> v(v_in);
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

// Robust peak/width: mode centre of the 5th–95th-percentile core, IQR/1.349
// width — not pulled by the straggling beam-dE tail.
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

// Per-channel beam-peak histogram shared by fitters and diagnostic plots:
// range = mode ± 4σ (~6 bins/sigma), never finer than 1 ADC.
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

// Mode-seeded ±2σ Gaussian fit; a narrow refit follows so bimodal spectra
// (e.g. Strip0 bump) lock onto the dominant peak, not its average.
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

  // Refit inside μ ± min(1.5σ_fit, 12% of μ): clean Gaussian reproduces the
  // first pass; contaminated spectrum excludes the secondary component.
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
  // Peak-like sanity gate: sigma > 50% of centroid means no real peak
  // (e.g. 37Cl Run84 Strip17 smeared spectrum).
  if (sigma_adc > 0.5 * peak_adc) {
    delete f;
    return kFALSE;
  }
  fit_out = f;
  return kTRUE;
}
