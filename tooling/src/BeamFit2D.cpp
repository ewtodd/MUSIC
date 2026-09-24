#include "BeamFit2D.hpp"
#include <TEllipse.h>
#include <TF1.h>
#include <TF2.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMath.h>

Bool_t BeamFitUtils::InEllipseXY(const BeamFit2D &b, Double_t x, Double_t y,
                                 Double_t n) {
  Double_t dx = x - b.mu_x;
  Double_t dy = y - b.mu_y;
  Double_t sx = b.sigma_x, sy = b.sigma_y, rho = b.rho;
  if (sx <= 0 || sy <= 0)
    return kFALSE;
  // Correlated 2D Gaussian (Mahalanobis) chi2: inside n sigma is chi2 < n^2.
  Double_t dx_s = dx / sx, dy_s = dy / sy;
  Double_t r2 = rho * rho;
  Double_t chi2 =
      (dx_s * dx_s + dy_s * dy_s - 2.0 * rho * dx_s * dy_s) / (1.0 - r2);
  return chi2 < n * n;
}

void BeamFitUtils::DrawEllipse(const BeamFit2D &b, Double_t n, Color_t color) {
  Double_t sxx = b.sigma_x * b.sigma_x;
  Double_t syy = b.sigma_y * b.sigma_y;
  Double_t sxy = b.rho * b.sigma_x * b.sigma_y;
  Double_t sum = sxx + syy;
  Double_t diff = sxx - syy;
  Double_t det = TMath::Sqrt(diff * diff + 4.0 * sxy * sxy);
  Double_t lambda1 = 0.5 * (sum + det);
  Double_t lambda2 = 0.5 * (sum - det);
  Double_t theta = 0.5 * TMath::ATan2(2.0 * sxy, diff) * 180.0 / TMath::Pi();
  TEllipse *e = new TEllipse(b.mu_x, b.mu_y, n * TMath::Sqrt(lambda1),
                             n * TMath::Sqrt(lambda2), 0, 360, theta);
  e->SetFillStyle(0);
  e->SetLineColor(color);
  e->SetLineWidth(2);
  e->Draw();
}

Moments2D BeamFitUtils::ComputeMoments(TH2F *h, Int_t lo_bx, Int_t hi_bx,
                                       Int_t lo_by, Int_t hi_by,
                                       Double_t thresh, Double_t bw_x,
                                       Double_t bw_y) {
  Moments2D m;
  Double_t W = 0, Mx = 0, My = 0, Cxx = 0, Cyy = 0, Cxy = 0;
  for (Int_t ix = lo_bx; ix <= hi_bx; ix++) {
    Double_t x = h->GetXaxis()->GetBinCenter(ix);
    for (Int_t iy = lo_by; iy <= hi_by; iy++) {
      Double_t w = h->GetBinContent(ix, iy);
      if (w < thresh)
        continue;
      Double_t y = h->GetYaxis()->GetBinCenter(iy);
      W += w;
      Mx += w * x;
      My += w * y;
    }
  }
  if (W <= 0)
    return m;
  Mx /= W;
  My /= W;
  for (Int_t ix = lo_bx; ix <= hi_bx; ix++) {
    Double_t x = h->GetXaxis()->GetBinCenter(ix);
    for (Int_t iy = lo_by; iy <= hi_by; iy++) {
      Double_t w = h->GetBinContent(ix, iy);
      if (w < thresh)
        continue;
      Double_t y = h->GetYaxis()->GetBinCenter(iy);
      Double_t dx = x - Mx, dy = y - My;
      Cxx += w * dx * dx;
      Cyy += w * dy * dy;
      Cxy += w * dx * dy;
    }
  }
  Cxx /= W;
  Cyy /= W;
  Cxy /= W;
  m.mu_x = Mx;
  m.mu_y = My;
  m.sigma_x = std::max(std::sqrt(std::max(Cxx, 0.0)), 2.0 * bw_x);
  m.sigma_y = std::max(std::sqrt(std::max(Cyy, 0.0)), 2.0 * bw_y);
  Double_t r = Cxy / (m.sigma_x * m.sigma_y);
  if (r > 0.95)
    r = 0.95;
  if (r < -0.95)
    r = -0.95;
  m.rho = r;
  m.weight = W;
  return m;
}

Moments2D BeamFitUtils::ClippedMoments(
    const std::vector<std::pair<Float_t, Float_t>> &pts, const Moments2D &seed,
    Double_t clip) {
  const Int_t kMaxPasses = 50;
  const Double_t kTolerance = 1.0e-4;
  // A 2-D Gaussian's r^2 is exponential with mean 2: inside radius c each
  // axis keeps 0.949 of its variance at c = 3.
  const Double_t e = std::exp(-0.5 * clip * clip);
  const Double_t kept_var = 1.0 - 0.5 * clip * clip * e / (1.0 - e);
  Moments2D m = seed;
  for (Int_t pass = 0; pass < kMaxPasses; pass++) {
    if (!(m.sigma_x > 0.0 && m.sigma_y > 0.0))
      return seed;
    Double_t W = 0, Sx = 0, Sy = 0, Sxx = 0, Syy = 0, Sxy = 0;
    const Double_t one_minus_r2 = 1.0 - m.rho * m.rho;
    for (size_t k = 0; k < pts.size(); k++) {
      const Double_t x = pts[k].first, y = pts[k].second;
      const Double_t u = (x - m.mu_x) / m.sigma_x;
      const Double_t v = (y - m.mu_y) / m.sigma_y;
      if ((u * u + v * v - 2.0 * m.rho * u * v) / one_minus_r2 >= clip * clip)
        continue;
      W += 1.0;
      Sx += x;
      Sy += y;
      Sxx += x * x;
      Syy += y * y;
      Sxy += x * y;
    }
    if (W < 10.0)
      return seed;
    Moments2D next;
    next.mu_x = Sx / W;
    next.mu_y = Sy / W;
    const Double_t cxx = std::max(Sxx / W - next.mu_x * next.mu_x, 0.0);
    const Double_t cyy = std::max(Syy / W - next.mu_y * next.mu_y, 0.0);
    const Double_t cxy = Sxy / W - next.mu_x * next.mu_y;
    next.sigma_x = std::sqrt(cxx / kept_var);
    next.sigma_y = std::sqrt(cyy / kept_var);
    if (!(next.sigma_x > 0.0 && next.sigma_y > 0.0))
      return seed;
    next.rho = std::max(-0.95, std::min(0.95, cxy / std::sqrt(cxx * cyy)));
    next.weight = W;
    const Bool_t done =
        std::fabs(next.sigma_x / m.sigma_x - 1.0) < kTolerance &&
        std::fabs(next.sigma_y / m.sigma_y - 1.0) < kTolerance &&
        std::fabs(next.rho - m.rho) < kTolerance;
    m = next;
    if (done)
      break;
  }
  return m;
}

Moments2D BeamFitUtils::SeedSpotMoments(TH2F *h, Int_t seed_half_bins) {
  // The seed threshold: a fraction of the peak bin.
  const Double_t kSeedFrac = 0.3;
  const Double_t bw_x = h->GetXaxis()->GetBinWidth(1);
  const Double_t bw_y = h->GetYaxis()->GetBinWidth(1);
  Int_t bx, by, bz;
  h->GetMaximumBin(bx, by, bz);
  const Double_t peak_val = h->GetBinContent(bx, by);
  Int_t lo_bx = std::max(1, bx - seed_half_bins);
  Int_t hi_bx = std::min(h->GetNbinsX(), bx + seed_half_bins);
  Int_t lo_by = std::max(1, by - seed_half_bins);
  Int_t hi_by = std::min(h->GetNbinsY(), by + seed_half_bins);
  return ComputeMoments(h, lo_bx, hi_bx, lo_by, hi_by, kSeedFrac * peak_val,
                        bw_x, bw_y);
}

BeamFitUtils::SpotFit
BeamFitUtils::FitSpotFromPoints(TH2F *h,
                                std::vector<std::pair<Float_t, Float_t>> &pts,
                                const Moments2D &seed) {
  // The core moments only seed the clip; the clip seeds the fit and is the
  // fallback.
  Moments2D m = ClippedMoments(pts, seed);

  // The width the fit reports is a fitted sigma on the inner kFitWindow,
  // finely binned so pileup and the reaction shoulder stay out of it.
  const Double_t kSpotWindow = 4.0;
  const Double_t kFitWindow = 2.0;
  const Int_t kSpotBins = 160;
  TH2F *hf = new TH2F(Form("%s_fit", h->GetName()), "", kSpotBins,
                      m.mu_x - kSpotWindow * m.sigma_x,
                      m.mu_x + kSpotWindow * m.sigma_x, kSpotBins,
                      m.mu_y - kSpotWindow * m.sigma_y,
                      m.mu_y + kSpotWindow * m.sigma_y);
  hf->SetDirectory(nullptr);
  for (size_t k = 0; k < pts.size(); k++)
    hf->Fill(pts[k].first, pts[k].second);
  std::vector<std::pair<Float_t, Float_t>>().swap(pts);
  SpotFit out;
  BeamFit2D fit = FitSpot(hf, m, kFitWindow, &out.chi2_ndf);
  delete hf;
  if (fit.ok) {
    out.fit = fit;
    out.fit_used = kTRUE;
    return out;
  }
  // A fit that ran away is dropped for the clipped moments.
  Int_t bx, by, bz;
  h->GetMaximumBin(bx, by, bz);
  out.fit.amp = h->GetBinContent(bx, by);
  out.fit.mu_x = m.mu_x;
  out.fit.mu_y = m.mu_y;
  out.fit.sigma_x = m.sigma_x;
  out.fit.sigma_y = m.sigma_y;
  out.fit.rho = m.rho;
  out.fit.ok = kTRUE;
  return out;
}

BeamFit2D BeamFitUtils::FitSpot(TH2F *h, const Moments2D &seed,
                                Double_t window_nsigma, Double_t *chi2_ndf) {
  // Off the limits by this much still counts as a parameter Minuit chose.
  const Double_t kOffLimit = 1.0e-3;
  const Double_t kSigmaLo = 0.2, kSigmaHi = 5.0;
  BeamFit2D out;
  if (chi2_ndf)
    *chi2_ndf = -1.0;
  if (!h || !(seed.sigma_x > 0.0 && seed.sigma_y > 0.0))
    return out;

  const Double_t x_lo = seed.mu_x - window_nsigma * seed.sigma_x;
  const Double_t x_hi = seed.mu_x + window_nsigma * seed.sigma_x;
  const Double_t y_lo = seed.mu_y - window_nsigma * seed.sigma_y;
  const Double_t y_hi = seed.mu_y + window_nsigma * seed.sigma_y;
  const Double_t peak = TMath::Max(h->GetMaximum(), 1.0);

  // bigaus in its parameter order, plus [6], a flat pedestal.
  TF2 f(Form("beamspot_fit_%s", h->GetName()),
        "[0]*exp(-0.5/(1.0-[5]*[5])*((x-[1])*(x-[1])/([2]*[2])"
        "+(y-[3])*(y-[3])/([4]*[4])"
        "-2.0*[5]*(x-[1])*(y-[3])/([2]*[4])))+[6]",
        x_lo, x_hi, y_lo, y_hi);
  f.SetParameters(peak, seed.mu_x, seed.sigma_x, seed.mu_y, seed.sigma_y,
                  TMath::Max(-0.9, TMath::Min(0.9, seed.rho)), 0.0);
  // Unbounded, Minuit runs the amplitude to a bound and walks the mean out.
  f.SetParLimits(0, 0.0, 1.0e3 * peak);
  f.SetParLimits(1, x_lo, x_hi);
  f.SetParLimits(2, kSigmaLo * seed.sigma_x, kSigmaHi * seed.sigma_x);
  f.SetParLimits(3, y_lo, y_hi);
  f.SetParLimits(4, kSigmaLo * seed.sigma_y, kSigmaHi * seed.sigma_y);
  f.SetParLimits(5, -0.95, 0.95);
  f.SetParLimits(6, 0.0, peak);
  // "R" fits the window only; the rest of h is the goodness check below.
  TFitResultPtr r = h->Fit(&f, "QNRSL");
  if (!r.Get() || !r->IsValid())
    return out;

  // Judged over the whole histogram: a shoulder the window kept out of the
  // width still shows here.
  if (chi2_ndf) {
    Double_t chi2 = 0.0;
    Int_t nbin = 0;
    for (Int_t ix = 1; ix <= h->GetNbinsX(); ix++)
      for (Int_t iy = 1; iy <= h->GetNbinsY(); iy++) {
        const Double_t e = f.Eval(h->GetXaxis()->GetBinCenter(ix),
                                  h->GetYaxis()->GetBinCenter(iy));
        if (!(e > 1.0))
          continue;
        const Double_t o = h->GetBinContent(ix, iy);
        chi2 += (o - e) * (o - e) / e;
        nbin++;
      }
    if (nbin > 7)
      *chi2_ndf = chi2 / Double_t(nbin - 7);
  }

  const Double_t sx = TMath::Abs(f.GetParameter(2));
  const Double_t sy = TMath::Abs(f.GetParameter(4));
  const Double_t mx = f.GetParameter(1), my = f.GetParameter(3);
  const Double_t rho = f.GetParameter(5);
  // Minuit reports success on fits that ran away: every parameter off its
  // limit, and finite.
  const Bool_t sane =
      TMath::Finite(sx) && TMath::Finite(sy) && TMath::Finite(rho) &&
      sx > (kSigmaLo + kOffLimit) * seed.sigma_x &&
      sx < (kSigmaHi - kOffLimit) * seed.sigma_x &&
      sy > (kSigmaLo + kOffLimit) * seed.sigma_y &&
      sy < (kSigmaHi - kOffLimit) * seed.sigma_y && TMath::Abs(rho) < 0.94 &&
      mx > x_lo && mx < x_hi && my > y_lo && my < y_hi;
  if (!sane)
    return out;

  out.amp = f.GetParameter(0);
  out.mu_x = mx;
  out.sigma_x = sx;
  out.mu_y = my;
  out.sigma_y = sy;
  out.rho = rho;
  out.ok = kTRUE;
  return out;
}

Bool_t BeamFitUtils::FitPeak(TH1F *h, Double_t seed_mu, Double_t seed_sigma,
                             Double_t window_nsigma, Bool_t fix_mu,
                             Double_t &mu, Double_t &sigma,
                             Double_t *chi2_ndf) {
  const Double_t kOffLimit = 1.0e-3;
  const Double_t kSigmaLo = 0.2, kSigmaHi = 5.0;
  if (chi2_ndf)
    *chi2_ndf = -1.0;
  if (!h || !(seed_sigma > 0.0))
    return kFALSE;
  const Double_t lo = seed_mu - window_nsigma * seed_sigma;
  const Double_t hi = seed_mu + window_nsigma * seed_sigma;
  const Double_t peak = TMath::Max(h->GetMaximum(), 1.0);
  TF1 f(Form("beampeak_fit_%s", h->GetName()), "gaus(0)+[3]", lo, hi);
  f.SetParameters(peak, seed_mu, seed_sigma, 0.0);
  f.SetParLimits(0, 0.0, 1.0e3 * peak);
  f.SetParLimits(1, lo, hi);
  f.SetParLimits(2, kSigmaLo * seed_sigma, kSigmaHi * seed_sigma);
  f.SetParLimits(3, 0.0, peak);
  if (fix_mu)
    f.FixParameter(1, seed_mu);
  TFitResultPtr r = h->Fit(&f, "QNRSL");
  if (!r.Get() || !r->IsValid())
    return kFALSE;
  if (chi2_ndf) {
    Double_t chi2 = 0.0;
    Int_t nbin = 0;
    for (Int_t ix = 1; ix <= h->GetNbinsX(); ix++) {
      const Double_t e = f.Eval(h->GetXaxis()->GetBinCenter(ix));
      if (!(e > 1.0))
        continue;
      const Double_t o = h->GetBinContent(ix);
      chi2 += (o - e) * (o - e) / e;
      nbin++;
    }
    if (nbin > 4)
      *chi2_ndf = chi2 / Double_t(nbin - 4);
  }
  const Double_t fit_mu = f.GetParameter(1);
  const Double_t fit_sigma = TMath::Abs(f.GetParameter(2));
  // Minuit reports success on fits that ran away: every parameter off its
  // limit, and finite.
  const Bool_t sane = TMath::Finite(fit_mu) && TMath::Finite(fit_sigma) &&
                      fit_sigma > (kSigmaLo + kOffLimit) * seed_sigma &&
                      fit_sigma < (kSigmaHi - kOffLimit) * seed_sigma &&
                      fit_mu > lo && fit_mu < hi;
  if (!sane)
    return kFALSE;
  mu = fit_mu;
  sigma = fit_sigma;
  return kTRUE;
}
