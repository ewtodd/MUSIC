#include "BeamFit2D.hpp"

void BeamFitUtils::DiagonalizeCov(const BeamFit2D &b, Double_t &major,
                                  Double_t &minor, Double_t &theta_deg) {
  Double_t vx = b.sigma_x * b.sigma_x;
  Double_t vy = b.sigma_y * b.sigma_y;
  Double_t cxy = b.rho * b.sigma_x * b.sigma_y;
  Double_t disc = std::sqrt(0.25 * (vx - vy) * (vx - vy) + cxy * cxy);
  Double_t lam1 = 0.5 * (vx + vy) + disc;
  Double_t lam2 = 0.5 * (vx + vy) - disc;
  theta_deg = 0.5 * std::atan2(2.0 * cxy, vx - vy) * 180.0 / TMath::Pi();
  major = std::sqrt(std::max(lam1, 0.0));
  minor = std::sqrt(std::max(lam2, 0.0));
}

Bool_t BeamFitUtils::InEllipse(const BeamFit2D &b, Double_t x, Double_t y,
                               Double_t n_sigma) {
  Double_t dx = (x - b.mu_x) / b.sigma_x;
  Double_t dy = (y - b.mu_y) / b.sigma_y;
  Double_t d2 =
      (dx * dx - 2.0 * b.rho * dx * dy + dy * dy) / (1.0 - b.rho * b.rho);
  return d2 < n_sigma * n_sigma;
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

BeamFit2D BeamFitUtils::FitBigausInWindow(TH2F *h, Double_t x_lo, Double_t x_hi,
                                          Double_t y_lo, Double_t y_hi,
                                          const Moments2D &seed,
                                          Double_t seed_amp,
                                          const TString &fname,
                                          const TString &tag) {
  BeamFit2D out;
  TF2 *f = new TF2(fname, "bigaus", x_lo, x_hi, y_lo, y_hi);
  f->SetParameters(seed_amp, seed.mu_x, seed.sigma_x, seed.mu_y, seed.sigma_y,
                   seed.rho);
  f->SetParLimits(5, -0.99, 0.99);
  TFitResultPtr r = h->Fit(f, "QSRN");
  if (r.Get() && r->IsValid()) {
    out.amp = f->GetParameter(0);
    out.mu_x = f->GetParameter(1);
    out.sigma_x = f->GetParameter(2);
    out.mu_y = f->GetParameter(3);
    out.sigma_y = f->GetParameter(4);
    out.rho = f->GetParameter(5);
    out.ok = kTRUE;
  } else {
    std::cerr << "  " << tag << ": bigaus fit invalid; using seed" << std::endl;
    out.amp = seed_amp;
    out.mu_x = seed.mu_x;
    out.mu_y = seed.mu_y;
    out.sigma_x = seed.sigma_x;
    out.sigma_y = seed.sigma_y;
    out.rho = seed.rho;
    out.ok = kTRUE;
  }
  delete f;
  return out;
}
