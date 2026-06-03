#include "BeamFit2D.hpp"
#include "Normalization.hpp"

namespace {
// Strip0-vs-Strip1-long-side beam gate geometry, used by delta-e-scatter to fit
// the calibration beam ellipse.
const Double_t kGateMin = 0.0;
const Double_t kGateMax = 12.0;
const Int_t kGateBins = 240;
const Double_t kGateNSigma = 2.0;
const Int_t kSeedHalfBins = 40;
const Double_t kSeedFrac = 0.30;
} // namespace

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

Bool_t BeamFitUtils::InEllipseXY(const BeamFit2D &b, Double_t x, Double_t y,
                                 Double_t nx, Double_t ny) {
  Double_t dx = (x - b.mu_x) / (nx * b.sigma_x);
  Double_t dy = (y - b.mu_y) / (ny * b.sigma_y);
  return dx * dx + dy * dy < 1.0;
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

Double_t BeamFitUtils::GateNSigma() { return kGateNSigma; }

// Pass 1: build the gate H2 (S0 vs L1) for the chain.
TH2F *BeamFitUtils::BuildGateHist(TChain *chain, const TString &name) {
  EnergyView ev;
  ev.Attach(chain);

  TH2F *h = new TH2F(name, ";#DeltaE S0 [MeV];#DeltaE L1 [MeV]", kGateBins,
                     kGateMin, kGateMax, kGateBins, kGateMin, kGateMax);
  h->SetDirectory(nullptr);

  Long64_t n = chain->GetEntries();
  for (Long64_t j = 0; j < n; j++) {
    chain->GetEntry(j);
    ev.Decode();
    Double_t x = ev.left[1];
    Double_t y = ev.right[2];
    if (x > 0 && y > 0)
      h->Fill(x, y);
  }
  return h;
}

// Moments-seeded bigaus fit on the single largest peak.
BeamFit2D BeamFitUtils::FitBigPeak(TH2F *h, const TString &tag) {
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

  Moments2D m = ComputeMoments(h, lo_bx, hi_bx, lo_by, hi_by,
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

  out = FitBigausInWindow(h, x_lo, x_hi, y_lo, y_hi, m, peak_val,
                          Form("f_gate_%s", tag.Data()),
                          TString("FitBigPeak ") + tag);

  std::cout << "  gate fit " << tag << ": mu=(" << out.mu_x << "," << out.mu_y
            << ") sigma=(" << out.sigma_x << "," << out.sigma_y
            << ") rho=" << out.rho << std::endl;
  return out;
}

void BeamFitUtils::DrawGateEllipse(const BeamFit2D &gate, Double_t n_sigma) {
  if (!gate.ok)
    return;
  Double_t major, minor, theta_deg;
  DiagonalizeCov(gate, major, minor, theta_deg);
  TEllipse *e = new TEllipse(gate.mu_x, gate.mu_y, n_sigma * major,
                             n_sigma * minor, 0.0, 360.0, theta_deg);
  e->SetFillStyle(0);
  e->SetLineColor(kRed + 1);
  e->SetLineWidth(2);
  e->Draw();
}
