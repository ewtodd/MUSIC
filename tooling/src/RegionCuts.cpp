#include "RegionCuts.hpp"
#include "Paths.hpp"
#include "PlottingUtils.hpp"
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TParameter.h>
#include <TSystem.h>
#include <cmath>
#include <iostream>
#include <vector>

// ---------------------------------------------------------------------------
// Store
// ---------------------------------------------------------------------------

namespace RegionCutStore {

static TString Key(const char *name, Int_t reac) {
  return Form("%s_reac%d", name, reac);
}

TString Dir() { return Paths::ResultsDir() + "/root_files/region_cuts"; }

TString Path(const char *name, Int_t reac) {
  return Dir() + "/" + Key(name, reac) + ".root";
}

static TString LegacyPath() {
  return Paths::ResultsDir() + "/root_files/RegionCuts.root";
}

static void WriteOne(const char *name, Int_t reac, TCutG *cut,
                     Double_t n_assigned) {
  if (!cut)
    return;
  TString path = Path(name, reac);
  TFile f(path, "RECREATE");
  if (f.IsZombie()) {
    std::cerr << "  [region] cannot open " << path << " to save cut"
              << std::endl;
    return;
  }
  f.cd();
  cut->Write(name);
  if (n_assigned >= 0.0)
    TParameter<Double_t>("n_assigned", n_assigned).Write();
  f.Close();
  std::cout << "  [region] saved " << name << " reac " << reac << " -> " << path
            << std::endl;
}

void Save(Int_t reac, TCutG *cut_an, TCutG *cut_aa, Double_t n_an_assigned) {
  gSystem->mkdir(Dir(), kTRUE);
  WriteOne("region_an", reac, cut_an, n_an_assigned);
  WriteOne("region_aa", reac, cut_aa, -1.0);
}

static const char *kFitKeys[14] = {
    "fit_beam_amp", "fit_beam_mx",  "fit_beam_sx", "fit_beam_my", "fit_beam_sy",
    "fit_beam_rho", "fit_reac_amp", "fit_reac_mx", "fit_reac_sx", "fit_reac_my",
    "fit_reac_sy",  "fit_reac_rho", "fit_n_beam",  "fit_n_reac"};

static void FitValues(const RegionFit &fit, Double_t *v) {
  const Gauss2D *g[2] = {&fit.beam, &fit.reac};
  for (Int_t k = 0; k < 2; k++) {
    v[6 * k + 0] = g[k]->amp;
    v[6 * k + 1] = g[k]->mx;
    v[6 * k + 2] = g[k]->sx;
    v[6 * k + 3] = g[k]->my;
    v[6 * k + 4] = g[k]->sy;
    v[6 * k + 5] = g[k]->rho;
  }
  v[12] = fit.n_beam;
  v[13] = fit.n_reac;
}

void SaveFit(Int_t reac, const RegionFit &fit) {
  TString path = Path("region_an", reac);
  TFile f(path, "UPDATE");
  if (f.IsZombie()) {
    std::cerr << "  [region] cannot open " << path << " to save the fit"
              << std::endl;
    return;
  }
  Double_t v[14];
  FitValues(fit, v);
  f.cd();
  for (Int_t k = 0; k < 14; k++)
    TParameter<Double_t>(kFitKeys[k], v[k])
        .Write(kFitKeys[k], TObject::kOverwrite);
  f.Close();
}

Bool_t LoadFit(Int_t reac, RegionFit &fit) {
  TString path = Path("region_an", reac);
  if (gSystem->AccessPathName(path))
    return kFALSE;
  TFile f(path, "READ");
  if (f.IsZombie())
    return kFALSE;
  Double_t v[14];
  for (Int_t k = 0; k < 14; k++) {
    TParameter<Double_t> *p =
        dynamic_cast<TParameter<Double_t> *>(f.Get(kFitKeys[k]));
    if (!p)
      return kFALSE;
    v[k] = p->GetVal();
  }
  Gauss2D *g[2] = {&fit.beam, &fit.reac};
  for (Int_t k = 0; k < 2; k++) {
    g[k]->amp = v[6 * k + 0];
    g[k]->mx = v[6 * k + 1];
    g[k]->sx = v[6 * k + 2];
    g[k]->my = v[6 * k + 3];
    g[k]->sy = v[6 * k + 4];
    g[k]->rho = v[6 * k + 5];
  }
  fit.n_beam = v[12];
  fit.n_reac = v[13];
  fit.ok = kTRUE;
  return kTRUE;
}

Double_t LoadAssigned(const char *name, Int_t reac) {
  TString path = Path(name, reac);
  if (gSystem->AccessPathName(path))
    return -1.0;
  TFile f(path, "READ");
  if (f.IsZombie())
    return -1.0;
  TParameter<Double_t> *p =
      dynamic_cast<TParameter<Double_t> *>(f.Get("n_assigned"));
  const Double_t n = p ? p->GetVal() : -1.0;
  f.Close();
  return n;
}

static TCutG *ReadFrom(const TString &path, const char *key, const char *name) {
  if (gSystem->AccessPathName(path))
    return nullptr;
  TFile f(path, "READ");
  if (f.IsZombie())
    return nullptr;
  TCutG *stored = dynamic_cast<TCutG *>(f.Get(key));
  if (!stored) {
    f.Close();
    return nullptr;
  }
  // The file owns the object; hand back a copy that outlives the close.
  TCutG *cut = static_cast<TCutG *>(stored->Clone(name));
  f.Close();
  return cut;
}

TCutG *Load(const char *name, Int_t reac) {
  TCutG *cut = ReadFrom(Path(name, reac), name, name);
  if (cut)
    return cut;
  return ReadFrom(LegacyPath(), Key(name, reac), name);
}

// Hand-drawn only. A per-cut file written by the interactive draw carries no
// n_assigned; one written by compute-regions does and is skipped. Then the
// pre-split RegionCuts.root, which only ever held drawn cuts.
TCutG *LoadDrawn(const char *name, Int_t reac) {
  TString path = Path(name, reac);
  if (!gSystem->AccessPathName(path)) {
    TFile f(path, "READ");
    const Bool_t fitted = !f.IsZombie() && f.Get("n_assigned") != nullptr;
    f.Close();
    if (!fitted)
      if (TCutG *cut = ReadFrom(path, name, name))
        return cut;
  }
  return ReadFrom(LegacyPath(), Key(name, reac), name);
}

} // namespace RegionCutStore

// ---------------------------------------------------------------------------
// Finder
// ---------------------------------------------------------------------------

namespace {

// ROOT's "bigaus": peak height `amp` in counts per bin, no normalisation.
Double_t Bigaus(const Gauss2D &g, Double_t x, Double_t y) {
  const Double_t dx = (x - g.mx) / g.sx, dy = (y - g.my) / g.sy;
  const Double_t r2 = 1.0 - g.rho * g.rho;
  return g.amp *
         std::exp(-0.5 / r2 * (dx * dx + dy * dy - 2.0 * g.rho * dx * dy));
}

// Squared Mahalanobis distance of (x, y) from a component.
Double_t Mahal2(const Gauss2D &g, Double_t x, Double_t y) {
  const Double_t dx = (x - g.mx) / g.sx, dy = (y - g.my) / g.sy;
  return (dx * dx + dy * dy - 2.0 * g.rho * dx * dy) / (1.0 - g.rho * g.rho);
}

// log of the bivariate normalisation, so densities of two components with
// different widths compare fairly.
Double_t LogNorm(const Gauss2D &g) {
  return std::log(2.0 * TMath::Pi() * g.sx * g.sy *
                  std::sqrt(1.0 - g.rho * g.rho));
}

// A component from weighted moments (sum w, sum wx, sum wy, sum wxx, sum wyy,
// sum wxy). Widths are floored at a bin, the correlation clamped, and the
// amplitude is the peak height in counts per bin of a Gaussian holding sw
// events, which is what ModelHist needs to draw it on the data's scale.
Gauss2D MomentsToGauss(Double_t sw, Double_t sx, Double_t sy, Double_t sxx,
                       Double_t syy, Double_t sxy, Double_t bwx, Double_t bwy) {
  Gauss2D g;
  g.mx = sx / sw;
  g.my = sy / sw;
  const Double_t vx = sxx / sw - g.mx * g.mx, vy = syy / sw - g.my * g.my;
  const Double_t cv = sxy / sw - g.mx * g.my;
  g.sx = std::sqrt(TMath::Max(vx, bwx * bwx));
  g.sy = std::sqrt(TMath::Max(vy, bwy * bwy));
  g.rho = TMath::Max(-0.95, TMath::Min(0.95, cv / (g.sx * g.sy)));
  g.amp = sw * bwx * bwy /
          (2.0 * TMath::Pi() * g.sx * g.sy * std::sqrt(1.0 - g.rho * g.rho));
  return g;
}

// FWHM-based width of a projection around its maximum (over its set range).
Double_t WidthAtMax(TH1D *h, Double_t &at) {
  Int_t b = h->GetMaximumBin();
  Double_t c = h->GetBinContent(b);
  at = h->GetBinCenter(b);
  Int_t lo = b, hi = b;
  while (lo > 1 && h->GetBinContent(lo) > 0.5 * c)
    lo--;
  while (hi < h->GetNbinsX() && h->GetBinContent(hi) > 0.5 * c)
    hi++;
  return TMath::Max(h->GetBinWidth(1),
                    (h->GetBinCenter(hi) - h->GetBinCenter(lo)) / 2.355);
}

// Height above the ridge: RegionCutFinder::AboveRidge, declared in the header
// so cross-section can rebuild a band. The reaction island lives at u ~ 4;
// the ridge itself, whatever its x, is at u ~ 0. This is the separation a
// plain y projection washes out.
using RegionCutFinder::AboveRidge;

// The fitted mixture evaluated per bin over the window, in a copy of the
// scatter, so it can be projected and contoured exactly like the data.
TH2F *ModelHist(TH2F *scatter, const RegionFit &fit, Int_t reac) {
  TH2F *m = static_cast<TH2F *>(scatter->Clone(Form("model_reac%d", reac)));
  m->SetDirectory(nullptr);
  m->Reset();
  TAxis *ax = m->GetXaxis(), *ay = m->GetYaxis();
  for (Int_t i = ax->FindBin(fit.x_lo); i <= ax->FindBin(fit.x_hi); i++)
    for (Int_t j = ay->FindBin(fit.y_lo); j <= ay->FindBin(fit.y_hi); j++) {
      Double_t x = ax->GetBinCenter(i), y = ay->GetBinCenter(j);
      m->SetBinContent(i, j,
                       Bigaus(fit.beam, x, y) +
                           (fit.has_reac ? Bigaus(fit.reac, x, y) : 0.0));
    }
  return m;
}

// The beam-like component from its core: the maximum and FWHM of each
// projection seed a +-2 sigma box, and the box's moments give the component.
// No minimiser: on a peak of 1e5 counts per bin Minuit's bigaus fit runs the
// amplitude to its bound and walks the mean out of the window, while the
// moments are exact. Returns kFALSE (with why) on an empty core.
Bool_t BeamFromCore(TH2F *scatter, Int_t reac, Int_t bx0, Int_t bx1, Int_t by0,
                    Int_t by1, Gauss2D &beam, TString &why) {
  TAxis *ax = scatter->GetXaxis(), *ay = scatter->GetYaxis();
  const Double_t bwx = ax->GetBinWidth(1), bwy = ay->GetBinWidth(1);
  TH1D *px = scatter->ProjectionX(Form("rcf_px_%d", reac), by0, by1);
  TH1D *py = scatter->ProjectionY(Form("rcf_py_%d", reac), bx0, bx1);
  px->SetDirectory(nullptr);
  py->SetDirectory(nullptr);
  px->GetXaxis()->SetRange(bx0, bx1);
  py->GetXaxis()->SetRange(by0, by1);
  Double_t mx0 = 0.0, my0 = 0.0;
  const Double_t sx0 = WidthAtMax(px, mx0), sy0 = WidthAtMax(py, my0);
  delete px;
  delete py;
  Double_t sw = 0, sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
  for (Int_t i = ax->FindBin(mx0 - 2 * sx0); i <= ax->FindBin(mx0 + 2 * sx0);
       i++)
    for (Int_t j = ay->FindBin(my0 - 2 * sy0); j <= ay->FindBin(my0 + 2 * sy0);
         j++) {
      Double_t w = scatter->GetBinContent(i, j);
      Double_t x = ax->GetBinCenter(i), y = ay->GetBinCenter(j);
      sw += w;
      sx += w * x;
      sy += w * y;
      sxx += w * x * x;
      syy += w * y * y;
      sxy += w * x * y;
    }
  if (!(sw > 0)) {
    why = "empty beam core";
    return kFALSE;
  }
  beam = MomentsToGauss(sw, sx, sy, sxx, syy, sxy, bwx, bwy);
  return kTRUE;
}

} // namespace

namespace RegionCutFinder {

RegionFit FitMixture(TH2F *scatter, Int_t reac, Double_t x_lo, Double_t x_hi,
                     Double_t y_lo, Double_t y_hi) {
  RegionFit fit;
  fit.x_lo = x_lo;
  fit.x_hi = x_hi;
  fit.y_lo = y_lo;
  fit.y_hi = y_hi;
  TAxis *ax = scatter->GetXaxis(), *ay = scatter->GetYaxis();
  const Int_t bx0 = ax->FindBin(x_lo), bx1 = ax->FindBin(x_hi);
  const Int_t by0 = ay->FindBin(y_lo), by1 = ay->FindBin(y_hi);
  const Double_t bwx = ax->GetBinWidth(1), bwy = ay->GetBinWidth(1);

  Gauss2D beam;
  if (!BeamFromCore(scatter, reac, bx0, bx1, by0, by1, beam, fit.why))
    return fit;

  // Reaction seed. The island is both above the ridge (u ~ 4) and beyond the
  // beam in x, while the beam's own tail above the ridge spreads over the
  // ridge's x range. So: in the band 3..15 conditional sigma above the ridge,
  // take the x maximum beyond +3 sigma_x of the beam, and seed at the centroid
  // of that neighbourhood. A local maximum in any 1D projection is not
  // required -- the real beam tail is heavier than Gaussian and the island
  // usually sits on a slope, not a bump.
  TH1D hx(Form("rcf_hx_%d", reac), "", 80, x_lo, x_hi);
  hx.SetDirectory(nullptr);
  for (Int_t i = bx0; i <= bx1; i++)
    for (Int_t j = by0; j <= by1; j++) {
      const Double_t w = scatter->GetBinContent(i, j);
      if (!(w > 0))
        continue;
      const Double_t x = ax->GetBinCenter(i), y = ay->GetBinCenter(j);
      const Double_t u = AboveRidge(beam, x, y);
      if (u > 3.0 && u < 15.0 && x > beam.mx + 3.0 * beam.sx)
        hx.Fill(x, w);
    }
  const Int_t bpk = hx.GetMaximumBin();
  const Double_t xpk = hx.GetBinCenter(bpk);
  if (hx.GetBinContent(bpk) < 5.0) {
    fit.why = "nothing above the ridge beyond +3 sigma_x of the beam";
    return fit;
  }
  Double_t s[6] = {0, 0, 0, 0, 0, 0};
  for (Int_t i = bx0; i <= bx1; i++)
    for (Int_t j = by0; j <= by1; j++) {
      const Double_t w = scatter->GetBinContent(i, j);
      if (!(w > 0))
        continue;
      const Double_t x = ax->GetBinCenter(i), y = ay->GetBinCenter(j);
      const Double_t u = AboveRidge(beam, x, y);
      if (u > 3.0 && u < 15.0 && std::fabs(x - xpk) < 2.0 * beam.sx) {
        s[0] += w;
        s[1] += w * x;
        s[2] += w * y;
        s[3] += w * x * x;
        s[4] += w * y * y;
        s[5] += w * x * y;
      }
    }
  if (s[0] < 20.0) {
    fit.why = Form("only %.0f events in the seed neighbourhood", s[0]);
    return fit;
  }
  Gauss2D reac_g = beam; // shape starts as the beam's
  reac_g.mx = s[1] / s[0];
  reac_g.my = s[2] / s[0];

  // Classification EM over the window, trimmed and capped: a bin joins the
  // reaction component only if it is above the ridge (u > 2), within 3 sigma
  // of the component, and the component's log-density with its prior beats
  // the beam's. Each component is then re-estimated from the moments of its
  // bins, the reaction's widths capped at twice the beam's. Without the trim
  // and cap the reaction component swallows the beam's non-Gaussian halo and
  // ends up as a wide blob on the ridge; the prior is what keeps a bin four
  // sigma down the beam tail with the beam.
  Double_t n_beam = 0.0, n_reac = 0.0;
  Double_t total = 0.0;
  for (Int_t i = bx0; i <= bx1; i++)
    for (Int_t j = by0; j <= by1; j++)
      total += scatter->GetBinContent(i, j);
  Double_t prior_reac = 1.0e-3; // until the first assignment counts it
  for (Int_t iter = 0; iter < 20; iter++) {
    Double_t mb[6] = {0, 0, 0, 0, 0, 0}, mr[6] = {0, 0, 0, 0, 0, 0};
    const Double_t lp_b = std::log(1.0 - prior_reac) - LogNorm(beam);
    const Double_t lp_r = std::log(prior_reac) - LogNorm(reac_g);
    for (Int_t i = bx0; i <= bx1; i++)
      for (Int_t j = by0; j <= by1; j++) {
        const Double_t w = scatter->GetBinContent(i, j);
        if (!(w > 0))
          continue;
        const Double_t x = ax->GetBinCenter(i), y = ay->GetBinCenter(j);
        const Double_t d2r = Mahal2(reac_g, x, y);
        // Admission starts where the seed search did, 3 sigma above the
        // ridge: letting bins in from 2 sigma handed the component the
        // beam's upper tail sheet, which it then tilted along the ridge to
        // absorb, and at strips where the island is only ~3 sigma_x off the
        // ridge that doubled the count with background.
        const Bool_t to_reac =
            AboveRidge(beam, x, y) > 3.0 && d2r < 9.0 &&
            (lp_r - 0.5 * d2r) > (lp_b - 0.5 * Mahal2(beam, x, y));
        Double_t *m = to_reac ? mr : mb;
        m[0] += w;
        m[1] += w * x;
        m[2] += w * y;
        m[3] += w * x * x;
        m[4] += w * y * y;
        m[5] += w * x * y;
      }
    if (mr[0] < 20.0) {
      fit.why = Form("reaction component starved (%.0f events) at iteration %d",
                     mr[0], iter);
      return fit;
    }
    Gauss2D nb =
        MomentsToGauss(mb[0], mb[1], mb[2], mb[3], mb[4], mb[5], bwx, bwy);
    Gauss2D nr =
        MomentsToGauss(mr[0], mr[1], mr[2], mr[3], mr[4], mr[5], bwx, bwy);
    // The island is the beam blob displaced by the reaction. In the
    // simulation, which has no background to confuse the measurement, it is
    // 4-6% wider than the beam in x and 8-10% in y at every strip: the
    // vertex position within the strip and the kinematics add that little.
    // So the cap is 1.1x. Anything looser lets the component take in beam
    // tail at strips where the island sits only ~3 sigma_x off the ridge
    // (at 1.5x, strip 5 held 931 events against 781 at the beam's width).
    nr.sx = TMath::Min(nr.sx, 1.1 * nb.sx);
    nr.sy = TMath::Min(nr.sy, 1.1 * nb.sy);
    // The island's x-y correlation comes from the same beam-energy spread as
    // the beam blob's, so it cannot exceed it: a reaction component more
    // elongated along the ridge than the beam is tracing the beam's tail.
    if (std::fabs(nr.rho) > std::fabs(nb.rho))
      nr.rho = nr.rho < 0 ? -std::fabs(nb.rho) : std::fabs(nb.rho);
    const Bool_t moved = std::fabs(nr.mx - reac_g.mx) > 0.1 * bwx ||
                         std::fabs(nr.my - reac_g.my) > 0.1 * bwy;
    beam = nb;
    reac_g = nr;
    n_beam = mb[0];
    n_reac = mr[0];
    prior_reac = TMath::Max(1.0e-6, n_reac / total);
    if (!moved)
      break;
  }
  fit.beam = beam;
  fit.reac = reac_g;

  if (fit.reac.my < fit.beam.my + 2.0 * fit.beam.sy) {
    fit.why = "reaction component fell back onto the ridge";
    return fit;
  }
  if (n_reac < 20.0) {
    fit.why = Form("reaction component holds only %.0f events", n_reac);
    return fit;
  }
  fit.n_beam = n_beam;
  fit.n_reac = n_reac;
  fit.ok = kTRUE;
  return fit;
}

Double_t AboveRidge(const Gauss2D &b, Double_t x, Double_t y) {
  return ((y - b.my) - b.rho * (b.sy / b.sx) * (x - b.mx)) /
         (b.sy * std::sqrt(1.0 - b.rho * b.rho));
}

RegionFit FitBeam(TH2F *scatter, Int_t reac, Double_t x_lo, Double_t x_hi,
                  Double_t y_lo, Double_t y_hi) {
  RegionFit fit;
  fit.x_lo = x_lo;
  fit.x_hi = x_hi;
  fit.y_lo = y_lo;
  fit.y_hi = y_hi;
  fit.has_reac = kFALSE;
  TAxis *ax = scatter->GetXaxis(), *ay = scatter->GetYaxis();
  if (!BeamFromCore(scatter, reac, ax->FindBin(x_lo), ax->FindBin(x_hi),
                    ay->FindBin(y_lo), ay->FindBin(y_hi), fit.beam, fit.why))
    return fit;
  fit.reac = fit.beam;
  fit.reac.amp = 0.0;
  fit.ok = kTRUE;
  return fit;
}

TCutG *RidgeBandCut(const char *name, const Gauss2D &beam, Double_t nsig_lo,
                    Double_t nsig_hi, Double_t x_lo, Double_t x_hi,
                    Double_t y_lo, Double_t y_hi) {
  // y on the ridge line at x, plus n conditional sigma; clipped to the window.
  auto edge = [&](Double_t x, Double_t n) {
    const Double_t y = beam.my +
                       beam.rho * (beam.sy / beam.sx) * (x - beam.mx) +
                       n * beam.sy * std::sqrt(1.0 - beam.rho * beam.rho);
    return TMath::Max(y_lo, TMath::Min(y_hi, y));
  };
  TCutG *c = new TCutG(name, 5);
  c->SetPoint(0, x_lo, edge(x_lo, nsig_lo));
  c->SetPoint(1, x_hi, edge(x_hi, nsig_lo));
  c->SetPoint(2, x_hi, edge(x_hi, nsig_hi));
  c->SetPoint(3, x_lo, edge(x_lo, nsig_hi));
  c->SetPoint(4, x_lo, edge(x_lo, nsig_lo));
  c->SetLineColor(kBlack);
  c->SetLineWidth(2);
  return c;
}

TCutG *EllipseCut(const char *name, const Gauss2D &g, Double_t nsigma,
                  Int_t npts) {
  // Mahalanobis contour at nsigma, via the Cholesky factor of the covariance.
  TCutG *c = new TCutG(name, npts + 1);
  const Double_t q = std::sqrt(1.0 - g.rho * g.rho);
  for (Int_t i = 0; i <= npts; i++) {
    Double_t t = 2.0 * TMath::Pi() * i / npts;
    Double_t u = std::cos(t), v = std::sin(t);
    c->SetPoint(i, g.mx + nsigma * g.sx * u,
                g.my + nsigma * g.sy * (g.rho * u + q * v));
  }
  c->SetLineColor(kBlack);
  c->SetLineWidth(2);
  return c;
}

Double_t CountInside(TH2F *scatter, TCutG *cut, const RegionFit &fit) {
  if (!cut)
    return 0.0;
  TAxis *ax = scatter->GetXaxis(), *ay = scatter->GetYaxis();
  Double_t n = 0.0;
  for (Int_t i = ax->FindBin(fit.x_lo); i <= ax->FindBin(fit.x_hi); i++)
    for (Int_t j = ay->FindBin(fit.y_lo); j <= ay->FindBin(fit.y_hi); j++)
      if (cut->IsInside(ax->GetBinCenter(i), ay->GetBinCenter(j)))
        n += scatter->GetBinContent(i, j);
  return n;
}

void SaveFigures(TH2F *scatter, Int_t reac, const RegionFit &fit, TCutG *an,
                 TCutG *aa, const TString &subdir) {
  scatter->GetXaxis()->SetRangeUser(fit.x_lo, fit.x_hi);
  scatter->GetYaxis()->SetRangeUser(fit.y_lo, fit.y_hi);
  const Int_t cReac = kRed + 1, cBeam = kAzure + 1;

  // 1. The regions on the scatter, with each fitted component's 1/2/3 sigma
  //    contours dashed behind them, so the fit and the cut it produced are
  //    both on the page.
  {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    c->SetLogz(kTRUE);
    PlottingUtils::ConfigureAndDraw2DHistogram(
        scatter, c,
        Form("reac %d: fitted components (dashed 1,2,3#sigma) and regions "
             "(filled)",
             reac));
    std::vector<TCutG *> tmp;
    for (Int_t k = 1; k <= 3; k++) {
      TCutG *eb = EllipseCut(Form("cb%d_%d", k, reac), fit.beam, k);
      eb->SetLineColor(cBeam);
      eb->SetLineStyle(2);
      eb->SetLineWidth(1);
      eb->Draw("L SAME");
      tmp.push_back(eb);
      if (!fit.has_reac)
        continue;
      TCutG *er = EllipseCut(Form("cr%d_%d", k, reac), fit.reac, k);
      er->SetLineColor(cReac);
      er->SetLineStyle(2);
      er->SetLineWidth(1);
      er->Draw("L SAME");
      tmp.push_back(er);
    }
    if (aa) {
      aa->SetFillColorAlpha(cBeam, 0.15);
      aa->SetLineColor(cBeam);
      aa->Draw("F SAME");
      aa->Draw("L SAME");
    }
    if (an) {
      an->SetFillColorAlpha(cReac, 0.25);
      an->SetLineColor(cReac);
      an->Draw("F SAME");
      an->Draw("L SAME");
    }
    // 2D: z is already log from the utility; kLOG means a log y here.
    PlottingUtils::SaveFigure(c, Form("regions_reac%d", reac), subdir,
                              PlotSaveOptions::kLINEAR);
    delete c;
    for (Int_t k = 0; k < Int_t(tmp.size()); k++)
      delete tmp[k];
  }

  // 2. The fitted density itself as contours over the data.
  TH2F *model = ModelHist(scatter, fit, reac);
  {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    c->SetLogz(kTRUE);
    PlottingUtils::ConfigureAndDraw2DHistogram(
        scatter, c,
        Form("reac %d: data with the fitted mixture's density", reac));
    model->SetContour(10);
    model->SetLineColor(cReac);
    model->SetLineWidth(1);
    model->Draw("CONT3 SAME");
    PlottingUtils::SaveFigure(c, Form("regions_reac%d_model", reac), subdir,
                              PlotSaveOptions::kLINEAR);
    delete c;
  }

  // 3./4. Each projection with the model's projection overlaid.
  TAxis *ax = scatter->GetXaxis(), *ay = scatter->GetYaxis();
  const Int_t bx0 = ax->FindBin(fit.x_lo), bx1 = ax->FindBin(fit.x_hi);
  const Int_t by0 = ay->FindBin(fit.y_lo), by1 = ay->FindBin(fit.y_hi);
  for (Int_t which = 0; which < 2; which++) {
    TH1D *d = which == 0 ? scatter->ProjectionX(Form("dpx_%d", reac), by0, by1)
                         : scatter->ProjectionY(Form("dpy_%d", reac), bx0, bx1);
    TH1D *m = which == 0 ? model->ProjectionX(Form("mpx_%d", reac), by0, by1)
                         : model->ProjectionY(Form("mpy_%d", reac), bx0, bx1);
    d->SetDirectory(nullptr);
    m->SetDirectory(nullptr);
    d->GetXaxis()->SetRange(which == 0 ? bx0 : by0, which == 0 ? bx1 : by1);
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kTRUE);
    PlottingUtils::ConfigureAndDrawHistogram(
        d, kBlack,
        Form("reac %d: %s projection, data (black) and fitted mixture (red)",
             reac, which == 0 ? "x" : "y"));
    m->SetLineColor(cReac);
    m->SetLineWidth(2);
    m->Draw("HIST SAME");
    PlottingUtils::SaveFigure(
        c, Form("regions_reac%d_proj%s", reac, which == 0 ? "x" : "y"), subdir,
        PlotSaveOptions::kLOG);
    delete c;
    delete d;
    delete m;
  }
  delete model;
}

} // namespace RegionCutFinder
