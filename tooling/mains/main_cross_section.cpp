// cross-section: absolute (a,xn) cross section per reaction strip.
//
// A read-only pass over the scatter cache and the region cuts. For each strip
//
//   sigma = N_reac / (N_beam * n_gas * L_strip)
//
// where N_reac is the number of tagged events inside that strip's (a,n)
// region, corrected for the fraction of the fitted component the region
// encloses, and N_beam is the number of beam particles that reached that strip
// under every cut a reaction there also had to pass, so those cuts cancel in
// the ratio.
//
// Beam energies come from the simulated unreacted beam, which is calibrated
// against the measured per-strip energy loss, so the centre-of-mass energy of
// each strip is the simulation's rather than a nominal dE/dx table's.
#include "BeamEnergies.hpp"
#include "Constants.hpp"
#include "InitUtils.hpp"
#include "Paths.hpp"
#include "PlottingUtils.hpp"
#include "RegionCuts.hpp"
#include "StripSumScatter.hpp"
#include <TAxis.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TKey.h>
#include <TLegend.h>
#include <TMath.h>
#include <TParameter.h>
#include <TROOT.h>
#include <TSystem.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

// Boltzmann constant [J/K] and Torr in Pa, for the ideal-gas number density.
const Double_t kBoltzmann = 1.380649e-23;
const Double_t kPaPerTorr = 133.322368;
// MUSIC's geometry and fill temperature, the same for every dataset: the
// active length of one anode strip along the beam, and room temperature.
const Double_t kStripLengthCm = 1.578;
const Double_t kGasTemperatureK = 293.0;
// A barn is 1e-24 cm^2, so a millibarn is 1e-27 cm^2.
const Double_t kCm2PerMb = 1.0e-27;

Long64_t ReadCount(TFile &f, const char *name, Bool_t &ok) {
  TParameter<Long64_t> *p = static_cast<TParameter<Long64_t> *>(f.Get(name));
  if (!p) {
    ok = kFALSE;
    return 0;
  }
  return p->GetVal();
}

// One TALYS model's (a,xn) excitation function: the sum of the
// residual-production graphs talys-xs wrote for every residue with
// Z = Z_beam + 2 and A <= A_beam + 3, i.e. all xn exits and no (a,gamma).
struct TalysCurve {
  TString label;
  TGraph *axn;
};

// Every model in root_files/talys/talys_xs.root, in the config's order.
// Empty when there is no file; the plot goes on without them.
std::vector<TalysCurve> LoadTalys() {
  const CrossSectionConfig &X = Constants::cfg.CROSS_SECTION_CONFIG;
  std::vector<TalysCurve> curves;
  const TString path = Paths::ResultsDir() + "/root_files/talys/talys_xs.root";
  if (gSystem->AccessPathName(path))
    return curves;
  TFile f(path, "READ");
  if (f.IsZombie())
    return curves;
  for (Int_t m = 0;; m++) {
    TDirectory *d = f.GetDirectory(Form("m%d", m));
    if (!d)
      break;
    std::map<Double_t, Double_t> sum;
    Int_t n_channels = 0;
    for (TIter it(d->GetListOfKeys()); TObject *k = it();) {
      Int_t z = 0, a = 0;
      if (sscanf(k->GetName(), "rp%3d%3d", &z, &a) != 2)
        continue;
      if (z != X.BEAM_Z + 2 || a > X.BEAM_A + 3)
        continue;
      TGraph *g = dynamic_cast<TGraph *>(d->Get(k->GetName()));
      if (!g)
        continue;
      for (Int_t p = 0; p < g->GetN(); p++)
        sum[g->GetX()[p]] += g->GetY()[p];
      n_channels++;
    }
    std::vector<Double_t> tx, ty;
    for (std::map<Double_t, Double_t>::const_iterator it = sum.begin();
         it != sum.end(); ++it)
      if (it->second > 0.0) {
        tx.push_back(it->first);
        ty.push_back(it->second);
      }
    if (tx.empty())
      continue;
    // The label is presentation, so the config's current one wins over the
    // one stamped in the file; relabelling never needs a TALYS rerun.
    TalysCurve c;
    TNamed *label = dynamic_cast<TNamed *>(d->Get("label"));
    c.label = m < Int_t(X.TALYS_MODELS.size()) ? X.TALYS_MODELS[m].label
              : label                          ? TString(label->GetTitle())
                      : TString(Form("TALYS model %d", m));
    c.axn = new TGraph(Int_t(tx.size()), &tx[0], &ty[0]);
    std::cout << "cross-section: " << c.label << ": (a,xn) from " << n_channels
              << " residual channels" << std::endl;
    curves.push_back(c);
  }
  return curves;
}

// The curve clipped to an energy window, for drawing.
TGraph *Clipped(TGraph *g, Double_t e_lo, Double_t e_hi) {
  std::vector<Double_t> x, y;
  for (Int_t p = 0; p < g->GetN(); p++)
    if (g->GetX()[p] >= e_lo && g->GetX()[p] <= e_hi) {
      x.push_back(g->GetX()[p]);
      y.push_back(g->GetY()[p]);
    }
  return x.empty() ? nullptr : new TGraph(Int_t(x.size()), &x[0], &y[0]);
}

// The effective energy of a strip spanning e_out..e_in: the energy at which
// the model's cross section equals its average over the strip, so a thin
// target at e_eff would give what the strip gives (Szegedi et al. 2021).
// The beam's dE/dx varies by a few percent within one strip, so the average
// is taken uniform in energy. Only the model's shape enters: a constant
// factor on sigma cancels. Bisection on a rising curve; the midpoint when
// there is no curve or the strip lies outside it.
Double_t EffectiveEnergy(TGraph *axn, Double_t e_out, Double_t e_in) {
  const Double_t mid = 0.5 * (e_out + e_in);
  if (!axn || axn->GetN() < 2)
    return mid;
  const Double_t lo = TMath::Min(e_out, e_in), hi = TMath::Max(e_out, e_in);
  if (lo < axn->GetX()[0] || hi > axn->GetX()[axn->GetN() - 1])
    return mid;
  const Int_t kSteps = 200;
  Double_t mean = 0.0;
  for (Int_t k = 0; k < kSteps; k++)
    mean += axn->Eval(lo + (k + 0.5) * (hi - lo) / kSteps);
  mean /= Double_t(kSteps);
  Double_t a = lo, b = hi;
  if ((axn->Eval(a) - mean) * (axn->Eval(b) - mean) > 0.0)
    return mid;
  for (Int_t k = 0; k < 40; k++) {
    const Double_t c = 0.5 * (a + b);
    if ((axn->Eval(a) - mean) * (axn->Eval(c) - mean) <= 0.0)
      b = c;
    else
      a = c;
  }
  return 0.5 * (a + b);
}

// Events in the scatter whose bin centre falls inside the cut, which is first
// scaled about its own centroid by `scale`. Scaling a region drawn at N sigma
// by k puts it at kN sigma, so counting at a few scales and correcting each by
// the fraction of a Gaussian it should enclose measures how much of what the
// region holds is not the component: a clean peak gives the same number at
// every scale, contamination grows with area.
// The cut scaled about its own centroid; caller owns the copy.
TCutG *ScaledCut(TCutG *cut, Double_t scale) {
  Double_t cx = 0.0, cy = 0.0;
  for (Int_t p = 0; p < cut->GetN(); p++) {
    cx += cut->GetX()[p];
    cy += cut->GetY()[p];
  }
  cx /= Double_t(cut->GetN());
  cy /= Double_t(cut->GetN());
  TCutG *scaled = new TCutG(*cut);
  for (Int_t p = 0; p < scaled->GetN(); p++)
    scaled->SetPoint(p, cx + scale * (cut->GetX()[p] - cx),
                     cy + scale * (cut->GetY()[p] - cy));
  return scaled;
}

Double_t CountInCut(TH2F *scatter, TCutG *cut, Double_t scale) {
  TCutG *scaled = ScaledCut(cut, scale);
  TAxis *ax = scatter->GetXaxis();
  TAxis *ay = scatter->GetYaxis();
  Double_t n = 0.0;
  for (Int_t i = 1; i <= scatter->GetNbinsX(); i++)
    for (Int_t j = 1; j <= scatter->GetNbinsY(); j++)
      if (scaled->IsInside(ax->GetBinCenter(i), ay->GetBinCenter(j)))
        n += scatter->GetBinContent(i, j);
  delete scaled;
  return n;
}

// Fraction of a bivariate Gaussian inside its n-sigma Mahalanobis ellipse.
Double_t Enclosed(Double_t nsigma) {
  return 1.0 - std::exp(-0.5 * nsigma * nsigma);
}

} // namespace

int main() {
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                Paths::ResultsDir() + "/plots",
                                Paths::ResultsDir() + "/root_files");
  gROOT->SetBatch(kTRUE);
  const CrossSectionConfig &X = Constants::cfg.CROSS_SECTION_CONFIG;
  const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;

  if (!(X.GAS_PRESSURE_TORR > 0.0) || X.BEAM_A <= 0) {
    std::cerr << "cross-section: this dataset has no CROSS_SECTION_CONFIG "
                 "(gas pressure and beam are both required)"
              << std::endl;
    return 1;
  }

  TString cache =
      IO::GetRootFilesBaseDir() + "/" + StripSumScatter::CacheName();
  TFile f(cache, "READ");
  if (f.IsZombie()) {
    std::cerr << "cross-section: no scatter cache at " << cache
              << "; run strip-sum-scatter first" << std::endl;
    return 1;
  }
  Bool_t have = kTRUE;
  const Long64_t n_seen = ReadCount(f, "n_seen", have);
  const Long64_t n_beam = ReadCount(f, "n_normed", have);
  if (!have || n_beam <= 0) {
    std::cerr << "cross-section: the cache carries no normalization counts, so "
                 "it predates them; re-run strip-sum-scatter to refill"
              << std::endl;
    return 1;
  }

  // Areal density of one strip: ideal gas at the fill conditions.
  const Double_t n_gas = X.GAS_PRESSURE_TORR * kPaPerTorr /
                         (kBoltzmann * kGasTemperatureK) * 1.0e-6 *
                         TargetGasAtomsPerMolecule(X.TARGET_GAS);
  const Double_t areal = n_gas * kStripLengthCm;

  Double_t dE[18];
  Double_t e_strip0 = 0.0;
  if (!BeamEnergies::Profile(BeamEnergies::SimPath(), dE, e_strip0)) {
    std::cerr << "cross-section: cannot read the simulated beam at "
              << BeamEnergies::SimPath() << std::endl;
    return 1;
  }
  const Double_t cm_frac = BeamEnergies::CmFraction();

  // Centre-of-mass energy at the middle of every strip, so the alignment to a
  // published table can be found over the whole detector rather than only the
  // strips a cross section is being reported for.
  Double_t e_mid[17];
  {
    Double_t e = e_strip0;
    for (Int_t s = 0; s <= 16; s++) {
      e_mid[s] = (e - 0.5 * dE[s]) * cm_frac;
      e -= dE[s];
    }
  }

  std::cout << "cross-section: " << n_beam << " beam events past every pre-tag "
            << "cut, of " << n_seen << " seen" << std::endl;
  std::cout << Form("  gas %.0f Torr at %.0f K -> %.4e atoms/cm^3; strip "
                    "%.3f cm -> %.4e atoms/cm^2",
                    X.GAS_PRESSURE_TORR, kGasTemperatureK, n_gas,
                    kStripLengthCm, areal)
            << std::endl;
  std::cout << Form("  regions at %.1f sigma enclose %.1f%% of the fitted "
                    "component; counts are corrected for that",
                    C.AN_REGION_NSIGMA,
                    100.0 * (1.0 - std::exp(-0.5 * C.AN_REGION_NSIGMA *
                                            C.AN_REGION_NSIGMA)))
            << std::endl;

  std::cout << std::endl;
  std::cout << "  strip   E_cm range [MeV]    E_cm,eff   N_reac    N_beam    "
               "sigma [mb]"
            << std::endl;

  // The TALYS models, loaded now because the first one's shape sets each
  // strip's effective energy.
  std::vector<TalysCurve> talys = LoadTalys();
  if (talys.empty())
    std::cout << "cross-section: no TALYS graphs in root_files (run talys-xs); "
                 "midpoint energies, no overlay"
              << std::endl;
  std::vector<Double_t> vx, vy, vexl, vexh, vey;
  std::vector<Int_t> vs;
  for (Int_t reac = X.XS_STRIP_MIN; reac <= X.XS_STRIP_MAX; reac++) {
    TH2F *h = static_cast<TH2F *>(f.Get(Form("scatter_r%d", reac)));
    TCutG *an = RegionCutStore::Load("region_an", reac);
    if (!h || !an) {
      std::cout << Form("   %2d    no %s", reac,
                        h ? "region cut" : "scatter in the cache")
                << std::endl;
      delete an;
      continue;
    }
    // Beam particles that reached this strip under the same conditions a
    // reaction here had to satisfy. Falls back to the flat count only if the
    // cache predates the per-strip one.
    Bool_t at_ok = kTRUE;
    const Long64_t n_at = ReadCount(f, Form("n_normed_r%d", reac), at_ok);
    const Long64_t n_denom = at_ok && n_at > 0 ? n_at : n_beam;

    // Beam energy at the strip's entrance and exit: what entered strip 0,
    // less each earlier strip's mean loss.
    Double_t e_lab = e_strip0;
    for (Int_t s = 0; s < reac; s++)
      e_lab -= dE[s];
    const Double_t e_cm_in = e_lab * cm_frac;
    const Double_t e_cm_out = (e_lab - dE[reac]) * cm_frac;
    // Effective energy from the first model's shape; the others say how much
    // it depends on the shape.
    const Double_t e_cm_eff = EffectiveEnergy(
        talys.empty() ? nullptr : talys[0].axn, e_cm_out, e_cm_in);
    Double_t e_eff_lo = e_cm_eff, e_eff_hi = e_cm_eff;
    for (Int_t m = 1; m < Int_t(talys.size()); m++) {
      const Double_t e = EffectiveEnergy(talys[m].axn, e_cm_out, e_cm_in);
      e_eff_lo = TMath::Min(e_eff_lo, e);
      e_eff_hi = TMath::Max(e_eff_hi, e);
    }
    if (talys.size() > 1)
      std::cout << Form(
                       "         E_cm,eff %.3f (midpoint %.3f); across the %zu "
                       "model shapes %.3f..%.3f -> energy systematic "
                       "+%.3f/-%.3f",
                       e_cm_eff, 0.5 * (e_cm_in + e_cm_out), talys.size(),
                       e_eff_lo, e_eff_hi, e_eff_hi - e_cm_eff,
                       e_cm_eff - e_eff_lo)
                << std::endl;

    // Two estimates of the reaction count. The geometric one is everything
    // inside the region, corrected for the fraction of a Gaussian it should
    // enclose; the attributed one is what the mixture fit assigned to the
    // reaction component (trimmed at 3 sigma, so corrected for that). Where
    // the island is well off the beam ridge they agree; where the region also
    // covers beam tail the geometric count runs away, since the beam's real
    // tail is far heavier than any Gaussian and cannot be subtracted by
    // model. The attributed count is what a person drawing a tight region
    // gets, so it is the value; the two estimates' disagreement is the region
    // systematic, because that disagreement is exactly the overlap ambiguity.
    const Double_t n_raw = CountInCut(h, an, 1.0);
    const Double_t est_geo = n_raw / Enclosed(C.AN_REGION_NSIGMA);
    const Double_t n_assigned = RegionCutStore::LoadAssigned(reac);
    const Bool_t have_fit = n_assigned >= 0.0;
    const Double_t est_fit = have_fit ? n_assigned / Enclosed(3.0) : est_geo;
    const Double_t n_reac = est_fit;
    const Double_t norm = Double_t(n_denom) * areal * kCm2PerMb;
    const Double_t sigma = n_reac / norm;
    const Double_t sigma_stat = n_reac > 0.0 ? sigma / std::sqrt(n_reac) : 0.0;
    // Without an attributed count (a drawn region) fall back to scaling the
    // region by 0.75x and 1.25x, each corrected for its enclosed fraction.
    Double_t sigma_sys = 0.0;
    if (have_fit) {
      // Against the region's core: the 1-sigma ellipse (half the drawn
      // region), corrected for the 39% it encloses. The full region's count
      // is no check where it also covers tail, but the core is where the
      // island dominates, so core and attributed counts should agree, and
      // their disagreement is the honest size of the overlap ambiguity.
      const Double_t est_core =
          CountInCut(h, an, 1.0 / C.AN_REGION_NSIGMA) / Enclosed(1.0);
      sigma_sys = std::fabs(est_core - est_fit) / norm;
    } else {
      Double_t lo = sigma, hi = sigma;
      const Double_t kScale[2] = {0.75, 1.25};
      for (Int_t k = 0; k < 2; k++) {
        const Double_t s = CountInCut(h, an, kScale[k]) /
                           Enclosed(kScale[k] * C.AN_REGION_NSIGMA) / norm;
        lo = TMath::Min(lo, s);
        hi = TMath::Max(hi, s);
      }
      sigma_sys = 0.5 * (hi - lo);
    }
    const Double_t sigma_err =
        std::sqrt(sigma_stat * sigma_stat + sigma_sys * sigma_sys);

    std::cout << Form("   %2d    [%5.2f, %5.2f]      %6.2f   %6.0f  %9lld   "
                      "%7.1f +- %.1f (%.1f stat, %.1f region; in-region "
                      "%.0f, attributed %s)",
                      reac, e_cm_in, e_cm_out, e_cm_eff, n_reac, n_denom, sigma,
                      sigma_err, sigma_stat, sigma_sys, est_geo,
                      have_fit ? Form("%.0f", est_fit) : "none")
              << std::endl;
    // The energy error is the strip's extent about the effective energy,
    // asymmetric since that energy sits above the midpoint on a rising
    // curve (the reference table's convention), with the effective energy's
    // dependence on the model shape -- how far the other models move it --
    // added in quadrature on each side as a systematic.
    const Double_t e_sys_lo = e_cm_eff - e_eff_lo,
                   e_sys_hi = e_eff_hi - e_cm_eff;
    vs.push_back(reac);
    vx.push_back(e_cm_eff);
    vy.push_back(sigma);
    vexl.push_back(
        std::hypot(e_cm_eff - TMath::Min(e_cm_in, e_cm_out), e_sys_lo));
    vexh.push_back(
        std::hypot(TMath::Max(e_cm_in, e_cm_out) - e_cm_eff, e_sys_hi));
    vey.push_back(sigma_err);
    delete an;
  }

  if (vx.empty()) {
    std::cerr << "cross-section: no strip produced a cross section"
              << std::endl;
    return 1;
  }

  // Against the published values. A published table is a list of strips, so
  // the two are aligned strip to row by the single integer offset that best
  // matches the energies over the whole table, not row by row: pairing each
  // strip with its nearest published energy silently pairs two strips with
  // one row when a dataset reaches energies the reference never reported.
  if (!X.REFERENCE_XS.empty()) {
    const Int_t n_ref = Int_t(X.REFERENCE_XS.size());
    Int_t best_off = 0;
    Double_t best_rms = 1.0e9;
    for (Int_t off = 0; off + n_ref - 1 <= 16; off++) {
      Double_t s = 0.0;
      Int_t used = 0;
      for (Int_t k = 0; k < n_ref; k++) {
        const Int_t strip = off + k;
        if (strip < 0 || strip > 16)
          continue;
        const Double_t d = e_mid[strip] - X.REFERENCE_XS[k][0];
        s += d * d;
        used++;
      }
      if (used < n_ref)
        continue;
      const Double_t rms = std::sqrt(s / Double_t(used));
      if (rms < best_rms) {
        best_rms = rms;
        best_off = off;
      }
    }
    std::cout << std::endl;
    std::cout << Form("  vs %s: first row is strip %d from the energies, rms "
                      "%.3f MeV",
                      X.REFERENCE_LABEL.Data(), best_off, best_rms)
              << std::endl;
    std::cout << "  strip   E_cm,eff   published E   this work [mb]      "
                 "published [mb]   ratio"
              << std::endl;
    for (Int_t reac = X.XS_STRIP_MIN; reac <= X.XS_STRIP_MAX; reac++) {
      const Int_t row = reac - best_off;
      Int_t idx = -1;
      for (Int_t i = 0; i < Int_t(vs.size()); i++)
        if (vs[i] == reac)
          idx = i;
      if (idx < 0)
        continue;
      if (row < 0 || row >= n_ref) {
        std::cout << Form("   %2d     %6.2f      --          %7.1f          "
                          "   -- (upstream of the published range)",
                          reac, vx[idx], vy[idx])
                  << std::endl;
        continue;
      }
      const Double_t ref = X.REFERENCE_XS[row][3];
      std::cout << Form("   %2d     %6.2f    %6.2f        %7.1f +- %-5.1f  "
                        "%7.1f (%.1f)   %5.2f",
                        reac, vx[idx], X.REFERENCE_XS[row][0], vy[idx],
                        vey[idx], ref, X.REFERENCE_XS[row][4],
                        ref > 0.0 ? vy[idx] / ref : 0.0)
                << std::endl;
    }
  }

  // Excitation function, this work against the published points. The frame
  // spans both sets so the measured strips are seen in the context of the
  // whole published curve, not just the rows they happen to sit beside.
  Double_t fx_lo = vx[0], fx_hi = vx[0], fy_lo = vy[0], fy_hi = vy[0];
  for (Int_t i = 0; i < Int_t(vx.size()); i++) {
    fx_lo = TMath::Min(fx_lo, vx[i] - vexl[i]);
    fx_hi = TMath::Max(fx_hi, vx[i] + vexh[i]);
    fy_lo = TMath::Min(fy_lo, vy[i]);
    fy_hi = TMath::Max(fy_hi, vy[i]);
  }
  for (Int_t k = 0; k < Int_t(X.REFERENCE_XS.size()); k++) {
    fx_lo = TMath::Min(fx_lo, X.REFERENCE_XS[k][0]);
    fx_hi = TMath::Max(fx_hi, X.REFERENCE_XS[k][0]);
    fy_lo = TMath::Min(fy_lo, X.REFERENCE_XS[k][3]);
    fy_hi = TMath::Max(fy_hi, X.REFERENCE_XS[k][3]);
  }
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kTRUE);
  TH1F *frame =
      c->DrawFrame(fx_lo - 0.4, 0.5 * fy_lo, fx_hi + 0.4, 3.0 * fy_hi);
  frame->SetTitle(Form("%s(#alpha, xn);E_{c.m.,eff} [MeV];#sigma [mb]",
                       Paths::DatasetName().Data()));
  // Hauser-Feshbach prediction, if talys-xs has written one: drawn first so
  // the measured points sit on top of it, and only over the frame's energies.
  // The first model in full, the shape checks dashed.
  std::vector<TGraph *> hf;
  for (Int_t m = 0; m < Int_t(talys.size()); m++) {
    TGraph *c = Clipped(talys[m].axn, fx_lo - 0.4, fx_hi + 0.4);
    if (!c)
      continue;
    c->SetLineColor(kAzure + 1);
    c->SetLineWidth(m == 0 ? 2 : 1);
    c->SetLineStyle(m == 0 ? 1 : 2);
    c->Draw("L SAME");
    hf.push_back(c);
  }
  TGraphAsymmErrors *g = new TGraphAsymmErrors(
      Int_t(vx.size()), &vx[0], &vy[0], &vexl[0], &vexh[0], &vey[0], &vey[0]);
  g->SetMarkerStyle(20);
  g->SetMarkerColor(kBlack);
  g->SetLineColor(kBlack);
  // The reference's energy errors are the strip's extent about the effective
  // energy, and asymmetric where the table gives them so.
  TGraphAsymmErrors *r = nullptr;
  if (!X.REFERENCE_XS.empty()) {
    std::vector<Double_t> rx, ry, rexl, rexh, rey;
    for (Int_t k = 0; k < Int_t(X.REFERENCE_XS.size()); k++) {
      const std::vector<Double_t> &row = X.REFERENCE_XS[k];
      rx.push_back(row[0]);
      rexh.push_back(row[1]);
      rexl.push_back(row[2]);
      ry.push_back(row[3]);
      rey.push_back(row[4]);
    }
    r = new TGraphAsymmErrors(Int_t(rx.size()), &rx[0], &ry[0], &rexl[0],
                              &rexh[0], &rey[0], &rey[0]);
    r->SetMarkerStyle(24);
    r->SetMarkerColor(kRed + 1);
    r->SetLineColor(kRed + 1);
    r->Draw("P SAME");
  }
  g->Draw("P SAME");
  // Bottom right is the one empty corner: the excitation function climbs to
  // the upper right and the reference table starts at the lower left.
  TLegend *leg = PlottingUtils::AddLegend(0.42, 0.89, 0.16,
                                          0.16 + 0.07 * (2 + Int_t(hf.size())));
  leg->AddEntry(g, "Present Work", "pe");
  if (r)
    leg->AddEntry(r, X.REFERENCE_LABEL, "pe");
  for (Int_t m = 0; m < Int_t(hf.size()); m++)
    leg->AddEntry(hf[m], talys[m].label, "l");
  leg->Draw();
  PlottingUtils::SaveFigure(c, "cross_section", "cross_section",
                            PlotSaveOptions::kLOG);
  delete c;
  f.Close();
  return 0;
}
