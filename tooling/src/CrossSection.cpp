#include "CrossSection.hpp"
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
#include <TFeldmanCousins.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TNamed.h>
#include <TPad.h>
#include <TParameter.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVectorD.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>

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

// A tagged epoch's figures carry its tag, beside the untagged eras'.
TString TagSuffix() {
  const TString &tag = Constants::ActiveFileTag();
  return tag.Length() > 0 ? "_" + tag : TString("");
}

// The curve's value at e, interpolated between the neighbouring grid points:
// log-linear where both are positive (as the curve is drawn), else linear.
Bool_t CurveAt(const TGraph *g, Double_t e, Double_t &v) {
  for (Int_t p = 1; p < g->GetN(); p++) {
    const Double_t x0 = g->GetX()[p - 1], x1 = g->GetX()[p];
    if (!(x0 < e && e < x1))
      continue;
    const Double_t y0 = g->GetY()[p - 1], y1 = g->GetY()[p];
    const Double_t t = (e - x0) / (x1 - x0);
    v = (y0 > 0.0 && y1 > 0.0)
            ? std::exp(std::log(y0) + t * (std::log(y1) - std::log(y0)))
            : y0 + t * (y1 - y0);
    return kTRUE;
  }
  return kFALSE;
}

// Marker styles and colours for the channels on a combined figure.
const Int_t kChannelMarker[4] = {20, 21, 22, 23};
const Int_t kChannelColor[4] = {kBlack, kGreen + 2, kMagenta + 2, kOrange + 7};
// Colours for the model curves on a single-channel figure, one per TALYS
// model in TALYS_MODELS order.
const Int_t kModelColor[8] = {kAzure + 1,   kRed + 1,  kGreen + 2, kOrange + 7,
                              kMagenta + 1, kCyan + 2, kGray + 2,  kViolet - 3};
// Colours for the per-exit components on the fit figure, in exit order.
const Int_t kComponentColor[4] = {kBlack, kBlue + 1, kGreen + 2, kMagenta + 1};

} // namespace

Double_t CrossSection::Point::ErrLo() const {
  return std::sqrt(stat_lo * stat_lo + sys * sys);
}

Double_t CrossSection::Point::ErrHi() const {
  return std::sqrt(stat_hi * stat_hi + sys * sys);
}

void CrossSection::PoissonInterval(Double_t count, Double_t &below,
                                   Double_t &above) {
  const Int_t n = Int_t(std::lround(TMath::Max(count, 0.0)));
  if (n >= Constants::cfg.CROSS_SECTION_CONFIG.FELDMAN_COUSINS_MAX_COUNT) {
    below = above = std::sqrt(count);
    return;
  }
  // Feldman & Cousins, PRD 57 (1998) 3873, no background: the interval on
  // mu whose likelihood-ratio-ordered 68.27% acceptance region holds n.
  TFeldmanCousins fc(0.6827);
  fc.SetMuMax(n + 10.0 * std::sqrt(n + 1.0) + 10.0);
  fc.SetMuStep(0.005);
  const Double_t hi = fc.CalculateUpperLimit(n, 0.0);
  const Double_t lo = fc.GetLowerLimit();
  below = n - lo;
  above = hi - n;
}

/// "n", "2n", "pn", "a", "g": light particles leaving the compound nucleus,
/// each with an optional multiplicity digit in front. The residue is the
/// compound (beam + alpha) minus what left.
Bool_t CrossSection::ExitResidue(const TString &exit, Int_t z_beam,
                                 Int_t a_beam, Int_t &z, Int_t &a) {
  Int_t dz = 0, da = 0, mult = 0;
  for (Int_t i = 0; i < exit.Length(); i++) {
    const char c = exit[i];
    if (c >= '0' && c <= '9') {
      mult = mult * 10 + (c - '0');
      continue;
    }
    Int_t pz = 0, pa = 0;
    switch (c) {
    case 'n':
      pz = 0;
      pa = 1;
      break;
    case 'p':
      pz = 1;
      pa = 1;
      break;
    case 'd':
      pz = 1;
      pa = 2;
      break;
    case 't':
      pz = 1;
      pa = 3;
      break;
    case 'h':
      pz = 2;
      pa = 3;
      break;
    case 'a':
      pz = 2;
      pa = 4;
      break;
    case 'g':
      pz = 0;
      pa = 0;
      break;
    default:
      return kFALSE;
    }
    const Int_t m = mult > 0 ? mult : 1;
    dz += m * pz;
    da += m * pa;
    mult = 0;
  }
  if (mult > 0 || exit.Length() == 0)
    return kFALSE;
  z = z_beam + 2 - dz;
  a = a_beam + 4 - da;
  return z > 0 && a > z;
}

TString CrossSection::Label(const CrossSectionChannel &ch) {
  if (ch.label.Length() > 0)
    return ch.label;
  if (ch.talys_exits.empty())
    return "";
  Bool_t all_neutron = kTRUE;
  for (Int_t k = 0; k < Int_t(ch.talys_exits.size()); k++) {
    TString e = ch.talys_exits[k];
    e.ReplaceAll("n", "");
    for (Int_t i = 0; i < e.Length(); i++)
      if (e[i] < '0' || e[i] > '9')
        all_neutron = kFALSE;
  }
  if (ch.talys_exits.size() == 1)
    return "(#alpha, " + ch.talys_exits[0] + ")";
  if (all_neutron)
    return "(#alpha, xn)";
  TString joined;
  for (Int_t k = 0; k < Int_t(ch.talys_exits.size()); k++)
    joined += (k ? ", " : "") + ch.talys_exits[k];
  return "(#alpha, " + joined + ")";
}

TString CrossSection::SubtractedLabel(const CrossSectionChannel &ch) {
  if (ch.subtract_exits.empty())
    return "";
  TString joined;
  for (Int_t k = 0; k < Int_t(ch.subtract_exits.size()); k++)
    joined += (k ? ", " : "") + ch.subtract_exits[k];
  return "(#alpha, " + joined + ")";
}

Long64_t CrossSection::ReadCount(TFile &f, const char *name, Bool_t &ok) {
  TParameter<Long64_t> *p = static_cast<TParameter<Long64_t> *>(f.Get(name));
  if (!p) {
    ok = kFALSE;
    return 0;
  }
  return p->GetVal();
}

// The curve at e, log-linear between grid points where both are positive
// (as it is drawn), zero outside the grid.
Double_t CrossSection::Interpolated(TGraph *g, Double_t e) {
  if (!g || g->GetN() < 2)
    return 0.0;
  for (Int_t p = 1; p < g->GetN(); p++) {
    const Double_t x0 = g->GetX()[p - 1], x1 = g->GetX()[p];
    if (!(x0 <= e && e <= x1))
      continue;
    const Double_t y0 = g->GetY()[p - 1], y1 = g->GetY()[p];
    const Double_t t = (e - x0) / (x1 - x0);
    return (y0 > 0.0 && y1 > 0.0)
               ? std::exp(std::log(y0) + t * (std::log(y1) - std::log(y0)))
               : y0 + t * (y1 - y0);
  }
  return 0.0;
}

// The curve clipped to an energy window, for drawing.
TGraph *CrossSection::Clipped(TGraph *g, Double_t e_lo, Double_t e_hi) {
  // The grid points inside the window, carried to its edges by
  // interpolation so the curve reaches the frame.
  std::vector<Double_t> x, y;
  Double_t v = 0.0;
  if (CurveAt(g, e_lo, v)) {
    x.push_back(e_lo);
    y.push_back(v);
  }
  for (Int_t p = 0; p < g->GetN(); p++)
    if (g->GetX()[p] >= e_lo && g->GetX()[p] <= e_hi) {
      x.push_back(g->GetX()[p]);
      y.push_back(g->GetY()[p]);
    }
  if (CurveAt(g, e_hi, v)) {
    x.push_back(e_hi);
    y.push_back(v);
  }
  return x.empty() ? nullptr : new TGraph(Int_t(x.size()), &x[0], &y[0]);
}

/// The effective energy of a strip spanning e_out..e_in: the energy below
/// which half of the strip's yield is produced, with the model's cross
/// section taken linear between its values at the strip entrance (sigma_1)
/// and exit (sigma_2) and the beam's dE/dx uniform across the strip,
///
///   E_eff = E_0 - dE + dE [ -s2/(s1 - s2) + sqrt((s1^2 + s2^2) / (2 (s1 -
///   s2)^2)) ]
///
/// a good approximation for sigma_1/sigma_2 up to about 10. Only the model's
/// shape enters: a constant factor on sigma cancels. Evaluated in the
/// equivalent form (s1 + s2) / (2 (sqrt((s1^2 + s2^2)/2) + s2)) for the
/// fraction of dE, which has no 0/0 at sigma_1 = sigma_2 and gives the
/// midpoint there. The midpoint when there is no curve, the strip lies
/// outside it, or the model gives no positive cross section at either edge.
Double_t CrossSection::EffectiveEnergy(TGraph *axn, Double_t e_out,
                                       Double_t e_in) {
  const Double_t mid = 0.5 * (e_out + e_in);
  if (!axn || axn->GetN() < 2)
    return mid;
  const Double_t lo = TMath::Min(e_out, e_in), hi = TMath::Max(e_out, e_in);
  if (lo < axn->GetX()[0] || hi > axn->GetX()[axn->GetN() - 1])
    return mid;
  const Double_t s1 = axn->Eval(hi), s2 = axn->Eval(lo);
  if (s1 < 0.0 || s2 < 0.0 || s1 + s2 <= 0.0)
    return mid;
  const Double_t frac =
      (s1 + s2) / (2.0 * (TMath::Sqrt(0.5 * (s1 * s1 + s2 * s2)) + s2));
  return lo + frac * (hi - lo);
}

// The cut scaled about its own centroid; caller owns the copy. Scaling a
// region drawn at N sigma by k puts it at kN sigma.
TCutG *CrossSection::ScaledCut(TCutG *cut, Double_t scale) {
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

/// Events in the scatter whose bin centre falls inside the cut scaled by
/// `scale`. Counting at a few scales and correcting each by the fraction of a
/// Gaussian it should enclose measures how much of what the region holds is
/// not the component: a clean peak gives the same number at every scale,
/// contamination grows with area.
Double_t CrossSection::CountInCut(TH2F *scatter, TCutG *cut, Double_t scale) {
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
Double_t CrossSection::Enclosed(Double_t nsigma) {
  return 1.0 - std::exp(-0.5 * nsigma * nsigma);
}

Double_t CrossSection::GasPressureSys(Double_t sigma) const {
  const CrossSectionConfig &X = Constants::cfg.CROSS_SECTION_CONFIG;
  if (!(X.GAS_PRESSURE_TORR_ERR > 0.0) || !(X.GAS_PRESSURE_TORR > 0.0))
    return 0.0;
  return sigma * X.GAS_PRESSURE_TORR_ERR / X.GAS_PRESSURE_TORR;
}

// The cut-variation systematic: per threshold the larger change of the
// count under its up and down shift, in quadrature; `detail` lists them.
Double_t CrossSection::CutVariation(Int_t reac, Double_t per_count,
                                    TString &detail) const {
  detail = "";
  Bool_t at_ok = kTRUE;
  const Long64_t n_nom = ReadCount(*cache_, Form("n_tagged_r%d", reac), at_ok);
  Bool_t d_ok = kTRUE;
  const Long64_t d_nom = ReadCount(*cache_, Form("n_normed_r%d", reac), d_ok);
  TNamed *vn = static_cast<TNamed *>(cache_->Get("cut_variants"));
  if (!at_ok || !d_ok || d_nom <= 0 || !vn)
    return 0.0;
  const Double_t s_nom = Double_t(n_nom) / Double_t(d_nom) * per_count;
  std::map<TString, Double_t> worst;
  std::map<TString, TString> sides;
  TString names = vn->GetTitle(), tok;
  Int_t from = 0;
  while (names.Tokenize(tok, from, ",")) {
    if (tok.IsNull())
      continue;
    Bool_t ok = kTRUE;
    const Long64_t n_v =
        ReadCount(*cache_, Form("n_tagged_r%d_%s", reac, tok.Data()), ok);
    if (!ok)
      continue;
    // The variant's own denominator; a cache written before the variants
    // carried one falls back to the nominal.
    Bool_t dv_ok = kTRUE;
    const Long64_t d_v =
        ReadCount(*cache_, Form("n_normed_r%d_%s", reac, tok.Data()), dv_ok);
    const Double_t denom = dv_ok && d_v > 0 ? Double_t(d_v) : Double_t(d_nom);
    const TString base = tok(0, tok.Length() - 1);
    const Double_t ds = Double_t(n_v) / denom * per_count - s_nom;
    const Double_t d = std::fabs(ds);
    if (d > worst[base])
      worst[base] = d;
    sides[base] += Form("%s%s %+.1f", sides[base].IsNull() ? "" : ", ",
                        tok.EndsWith("+") ? "up" : "down", ds);
  }
  Double_t var2 = 0.0;
  for (std::map<TString, Double_t>::const_iterator it = worst.begin();
       it != worst.end(); ++it) {
    var2 += it->second * it->second;
    detail += Form(" %s[%s]", it->first.Data(), sides[it->first].Data());
  }
  return std::sqrt(var2);
}

Bool_t CrossSection::LoadCache() {
  const CrossSectionConfig &X = Constants::cfg.CROSS_SECTION_CONFIG;
  TString path = IO::GetRootFilesBaseDir() + "/" + StripSumScatter::CacheName();
  cache_ = new TFile(path, "READ");
  if (cache_->IsZombie()) {
    std::cerr << "cross-section: no scatter cache at " << path
              << "; run strip-sum-scatter first" << std::endl;
    return kFALSE;
  }
  Bool_t have = kTRUE;
  n_seen_ = ReadCount(*cache_, "n_seen", have);
  n_beam_ = ReadCount(*cache_, "n_normed", have);
  if (!have || n_beam_ <= 0) {
    std::cerr << "cross-section: the cache carries no normalization counts, so "
                 "it predates them; re-run strip-sum-scatter to refill"
              << std::endl;
    return kFALSE;
  }
  // Areal density of one strip: ideal gas at the fill conditions.
  n_gas_ = X.GAS_PRESSURE_TORR * kPaPerTorr / (kBoltzmann * kGasTemperatureK) *
           1.0e-6 * TargetGasAtomsPerMolecule(X.TARGET_GAS);
  areal_ = n_gas_ * kStripLengthCm;
  return kTRUE;
}

/// The simulated beam's loss per strip and the energy entering strip 0, and
/// from them the centre-of-mass energy at the middle of every strip, so the
/// alignment to a published table can be found over the whole detector rather
/// than only the strips a cross section is reported for.
Bool_t CrossSection::LoadBeam() {
  if (!BeamEnergies::Profile(BeamEnergies::SimPath(), dE_, e_strip0_)) {
    std::cerr << "cross-section: cannot read the simulated beam at "
              << BeamEnergies::SimPath() << std::endl;
    return kFALSE;
  }
  cm_frac_ = BeamEnergies::CmFraction();
  Double_t e = e_strip0_;
  for (Int_t s = 0; s <= 16; s++) {
    e_mid_[s] = (e - 0.5 * dE_[s]) * cm_frac_;
    e -= dE_[s];
  }
  return kTRUE;
}

/// Every model in root_files/talys/talys_xs.root, in the config's order, with
/// every residual-production graph talys-xs wrote for it. Empty when there is
/// no file; the plot goes on without them. The label is presentation, so the
/// config's current one wins over the one stamped in the file; relabelling
/// never needs a TALYS rerun.
void CrossSection::LoadTalys() {
  const CrossSectionConfig &X = Constants::cfg.CROSS_SECTION_CONFIG;
  const TString path = Paths::ResultsDir() + "/root_files/talys/talys_xs.root";
  if (gSystem->AccessPathName(path))
    return;
  TFile f(path, "READ");
  if (f.IsZombie())
    return;
  for (Int_t m = 0;; m++) {
    TDirectory *d = f.GetDirectory(Form("m%d", m));
    if (!d)
      break;
    std::map<std::pair<Int_t, Int_t>, TGraph *> graphs;
    for (TIter it(d->GetListOfKeys()); TObject *k = it();) {
      Int_t z = 0, a = 0;
      if (sscanf(k->GetName(), "rp%3d%3d", &z, &a) != 2)
        continue;
      TGraph *g = dynamic_cast<TGraph *>(d->Get(k->GetName()));
      if (g)
        graphs[std::make_pair(z, a)] = static_cast<TGraph *>(g->Clone());
    }
    TNamed *label = dynamic_cast<TNamed *>(d->Get("label"));
    talys_labels_.push_back(m < Int_t(X.TALYS_MODELS.size())
                                ? X.TALYS_MODELS[m].label
                            : label ? TString(label->GetTitle())
                                    : TString(Form("TALYS model %d", m)));
    talys_raw_.push_back(graphs);
  }
}

std::vector<CrossSection::TalysCurve>
CrossSection::ExitCurves(const CrossSectionChannel &ch,
                         const std::vector<TString> &exits,
                         const char *what) const {
  const CrossSectionConfig &X = Constants::cfg.CROSS_SECTION_CONFIG;
  std::vector<TalysCurve> curves;
  if (talys_raw_.empty() || exits.empty())
    return curves;
  for (Int_t m = 0; m < Int_t(talys_raw_.size()); m++) {
    std::map<Double_t, Double_t> sum;
    TString used, missing;
    for (Int_t k = 0; k < Int_t(exits.size()); k++) {
      Int_t z = 0, a = 0;
      ExitResidue(exits[k], X.BEAM_Z, X.BEAM_A, z, a);
      std::map<std::pair<Int_t, Int_t>, TGraph *>::const_iterator it =
          talys_raw_[m].find(std::make_pair(z, a));
      if (it == talys_raw_[m].end()) {
        missing += Form(" %s(Z=%d,A=%d)", exits[k].Data(), z, a);
        continue;
      }
      for (Int_t p = 0; p < it->second->GetN(); p++)
        sum[it->second->GetX()[p]] += it->second->GetY()[p];
      used += Form(" %s(Z=%d,A=%d)", exits[k].Data(), z, a);
    }
    std::vector<Double_t> tx, ty;
    for (std::map<Double_t, Double_t>::const_iterator it = sum.begin();
         it != sum.end(); ++it)
      if (it->second > 0.0) {
        tx.push_back(it->first);
        ty.push_back(it->second);
      }
    std::cout << "cross-section: " << ch.name << ": " << talys_labels_[m]
              << ": summed " << what << (used.Length() ? used : " nothing")
              << (missing.Length() ? "; not in the file:" + missing : "")
              << std::endl;
    if (tx.empty())
      continue;
    TalysCurve c;
    c.label = talys_labels_[m];
    c.axn = new TGraph(Int_t(tx.size()), &tx[0], &ty[0]);
    curves.push_back(c);
  }
  return curves;
}

Double_t CrossSection::StripAverage(TGraph *g, Double_t e_out, Double_t e_in,
                                    Bool_t &covered) {
  const Double_t lo = TMath::Min(e_out, e_in), hi = TMath::Max(e_out, e_in);
  covered = g && g->GetN() >= 2 && lo >= g->GetX()[0] &&
            hi <= g->GetX()[g->GetN() - 1];
  const Int_t n_sub = 200;
  Double_t sum = 0.0;
  for (Int_t k = 0; k <= n_sub; k++) {
    const Double_t e = lo + (hi - lo) * Double_t(k) / Double_t(n_sub);
    sum += ((k == 0 || k == n_sub) ? 0.5 : 1.0) * Interpolated(g, e);
  }
  return sum / Double_t(n_sub);
}

Double_t
CrossSection::SubtractedFraction(const std::vector<TalysCurve> &talys,
                                 const std::vector<TalysCurve> &talys_sub,
                                 Double_t e_out, Double_t e_in,
                                 Double_t &spread) {
  spread = 0.0;
  if (talys_sub.empty() || talys.empty())
    return 0.0;
  Double_t sum = 0.0, lo = 1.0e9, hi = -1.0e9;
  Int_t used = 0, skipped = 0;
  // Pair kept and subtracted curves by model label.
  for (Int_t ms = 0; ms < Int_t(talys_sub.size()); ms++) {
    const TalysCurve *chan = nullptr;
    for (Int_t m = 0; m < Int_t(talys.size()); m++)
      if (talys[m].label == talys_sub[ms].label) {
        chan = &talys[m];
        break;
      }
    if (!chan) {
      skipped++;
      continue;
    }
    Bool_t cov_chan = kFALSE, cov_drop = kFALSE;
    const Double_t keep = StripAverage(chan->axn, e_out, e_in, cov_chan);
    const Double_t drop =
        StripAverage(talys_sub[ms].axn, e_out, e_in, cov_drop);
    if (!cov_chan || !cov_drop) {
      skipped++;
      continue;
    }
    if (!(keep + drop > 0.0))
      continue;
    const Double_t f = drop / (keep + drop);
    sum += f;
    lo = TMath::Min(lo, f);
    hi = TMath::Max(hi, f);
    used++;
  }
  if (used == 0) {
    if (skipped > 0)
      std::cout << Form("         no model both has the exits to subtract and "
                        "covers [%.2f, %.2f] MeV: nothing subtracted",
                        TMath::Min(e_out, e_in), TMath::Max(e_out, e_in))
                << std::endl;
    return 0.0;
  }
  spread = used > 1 ? 0.5 * (hi - lo) : 0.0;
  return sum / Double_t(used);
}

Bool_t CrossSection::Strip(const CrossSectionChannel &ch,
                           const std::vector<TalysCurve> &talys,
                           const std::vector<TalysCurve> &talys_sub, Int_t reac,
                           Point &pt) {
  const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;
  const TString region = "region_" + ch.name;
  TH2F *h = static_cast<TH2F *>(cache_->Get(Form("scatter_r%d", reac)));
  TCutG *cut = RegionCutStore::Load(region, reac);
  if (!h || !cut) {
    std::cout << Form("   %2d    no %s", reac,
                      h ? "region cut" : "scatter in the cache")
              << std::endl;
    delete cut;
    return kFALSE;
  }
  pt.reac = reac;
  // Beam particles reached this strip under the conditions a reaction here had
  // to satisfy; flat-count fallback if the cache predates the per-strip one.
  Bool_t at_ok = kTRUE;
  const Long64_t n_at = ReadCount(*cache_, Form("n_normed_r%d", reac), at_ok);
  pt.n_denom = Double_t(at_ok && n_at > 0 ? n_at : n_beam_);

  // Strip entrance/exit energies: strip 0's energy, less each earlier strip's
  // mean loss. e_eff from the first model; the others bound shape dependence.
  Double_t e_lab = e_strip0_;
  for (Int_t s = 0; s < reac; s++)
    e_lab -= dE_[s];
  pt.e_in = e_lab * cm_frac_;
  pt.e_out = (e_lab - dE_[reac]) * cm_frac_;
  const Bool_t use_eff = Constants::cfg.CROSS_SECTION_CONFIG.EFFECTIVE_ENERGY;
  pt.e_eff = use_eff ? EffectiveEnergy(talys.empty() ? nullptr : talys[0].axn,
                                       pt.e_out, pt.e_in)
                     : 0.5 * (pt.e_in + pt.e_out);
  pt.e_eff_lo = pt.e_eff_hi = pt.e_eff;
  for (Int_t m = 1; use_eff && m < Int_t(talys.size()); m++) {
    const Double_t e = EffectiveEnergy(talys[m].axn, pt.e_out, pt.e_in);
    pt.e_eff_lo = TMath::Min(pt.e_eff_lo, e);
    pt.e_eff_hi = TMath::Max(pt.e_eff_hi, e);
  }
  if (use_eff && talys.size() > 1)
    std::cout << Form("         E_cm,eff %.3f (midpoint %.3f); across the %zu "
                      "model shapes %.3f..%.3f -> energy systematic "
                      "+%.3f/-%.3f",
                      pt.e_eff, 0.5 * (pt.e_in + pt.e_out), talys.size(),
                      pt.e_eff_lo, pt.e_eff_hi, pt.e_eff_hi - pt.e_eff,
                      pt.e_eff - pt.e_eff_lo)
              << std::endl;

  const Double_t norm = pt.n_denom * areal_ * kCm2PerMb;
  const Bool_t all_tagged =
      C.AN_REGION_MODE == StripSumScatterConfig::AN_REGION_ALL_TAGGED;
  // Geometric = in-region / enclosed fraction; attributed = 3-sigma fit share.
  // Disagreement = region systematic; geometric runs away over beam tail.
  const Double_t n_raw = CountInCut(h, cut, 1.0);
  // Remove the TALYS-predicted branching fraction from the tagged count.
  pt.sub_frac =
      SubtractedFraction(talys, talys_sub, pt.e_out, pt.e_in, pt.sub_frac_err);
  const Double_t n_kept = n_raw * (1.0 - pt.sub_frac);
  if (pt.sub_frac > 0.0)
    std::cout << Form("         %s: TALYS puts %.2f%% +- %.2f%% of what is "
                      "tagged here in it, so %.0f of %.0f events are removed",
                      SubtractedLabel(ch).Data(), 100.0 * pt.sub_frac,
                      100.0 * pt.sub_frac_err, n_raw - n_kept, n_raw)
              << std::endl;
  if (all_tagged) {
    // Every tagged event is the reaction: no enclosed-fraction correction,
    // and the systematic is the cut variation plus the gas pressure.
    pt.n_reac = n_kept;
    pt.sigma = pt.n_reac / norm;
    // Scale the count interval by the retained fraction.
    Double_t below = 0.0, above = 0.0;
    PoissonInterval(n_raw, below, above);
    pt.stat_lo = below * (1.0 - pt.sub_frac) / norm;
    pt.stat_hi = above * (1.0 - pt.sub_frac) / norm;
    TString detail;
    const Double_t sys_cuts =
        CutVariation(reac, (1.0 - pt.sub_frac) / (areal_ * kCm2PerMb), detail);
    const Double_t sys_gas = GasPressureSys(pt.sigma);
    // Propagate the model spread on the removed fraction.
    const Double_t sys_sub = pt.sub_frac_err * n_raw / norm;
    pt.sys =
        std::sqrt(sys_cuts * sys_cuts + sys_gas * sys_gas + sys_sub * sys_sub);
    if (pt.sub_frac > 0.0)
      std::cout << Form("         subtraction systematic %.2f mb (%.1f%% of "
                        "the point), against %.2f mb cut variation and %.2f "
                        "mb gas pressure",
                        sys_sub,
                        pt.sigma > 0.0 ? 100.0 * sys_sub / pt.sigma : 0.0,
                        sys_cuts, sys_gas)
                << std::endl;
    std::cout << Form("   %2d    [%5.2f, %5.2f]      %6.2f   %6.0f  %9.0f   "
                      "%7.1f +%.1f -%.1f (+%.1f -%.1f stat%s, %.1f cut "
                      "variation, %.2f gas pressure%s; all tagged events "
                      "counted%s%s)",
                      reac, pt.e_in, pt.e_out, pt.e_eff, pt.n_reac, pt.n_denom,
                      pt.sigma, pt.ErrHi(), pt.ErrLo(), pt.stat_hi, pt.stat_lo,
                      n_raw < Constants::cfg.CROSS_SECTION_CONFIG
                                  .FELDMAN_COUSINS_MAX_COUNT
                          ? " Feldman-Cousins"
                          : "",
                      sys_cuts, sys_gas,
                      pt.sub_frac > 0.0 ? Form(", %.2f subtraction", sys_sub)
                                        : "",
                      detail.IsNull() ? "; no variation counts in the cache"
                                      : "; per condition, change of sigma "
                                        "with its threshold up (+) and "
                                        "down (-):",
                      detail.Data())
              << std::endl;
    delete cut;
    return kTRUE;
  }
  // Apply the same correction in the fitted-region mode.
  const Double_t keep = 1.0 - pt.sub_frac;
  const Double_t est_geo = n_raw * keep / Enclosed(C.AN_REGION_NSIGMA);
  const Double_t n_assigned = RegionCutStore::LoadAssigned(region, reac);
  const Bool_t have_fit = n_assigned >= 0.0;
  const Double_t est_fit =
      have_fit ? n_assigned * keep / Enclosed(3.0) : est_geo;
  pt.n_reac = est_fit;
  pt.sigma = pt.n_reac / norm;
  {
    // The Poisson count is what was counted; the point scales it by the
    // enclosed fraction, so the interval scales with it.
    const Double_t counted = have_fit ? n_assigned : n_raw;
    const Double_t scale = counted > 0.0 ? pt.sigma / counted : 0.0;
    Double_t below = 0.0, above = 0.0;
    PoissonInterval(counted, below, above);
    pt.stat_lo = below * scale;
    pt.stat_hi = above * scale;
  }
  if (have_fit) {
    // Core check: 1-sigma ellipse (half the drawn region), corrected for the
    // Apply the retained fraction to the core count.
    const Double_t est_core =
        keep * CountInCut(h, cut, 1.0 / C.AN_REGION_NSIGMA) / Enclosed(1.0);
    pt.sys = std::fabs(est_core - est_fit) / norm;
  } else {
    // A drawn region without an attributed count: scale it by 0.75x and
    // Correct each variation for enclosure and subtraction.
    Double_t lo = pt.sigma, hi = pt.sigma;
    const Double_t kScale[2] = {0.75, 1.25};
    for (Int_t k = 0; k < 2; k++) {
      const Double_t s = keep * CountInCut(h, cut, kScale[k]) /
                         Enclosed(kScale[k] * C.AN_REGION_NSIGMA) / norm;
      lo = TMath::Min(lo, s);
      hi = TMath::Max(hi, s);
    }
    pt.sys = 0.5 * (hi - lo);
  }
  {
    const Double_t sys_sub =
        keep > 0.0 ? pt.sub_frac_err / keep * pt.sigma : 0.0;
    pt.sys = std::sqrt(pt.sys * pt.sys +
                       GasPressureSys(pt.sigma) * GasPressureSys(pt.sigma) +
                       sys_sub * sys_sub);
  }
  std::cout << Form("   %2d    [%5.2f, %5.2f]      %6.2f   %6.0f  %9.0f   "
                    "%7.1f +%.1f -%.1f (+%.1f -%.1f stat, %.1f region; "
                    "in-region %.0f, attributed %s)",
                    reac, pt.e_in, pt.e_out, pt.e_eff, pt.n_reac, pt.n_denom,
                    pt.sigma, pt.ErrHi(), pt.ErrLo(), pt.stat_hi, pt.stat_lo,
                    pt.sys, est_geo, have_fit ? Form("%.0f", est_fit) : "none")
            << std::endl;
  delete cut;
  return kTRUE;
}

/// Against the published values. A published table is a list of strips, so
/// the two are aligned strip to row by the single integer offset that best
/// matches the energies over the whole table, not row by row: pairing each
/// strip with its nearest published energy silently pairs two strips with one
/// row when a dataset reaches energies the reference never reported.
void CrossSection::CompareReference(const ChannelResult &r) const {
  const std::vector<std::vector<Double_t>> &ref = r.ch->reference_xs;
  if (ref.empty())
    return;
  const Int_t n_ref = Int_t(ref.size());
  Int_t best_off = 0;
  Double_t best_rms = 1.0e9;
  for (Int_t off = 0; off + n_ref - 1 <= 16; off++) {
    Double_t s = 0.0;
    for (Int_t k = 0; k < n_ref; k++) {
      const Double_t d = e_mid_[off + k] - ref[k][0];
      s += d * d;
    }
    const Double_t rms = std::sqrt(s / Double_t(n_ref));
    if (rms < best_rms) {
      best_rms = rms;
      best_off = off;
    }
  }
  std::cout << std::endl;
  std::cout << Form("  %s vs %s: first row is strip %d from the energies, rms "
                    "%.3f MeV",
                    r.ch->name.Data(), r.ch->reference_label.Data(), best_off,
                    best_rms)
            << std::endl;
  std::cout << "  strip   E_cm,eff   published E   this work [mb]      "
               "published [mb]   ratio"
            << std::endl;
  for (Int_t i = 0; i < Int_t(r.points.size()); i++) {
    const Point &pt = r.points[i];
    const Int_t row = pt.reac - best_off;
    if (row < 0 || row >= n_ref) {
      std::cout << Form("   %2d     %6.2f      --          %7.1f          "
                        "   -- (upstream of the published range)",
                        pt.reac, pt.e_eff, pt.sigma)
                << std::endl;
      continue;
    }
    const Double_t v = ref[row][3];
    std::cout << Form("   %2d     %6.2f    %6.2f        %7.1f +%-5.1f -%-5.1f  "
                      "%7.1f (%.1f)   %5.2f",
                      pt.reac, pt.e_eff, ref[row][0], pt.sigma, pt.ErrHi(),
                      pt.ErrLo(), v, ref[row][4], v > 0.0 ? pt.sigma / v : 0.0)
              << std::endl;
  }
}

/// Excitation function(s), this work against the published points. The frame
/// spans every set drawn so the measured strips are seen in the context of the
/// whole published curve, not just the rows they happen to sit beside. One
/// channel: its own title. Several: the dataset alone, every channel labelled.
void CrossSection::Draw(const std::vector<const ChannelResult *> &rs,
                        const TString &name) const {
  Double_t fx_lo = 1.0e9, fx_hi = -1.0e9, fy_lo = 1.0e9, fy_hi = -1.0e9;
  struct Series {
    TGraphAsymmErrors *g;
    TString label;
  };
  std::vector<Series> measured, published;
  std::vector<std::pair<TGraph *, TString>> curves;
  for (Int_t c = 0; c < Int_t(rs.size()); c++) {
    const ChannelResult &r = *rs[c];
    std::vector<Double_t> vx, vy, vexl, vexh, veyl, veyh;
    for (Int_t i = 0; i < Int_t(r.points.size()); i++) {
      const Point &pt = r.points[i];
      // Energy error = strip extent about e_eff, asymmetric (above the midpoint
      // on a rising curve); e_eff's shape dependence added in quadrature.
      vx.push_back(pt.e_eff);
      vy.push_back(pt.sigma);
      vexl.push_back(std::hypot(pt.e_eff - TMath::Min(pt.e_in, pt.e_out),
                                pt.e_eff - pt.e_eff_lo));
      vexh.push_back(std::hypot(TMath::Max(pt.e_in, pt.e_out) - pt.e_eff,
                                pt.e_eff_hi - pt.e_eff));
      veyl.push_back(pt.ErrLo());
      veyh.push_back(pt.ErrHi());
      fx_lo = TMath::Min(fx_lo, vx.back() - vexl.back());
      fx_hi = TMath::Max(fx_hi, vx.back() + vexh.back());
      fy_lo = TMath::Min(fy_lo, vy.back());
      fy_hi = TMath::Max(fy_hi, vy.back());
    }
    if (vx.empty())
      continue;
    Series s;
    s.g = new TGraphAsymmErrors(Int_t(vx.size()), &vx[0], &vy[0], &vexl[0],
                                &vexh[0], &veyl[0], &veyh[0]);
    s.g->SetMarkerStyle(kChannelMarker[c % 4]);
    s.g->SetMarkerColor(kChannelColor[c % 4]);
    s.g->SetLineColor(kChannelColor[c % 4]);
    s.label = rs.size() > 1 ? "Present Work " + Label(*r.ch) : "Present Work";
    measured.push_back(s);
    const std::vector<std::vector<Double_t>> &ref = r.ch->reference_xs;
    if (!ref.empty()) {
      std::vector<Double_t> rx, ry, rexl, rexh, rey;
      for (Int_t k = 0; k < Int_t(ref.size()); k++) {
        rx.push_back(ref[k][0]);
        rexh.push_back(ref[k][1]);
        rexl.push_back(ref[k][2]);
        ry.push_back(ref[k][3]);
        rey.push_back(ref[k][4]);
        fx_lo = TMath::Min(fx_lo, ref[k][0]);
        fx_hi = TMath::Max(fx_hi, ref[k][0]);
        fy_lo = TMath::Min(fy_lo, ref[k][3]);
        fy_hi = TMath::Max(fy_hi, ref[k][3]);
      }
      Series p;
      p.g = new TGraphAsymmErrors(Int_t(rx.size()), &rx[0], &ry[0], &rexl[0],
                                  &rexh[0], &rey[0], &rey[0]);
      p.g->SetMarkerStyle(24 + (c % 4));
      p.g->SetMarkerColor(kRed + 1);
      p.g->SetLineColor(kRed + 1);
      p.label = rs.size() > 1 ? r.ch->reference_label + " " + Label(*r.ch)
                              : r.ch->reference_label;
      published.push_back(p);
    }
    // Hauser-Feshbach prediction, if talys-xs wrote one, one colour per
    // model; combined figure: first per channel only, in its colour.
    for (Int_t m = 0; m < Int_t(r.talys.size()); m++) {
      if (rs.size() > 1 && m > 0)
        break;
      curves.push_back(std::make_pair(
          r.talys[m].axn, rs.size() > 1 ? r.talys[m].label + " " + Label(*r.ch)
                                        : r.talys[m].label));
    }
  }
  if (measured.empty())
    return;
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kTRUE);
  TH1F *frame =
      c->DrawFrame(fx_lo - 0.4, 0.5 * fy_lo, fx_hi + 0.4, 3.0 * fy_hi);
  frame->SetTitle(Form("%s%s;%s [MeV];#sigma [mb]", Paths::DatasetName().Data(),
                       rs.size() == 1 ? Label(*rs[0]->ch).Data() : "",
                       Constants::cfg.CROSS_SECTION_CONFIG.EFFECTIVE_ENERGY
                           ? "E_{c.m.,eff}"
                           : "E_{c.m.}"));
  // Curves first so the points sit on top, and only over the frame's
  // energies.
  std::vector<std::pair<TGraph *, TString>> drawn;
  for (Int_t k = 0; k < Int_t(curves.size()); k++) {
    TGraph *g = Clipped(curves[k].first, fx_lo - 0.4, fx_hi + 0.4);
    if (!g)
      continue;
    // Combined figure: the channel's colour. One channel: a colour per
    // model. All dashed: the colour is the distinction, the dash says model.
    const Int_t ch = rs.size() > 1 ? k : 0;
    g->SetLineColor(rs.size() > 1 ? kChannelColor[ch % 4] : kModelColor[k % 8]);
    g->SetLineWidth(2);
    g->SetLineStyle(2);
    g->Draw("L SAME");
    drawn.push_back(std::make_pair(g, curves[k].second));
  }
  for (Int_t k = 0; k < Int_t(published.size()); k++)
    published[k].g->Draw("P SAME");
  for (Int_t k = 0; k < Int_t(measured.size()); k++)
    measured[k].g->Draw("P SAME");
  // Bottom right is the one empty corner: the excitation function climbs to
  // the upper right and the reference table starts at the lower left.
  const Int_t n_entries =
      Int_t(measured.size() + published.size() + drawn.size());
  TLegend *leg =
      PlottingUtils::AddLegend(0.42, 0.89, 0.16, 0.18 + 0.048 * n_entries);
  for (Int_t k = 0; k < Int_t(measured.size()); k++)
    leg->AddEntry(measured[k].g, measured[k].label, "pe");
  for (Int_t k = 0; k < Int_t(published.size()); k++)
    leg->AddEntry(published[k].g, published[k].label, "pe");
  for (Int_t k = 0; k < Int_t(drawn.size()); k++)
    leg->AddEntry(drawn[k].first, drawn[k].second, "l");
  leg->Draw();
  // The stamp, drawn last so it sits over the frame.
  StampPreliminary(c);
  PlottingUtils::SaveFigure(c, name, "cross_section", PlotSaveOptions::kLOG);
  delete c;
}

void CrossSection::StampPreliminary(const TPad *pad) {
  // The stamp: red, bold, in the top-left corner of the plot area (just
  // inside the frame's margins).
  if (!Constants::cfg.CROSS_SECTION_CONFIG.PRELIMINARY)
    return;
  TLatex *stamp = new TLatex();
  stamp->SetNDC();
  stamp->SetTextAlign(13); // left, top
  stamp->SetTextSize(0.05);
  stamp->SetTextFont(62);
  stamp->SetTextColor(kRed + 1);
  stamp->DrawLatex(pad->GetLeftMargin() + 0.02,
                   1.0 - pad->GetTopMargin() - 0.02, "PRELIMINARY");
}

Bool_t CrossSection::CollectFitComponents(const ChannelResult &r, Int_t m,
                                          FitComponents &comp) const {
  const CrossSectionConfig &X = Constants::cfg.CROSS_SECTION_CONFIG;
  const CrossSectionChannel &ch = *r.ch;
  for (Int_t k = 0; k < Int_t(ch.talys_exits.size()); k++) {
    Int_t z = 0, a = 0;
    ExitResidue(ch.talys_exits[k], X.BEAM_Z, X.BEAM_A, z, a);
    std::map<std::pair<Int_t, Int_t>, TGraph *>::const_iterator it =
        talys_raw_[m].find(std::make_pair(z, a));
    if (it == talys_raw_[m].end())
      continue;
    comp.graph.push_back(it->second);
    comp.label.push_back("(#alpha, " + ch.talys_exits[k] + ")");
    comp.tag.push_back(ch.talys_exits[k]);
  }
  comp.n_keep = comp.Size();
  // Include subtracted exits as free fit components.
  for (Int_t k = 0; k < Int_t(ch.subtract_exits.size()); k++) {
    Int_t z = 0, a = 0;
    ExitResidue(ch.subtract_exits[k], X.BEAM_Z, X.BEAM_A, z, a);
    std::map<std::pair<Int_t, Int_t>, TGraph *>::const_iterator it =
        talys_raw_[m].find(std::make_pair(z, a));
    if (it == talys_raw_[m].end())
      continue;
    comp.graph.push_back(it->second);
    comp.label.push_back("(#alpha, " + ch.subtract_exits[k] + ") subtracted");
    comp.tag.push_back(ch.subtract_exits[k]);
  }
  if (comp.n_keep == 0) {
    std::cerr << "cross-section: channel " << ch.name << ": "
              << talys_labels_[m] << " has none of the channel's exits"
              << std::endl;
    return kFALSE;
  }
  return kTRUE;
}

void CrossSection::FitExitScales(
    const ChannelResult &r, const FitComponents &comp, ScaleFit &fs,
    const std::vector<Double_t> *point_scale) const {
  const Int_t nc = comp.Size();
  const Int_t np = Int_t(r.points.size());

  // Weighted linear least squares for the scales, dropping exits that are
  // zero at every point (their scale is undetermined).
  for (Int_t k = 0; k < nc; k++) {
    Bool_t any = kFALSE;
    for (Int_t i = 0; i < np && !any; i++)
      any = Interpolated(comp.graph[k], r.points[i].e_eff) > 0.0;
    if (any)
      fs.free_idx.push_back(k);
  }
  fs.scale.assign(nc, 1.0);
  fs.scale_err.assign(nc, 0.0);
  fs.scale_raw.assign(nc, 1.0);
  fs.scale_raw_err.assign(nc, 0.0);
  Int_t nf = Int_t(fs.free_idx.size());
  fs.cov.ResizeTo(nf, nf);
  Bool_t first_pass = kTRUE;
  while (nf > 0 && np >= nf) {
    TMatrixD mat(nf, nf);
    TVectorD rhs(nf);
    for (Int_t i = 0; i < np; i++) {
      const Point &pt = r.points[i];
      const Double_t k = point_scale ? (*point_scale)[i] : 1.0;
      const Double_t err = k * 0.5 * (pt.ErrLo() + pt.ErrHi());
      if (!(err > 0.0))
        continue;
      const Double_t w = 1.0 / (err * err);
      for (Int_t a = 0; a < nf; a++) {
        const Double_t ca = Interpolated(comp.graph[fs.free_idx[a]], pt.e_eff);
        rhs[a] += w * ca * k * pt.sigma;
        for (Int_t b = 0; b < nf; b++)
          mat[a][b] +=
              w * ca * Interpolated(comp.graph[fs.free_idx[b]], pt.e_eff);
      }
    }
    fs.cov.ResizeTo(nf, nf);
    fs.cov = mat;
    fs.cov.Invert();
    const TVectorD sol = fs.cov * rhs;
    // Preserve the raw solution before nonnegative pinning.
    if (first_pass) {
      for (Int_t a = 0; a < nf; a++) {
        fs.scale_raw[fs.free_idx[a]] = sol[a];
        fs.scale_raw_err[fs.free_idx[a]] =
            std::sqrt(TMath::Max(0.0, fs.cov[a][a]));
      }
      first_pass = kFALSE;
    }
    Int_t worst = -1;
    for (Int_t a = 0; a < nf; a++)
      if (sol[a] < 0.0 && (worst < 0 || sol[a] < sol[worst]))
        worst = a;
    if (worst >= 0) {
      fs.pinned.push_back(fs.free_idx[worst]);
      fs.scale[fs.free_idx[worst]] = 0.0;
      fs.free_idx.erase(fs.free_idx.begin() + worst);
      nf--;
      continue;
    }
    for (Int_t a = 0; a < nf; a++) {
      fs.scale[fs.free_idx[a]] = sol[a];
      fs.scale_err[fs.free_idx[a]] = std::sqrt(TMath::Max(0.0, fs.cov[a][a]));
    }
    break;
  }
  fs.nf = nf;
  for (Int_t i = 0; i < np; i++) {
    const Point &pt = r.points[i];
    const Double_t k = point_scale ? (*point_scale)[i] : 1.0;
    const Double_t err = k * 0.5 * (pt.ErrLo() + pt.ErrHi());
    Double_t f = 0.0;
    for (Int_t c = 0; c < nc; c++)
      f += fs.scale[c] * Interpolated(comp.graph[c], pt.e_eff);
    if (err > 0.0)
      fs.chi2 += std::pow((k * pt.sigma - f) / err, 2);
  }
}

void CrossSection::PrintScaleFit(const CrossSectionChannel &ch, Int_t m,
                                 const FitComponents &comp, const ScaleFit &fs,
                                 Int_t np) const {
  const Int_t nc = comp.Size();
  std::cout << std::endl
            << "  " << ch.name << " " << Label(ch) << ": " << talys_labels_[m]
            << " scaled to this work's points" << std::endl;
  for (Int_t k = 0; k < nc; k++) {
    const Bool_t is_pinned =
        std::find(fs.pinned.begin(), fs.pinned.end(), k) != fs.pinned.end();
    std::cout << Form("    %-26s scale %.4f +- %.4f%s", comp.label[k].Data(),
                      fs.scale[k], fs.scale_err[k],
                      is_pinned ? "  (went negative: pinned at 0)"
                      : fs.scale_err[k] > 0.0
                          ? ""
                          : "  (zero at every point: fixed)")
              << std::endl;
  }
  for (Int_t a = 0; a < fs.nf; a++)
    for (Int_t b = a + 1; b < fs.nf; b++)
      std::cout << Form("    correlation %s-%s %.3f",
                        comp.label[fs.free_idx[a]].Data(),
                        comp.label[fs.free_idx[b]].Data(),
                        fs.cov[a][b] / std::sqrt(fs.cov[a][a] * fs.cov[b][b]))
                << std::endl;
  std::cout << Form("    chi2 %.2f for %d points, %d free scale(s)", fs.chi2,
                    np, fs.nf)
            << std::endl;
}

/// Compare the subtracted-exit fit before and after correction.
void CrossSection::PrintSubtractionCheck(const ChannelResult &r,
                                         const FitComponents &comp,
                                         const ScaleFit &fs) const {
  if (comp.n_keep >= comp.Size())
    return;
  const Int_t np = Int_t(r.points.size());
  // Undo the correction for the comparison fit.
  std::vector<Double_t> undo(np, 1.0);
  Bool_t have_undo = kTRUE;
  for (Int_t i = 0; i < np; i++) {
    const Double_t keep = 1.0 - r.points[i].sub_frac;
    if (!(keep > 0.0)) {
      have_undo = kFALSE;
      break;
    }
    undo[i] = 1.0 / keep;
  }
  ScaleFit before;
  if (have_undo)
    FitExitScales(r, comp, before, &undo);

  std::cout << "    subtraction check: the subtracted exits' scales with the "
               "points as measured and as corrected"
            << std::endl;
  for (Int_t k = comp.n_keep; k < comp.Size(); k++) {
    const Double_t after = fs.scale_raw[k], after_e = fs.scale_raw_err[k];
    if (!have_undo) {
      std::cout << Form("      %-26s after %.3f +- %.3f", comp.label[k].Data(),
                        after, after_e)
                << std::endl;
      continue;
    }
    // The expected change is the kept-channel scale.
    const Double_t raw = before.scale_raw[k];
    const Double_t expected = comp.n_keep > 0 ? -before.scale_raw[0] : -1.0;
    std::cout << Form("      %-26s before %.3f, after %.3f (+- %.3f): moved "
                      "%+.3f, expected %+.3f",
                      comp.label[k].Data(), raw, after, after_e, after - raw,
                      expected)
              << std::endl;
  }
  // Report whether the component fit is informative.
  Double_t rho = 0.0;
  for (Int_t a = 0; a < fs.nf; a++)
    for (Int_t b = a + 1; b < fs.nf; b++)
      if (fs.cov[a][a] > 0.0 && fs.cov[b][b] > 0.0) {
        const Double_t c =
            fs.cov[a][b] / std::sqrt(fs.cov[a][a] * fs.cov[b][b]);
        if (std::fabs(c) > std::fabs(rho))
          rho = c;
      }
  const Double_t chi2_ndf =
      np > fs.nf ? fs.chi2 / Double_t(np - fs.nf) : fs.chi2;
  const Bool_t degenerate = std::fabs(rho) > 0.9;
  const Bool_t poor_fit = chi2_ndf > 5.0;
  if (degenerate || poor_fit)
    std::cout << Form("      the scale is not a usable check here: %s%s%s. A "
                      "component this close to the channel's own absorbs the "
                      "shape mismatch rather than its own yield, so its value "
                      "says nothing about how much %s is left. The "
                      "subtraction rests on the TALYS branching, and its "
                      "uncertainty is the model spread already in the "
                      "systematic.",
                      degenerate
                          ? Form("the two scales are %.2f correlated", rho)
                          : "",
                      degenerate && poor_fit ? " and " : "",
                      poor_fit ? Form("chi2/ndf is %.1f, so the model does not "
                                      "describe the shape",
                                      chi2_ndf)
                               : "",
                      SubtractedLabel(*r.ch).Data())
              << std::endl;
  else
    std::cout << Form("      chi2/ndf %.1f and correlation %.2f, so the move "
                      "is meaningful",
                      chi2_ndf, rho)
              << std::endl;
}

Bool_t CrossSection::BuildScaledCurves(const ChannelResult &r,
                                       const FitComponents &comp,
                                       const ScaleFit &fs, TGraph *&sum,
                                       TGraph *&fit,
                                       TGraphErrors *&band) const {
  const Int_t nc = comp.Size();
  // The scaled sum and its 3-sigma band on the model's grid, and the
  // The unscaled sum contains kept exits only.
  std::vector<Double_t> gx, gsum, gfit, gband;
  for (Int_t p = 0; p < comp.graph[0]->GetN(); p++) {
    const Double_t e = comp.graph[0]->GetX()[p];
    Double_t s = 0.0, f = 0.0, var = 0.0;
    std::vector<Double_t> c(nc);
    for (Int_t k = 0; k < nc; k++) {
      c[k] = Interpolated(comp.graph[k], e);
      if (k < comp.n_keep)
        s += c[k];
      f += fs.scale[k] * c[k];
    }
    for (Int_t a = 0; a < fs.nf; a++)
      for (Int_t b = 0; b < fs.nf; b++)
        var += c[fs.free_idx[a]] * c[fs.free_idx[b]] * fs.cov[a][b];
    if (!(s > 0.0))
      continue;
    gx.push_back(e);
    gsum.push_back(s);
    gfit.push_back(f);
    gband.push_back(3.0 * std::sqrt(TMath::Max(0.0, var)));
  }
  if (gx.empty())
    return kFALSE;
  sum = new TGraph(Int_t(gx.size()), &gx[0], &gsum[0]);
  fit = new TGraph(Int_t(gx.size()), &gx[0], &gfit[0]);
  std::vector<Double_t> zero(gx.size(), 0.0);
  band =
      new TGraphErrors(Int_t(gx.size()), &gx[0], &gfit[0], &zero[0], &gband[0]);
  return kTRUE;
}

void CrossSection::BuildMeasuredGraphs(const ChannelResult &r, TGraph *fit,
                                       TGraphAsymmErrors *&measured,
                                       TGraphAsymmErrors *&deviation,
                                       FitFrame &frame) const {
  const Int_t np = Int_t(r.points.size());
  // The measured points and the deviation of each point from the scaled
  // sum, printed as a row.
  std::vector<Double_t> vx, vy, vexl, vexh, veyl, veyh, dy, deyl, deyh;
  for (Int_t i = 0; i < np; i++) {
    const Point &pt = r.points[i];
    vx.push_back(pt.e_eff);
    vy.push_back(pt.sigma);
    vexl.push_back(std::hypot(pt.e_eff - TMath::Min(pt.e_in, pt.e_out),
                              pt.e_eff - pt.e_eff_lo));
    vexh.push_back(std::hypot(TMath::Max(pt.e_in, pt.e_out) - pt.e_eff,
                              pt.e_eff_hi - pt.e_eff));
    veyl.push_back(pt.ErrLo());
    veyh.push_back(pt.ErrHi());
    frame.fx_lo = TMath::Min(frame.fx_lo, vx.back() - vexl.back());
    frame.fx_hi = TMath::Max(frame.fx_hi, vx.back() + vexh.back());
    frame.fy_lo = TMath::Min(frame.fy_lo, vy.back());
    frame.fy_hi = TMath::Max(frame.fy_hi, vy.back());
    const Double_t f = Interpolated(fit, pt.e_eff);
    const Double_t err = pt.sigma > f ? pt.ErrLo() : pt.ErrHi();
    const Double_t pull = err > 0.0 ? (pt.sigma - f) / err : 0.0;
    dy.push_back(pull);
    deyl.push_back(1.0);
    deyh.push_back(1.0);
    frame.pull_max = TMath::Max(frame.pull_max, std::abs(pull) + 1.0);
    std::cout << Form("     %2d     %6.2f      %8.2f          %8.2f          "
                      "%+6.2f",
                      pt.reac, pt.e_eff, pt.sigma, f, pull)
              << std::endl;
  }
  measured = new TGraphAsymmErrors(np, &vx[0], &vy[0], &vexl[0], &vexh[0],
                                   &veyl[0], &veyh[0]);
  measured->SetMarkerStyle(kChannelMarker[0]);
  measured->SetMarkerColor(kChannelColor[0]);
  measured->SetLineColor(kChannelColor[0]);
  deviation = new TGraphAsymmErrors(np, &vx[0], &dy[0], &vexl[0], &vexh[0],
                                    &deyl[0], &deyh[0]);
  deviation->SetMarkerStyle(kChannelMarker[0]);
  deviation->SetMarkerColor(kChannelColor[0]);
  deviation->SetLineColor(kChannelColor[0]);
}

TGraphAsymmErrors *
CrossSection::BuildPublishedGraph(const CrossSectionChannel &ch,
                                  FitFrame &frame) const {
  TGraphAsymmErrors *published = nullptr;
  const std::vector<std::vector<Double_t>> &ref = ch.reference_xs;
  if (!ref.empty()) {
    std::vector<Double_t> rx, ry, rexl, rexh, rey;
    for (Int_t k = 0; k < Int_t(ref.size()); k++) {
      rx.push_back(ref[k][0]);
      rexh.push_back(ref[k][1]);
      rexl.push_back(ref[k][2]);
      ry.push_back(ref[k][3]);
      rey.push_back(ref[k][4]);
      frame.fx_lo = TMath::Min(frame.fx_lo, ref[k][0]);
      frame.fx_hi = TMath::Max(frame.fx_hi, ref[k][0]);
      frame.fy_lo = TMath::Min(frame.fy_lo, ref[k][3]);
      frame.fy_hi = TMath::Max(frame.fy_hi, ref[k][3]);
    }
    published = new TGraphAsymmErrors(Int_t(rx.size()), &rx[0], &ry[0],
                                      &rexl[0], &rexh[0], &rey[0], &rey[0]);
    published->SetMarkerStyle(24);
    published->SetMarkerColor(kRed + 1);
    published->SetLineColor(kRed + 1);
  }
  return published;
}

void CrossSection::FitLegendWidthPx(const ScaleFit &fs,
                                    const CrossSectionChannel &ch,
                                    const FitComponents &comp,
                                    const TString &model_label,
                                    Bool_t has_published, TString &scale_text,
                                    Double_t &legend_width_px) const {
  scale_text = "#times";
  for (Int_t k = 0; k < comp.Size(); k++)
    scale_text +=
        Form("%s %.2f_{%s}", k ? " &" : "", fs.scale[k], comp.tag[k].Data());
  std::vector<TString> legend_labels;
  legend_labels.push_back("Present Work");
  if (has_published)
    legend_labels.push_back(ch.reference_label);
  legend_labels.push_back(Label(ch) + " scaled");
  legend_labels.push_back(scale_text);
  legend_labels.push_back(Label(ch));
  if (comp.Size() > 1)
    for (Int_t k = 0; k < comp.Size(); k++)
      legend_labels.push_back(comp.label[k]);
  legend_width_px = PlottingUtils::LegendWidthPx(legend_labels, model_label);
}

void CrossSection::DrawFitTopPanel(
    const CrossSectionChannel &ch, const FitFrame &frame,
    const TString &scale_text, Double_t legend_width_px, TGraph *sum,
    TGraph *fit, TGraphErrors *band, const FitComponents &comp,
    TGraphAsymmErrors *published, TGraphAsymmErrors *measured,
    FitPanels &pads) const {
  const Double_t x_lo = frame.fx_lo - 0.4, x_hi = frame.fx_hi + 0.4;
  const Int_t nc = comp.Size();
  const Int_t panel_height_px = 205;
  const Int_t plot_right_margin_px =
      TMath::Nint(1200.0 * gStyle->GetPadRightMargin());
  const Double_t legend_gap_px = 40.0;
  const Int_t extra_width_px =
      TMath::Max(0, TMath::Nint(legend_width_px + 2.0 * legend_gap_px) -
                        plot_right_margin_px);
  TCanvas *c = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      pads.plot, pads.side, extra_width_px, 800 + panel_height_px, kTRUE);
  TPad *top = nullptr;
  TPad *bottom = nullptr;
  PlottingUtils::SplitPadForPanel(pads.plot, panel_height_px, top, bottom);
  pads.c = c;
  pads.top = top;
  pads.bottom = bottom;
  pads.plot_right_margin_px = plot_right_margin_px;
  pads.extra_width_px = extra_width_px;

  pads.top->cd();
  TH1F *frame_hist =
      pads.top->DrawFrame(x_lo, 0.5 * frame.fy_lo, x_hi, 3.0 * frame.fy_hi);
  frame_hist->SetTitle(";;#sigma [mb]");
  frame_hist->GetXaxis()->SetLabelSize(0);
  frame_hist->GetXaxis()->SetTitleSize(0);

  {
    std::vector<Double_t> bx, by, bz, be;
    for (Int_t p = 0; p < band->GetN(); p++)
      if (band->GetX()[p] >= x_lo && band->GetX()[p] <= x_hi) {
        bx.push_back(band->GetX()[p]);
        by.push_back(band->GetY()[p]);
        bz.push_back(0.0);
        be.push_back(band->GetEY()[p]);
      }
    if (!bx.empty()) {
      TGraphErrors *band_c =
          new TGraphErrors(Int_t(bx.size()), &bx[0], &by[0], &bz[0], &be[0]);
      band_c->SetFillColorAlpha(kRed + 1, 0.25);
      band_c->SetLineWidth(0);
      band_c->Draw("3 SAME");
    }
  }
  if (TGraph *g = Clipped(fit, x_lo, x_hi)) {
    PlottingUtils::ConfigureGraph(g, kRed + 1);
    g->Draw("L SAME");
    pads.entries.push_back(std::make_pair(g, Label(ch) + " scaled"));
    pads.entries.push_back(std::make_pair((TObject *)nullptr, scale_text));
  }
  if (TGraph *g = Clipped(sum, x_lo, x_hi)) {
    PlottingUtils::ConfigureGraph(g, kBlack);
    g->Draw("L SAME");
    pads.entries.push_back(std::make_pair(g, Label(ch)));
  }
  if (nc > 1)
    for (Int_t k = 0; k < nc; k++) {
      TGraph *g = Clipped(comp.graph[k], x_lo, x_hi);
      if (!g)
        continue;
      PlottingUtils::ConfigureGraph(g, kComponentColor[k % 4]);
      // Subtracted exits are dotted and grey.
      const Bool_t subtracted = k >= comp.n_keep;
      if (subtracted)
        g->SetLineColor(kGray + 1);
      g->SetLineStyle(subtracted ? 3 : 2);
      g->Draw("L SAME");
      pads.entries.push_back(std::make_pair(g, comp.label[k]));
    }
  if (published)
    published->Draw("P SAME");
  measured->Draw("P SAME");
  StampPreliminary(pads.top);
}

void CrossSection::DrawPullPanel(const FitFrame &frame,
                                 TGraphAsymmErrors *deviation,
                                 FitPanels &pads) const {
  const CrossSectionConfig &X = Constants::cfg.CROSS_SECTION_CONFIG;
  const Double_t x_lo = frame.fx_lo - 0.4, x_hi = frame.fx_hi + 0.4;
  pads.bottom->cd();
  const Double_t d_lim = TMath::Max(5.0, 2.0 * std::ceil(0.6 * frame.pull_max));
  TH1F *dframe = pads.bottom->DrawFrame(x_lo, -d_lim, x_hi, d_lim);
  dframe->SetTitle(Form(";%s [MeV];#delta/#sigma",
                        X.EFFECTIVE_ENERGY ? "E_{c.m.,eff}" : "E_{c.m.}"));
  dframe->GetYaxis()->SetNdivisions(-202);
  dframe->GetYaxis()->CenterTitle(kTRUE);
  TLine *zero_line = new TLine(x_lo, 0.0, x_hi, 0.0);
  zero_line->SetLineStyle(2);
  zero_line->SetLineColor(kBlack);
  zero_line->SetLineWidth(PlottingUtils::GetLineWidth());
  zero_line->Draw("SAME");
  TLine *plus3_line = new TLine(x_lo, 3.0, x_hi, 3.0);
  plus3_line->SetLineStyle(3);
  plus3_line->SetLineColor(kGray + 2);
  plus3_line->SetLineWidth(PlottingUtils::GetLineWidth());
  plus3_line->Draw("SAME");
  TLine *minus3_line = new TLine(x_lo, -3.0, x_hi, -3.0);
  minus3_line->SetLineStyle(3);
  minus3_line->SetLineColor(kGray + 2);
  minus3_line->SetLineWidth(PlottingUtils::GetLineWidth());
  minus3_line->Draw("SAME");
  deviation->Draw("P SAME");
}

/// The fit figure: the chosen model's per-exit curves and their sum, unscaled,
/// then the sum with one free scale per exit fitted to this work's points by
/// weighted linear least squares (each point weighted by its total error,
/// the asymmetric sides averaged), with the 3-sigma band from the fit's
/// covariance and a panel of the points' deviation from the scaled sum. An
/// exit whose curve is zero at every point cannot be scaled and keeps 1.
void CrossSection::DrawFit(const ChannelResult &r, const TString &name) const {
  const CrossSectionChannel &ch = *r.ch;
  const Int_t m = ch.fit_model;
  if (m < 0 || m >= Int_t(talys_raw_.size())) {
    std::cerr << "cross-section: channel " << ch.name << ": fit_model " << m
              << " but " << talys_raw_.size() << " TALYS model(s) in the file"
              << std::endl;
    return;
  }
  FitComponents comp;
  if (!CollectFitComponents(r, m, comp))
    return;
  const Int_t np = Int_t(r.points.size());

  ScaleFit fs;
  FitExitScales(r, comp, fs);
  PrintScaleFit(ch, m, comp, fs, np);
  PrintSubtractionCheck(r, comp, fs);
  std::cout << "    strip   E_cm,eff   this work [mb]   scaled model [mb]   "
               "pull"
            << std::endl;

  TGraph *sum = nullptr;
  TGraph *fit = nullptr;
  TGraphErrors *band = nullptr;
  if (!BuildScaledCurves(r, comp, fs, sum, fit, band))
    return;

  FitFrame frame;
  TGraphAsymmErrors *measured = nullptr;
  TGraphAsymmErrors *deviation = nullptr;
  BuildMeasuredGraphs(r, fit, measured, deviation, frame);

  // The published points extend the frame before x_lo/x_hi are formed.
  TGraphAsymmErrors *published = BuildPublishedGraph(ch, frame);

  TString scale_text;
  Double_t legend_width_px = 0.0;
  FitLegendWidthPx(fs, ch, comp, talys_labels_[m], published != nullptr,
                   scale_text, legend_width_px);

  FitPanels pads;
  DrawFitTopPanel(ch, frame, scale_text, legend_width_px, sum, fit, band, comp,
                  published, measured, pads);
  DrawPullPanel(frame, deviation, pads);

  pads.side->cd();
  const Int_t n_entries = 2 + (published ? 1 : 0) + Int_t(pads.entries.size());
  const Double_t side_width_px =
      pads.plot_right_margin_px + pads.extra_width_px;
  const Double_t side_gap_px = 0.5 * (side_width_px - legend_width_px);
  TLegend *leg = PlottingUtils::AddLegend(
      side_gap_px / side_width_px,
      (side_gap_px + legend_width_px) / side_width_px, 0.5 - 0.035 * n_entries,
      0.5 + 0.035 * n_entries);
  leg->SetHeader(talys_labels_[m]);
  leg->AddEntry(measured, "Present Work", "pe");
  if (published)
    leg->AddEntry(published, ch.reference_label, "pe");
  for (Int_t k = 0; k < Int_t(pads.entries.size()); k++)
    leg->AddEntry(pads.entries[k].first, pads.entries[k].second,
                  pads.entries[k].first ? "l" : "");

  PlottingUtils::ScaleFigure(pads.plot);
  PlottingUtils::DrawTitle(
      pads.top, Form("%s%s", Paths::DatasetName().Data(), Label(ch).Data()),
      PlottingUtils::FigureScale(pads.plot));
  PlottingUtils::SaveFigure(pads.c, name, "cross_section",
                            PlotSaveOptions::kLOG);
  delete pads.c;
}

Bool_t CrossSection::RunChannel(const CrossSectionChannel &ch,
                                ChannelResult &out) {
  const CrossSectionConfig &X = Constants::cfg.CROSS_SECTION_CONFIG;
  out.ch = &ch;
  for (Int_t k = 0; k < Int_t(ch.talys_exits.size()); k++) {
    Int_t z = 0, a = 0;
    if (!ExitResidue(ch.talys_exits[k], X.BEAM_Z, X.BEAM_A, z, a)) {
      std::cerr << "cross-section: channel " << ch.name << ": exit \""
                << ch.talys_exits[k]
                << "\" does not parse (light particles n p d t h a g with an "
                   "optional multiplicity digit)"
                << std::endl;
      return kFALSE;
    }
  }
  for (Int_t k = 0; k < Int_t(ch.subtract_exits.size()); k++) {
    Int_t z = 0, a = 0;
    if (!ExitResidue(ch.subtract_exits[k], X.BEAM_Z, X.BEAM_A, z, a)) {
      std::cerr << "cross-section: channel " << ch.name
                << ": subtracted exit \"" << ch.subtract_exits[k]
                << "\" does not parse (light particles n p d t h a g with an "
                   "optional multiplicity digit)"
                << std::endl;
      return kFALSE;
    }
  }
  if (!talys_raw_.empty() && ch.talys_exits.empty()) {
    std::cerr << "cross-section: channel " << ch.name
              << ": TALYS models are declared but the channel names no exits"
              << std::endl;
    return kFALSE;
  }
  // Subtraction requires the TALYS curves.
  if (!ch.subtract_exits.empty() && talys_raw_.empty()) {
    std::cerr << "cross-section: channel " << ch.name
              << ": subtract_exits needs the TALYS curves it subtracts; run "
                 "talys-xs first"
              << std::endl;
    return kFALSE;
  }
  out.talys = ExitCurves(ch, ch.talys_exits, "the channel:");
  out.talys_sub = ExitCurves(ch, ch.subtract_exits, "to subtract:");
  if (!ch.subtract_exits.empty() && out.talys_sub.empty()) {
    std::cerr << "cross-section: channel " << ch.name
              << ": no model has any of the exits to subtract" << std::endl;
    return kFALSE;
  }
  std::cout << std::endl
            << "  " << ch.name << " " << Label(ch) << ", region_" << ch.name
            << (out.talys.empty() ? " (no TALYS curve: midpoint energies)" : "")
            << std::endl;
  if (!out.talys_sub.empty())
    std::cout << "  " << SubtractedLabel(ch)
              << " is not separable from it here, so TALYS's predicted share "
                 "of each strip's count is removed; the models' spread on "
                 "that share is a systematic"
              << std::endl;
  std::cout << "  strip   E_cm range [MeV]    "
            << (X.EFFECTIVE_ENERGY ? "E_cm,eff" : "E_cm,mid")
            << "   N_reac    N_beam    sigma [mb]" << std::endl;
  for (Int_t reac = X.XS_STRIP_MIN; reac <= X.XS_STRIP_MAX; reac++) {
    Point pt;
    if (Strip(ch, out.talys, out.talys_sub, reac, pt))
      out.points.push_back(pt);
  }
  if (out.points.empty()) {
    std::cerr << "cross-section: channel " << ch.name
              << ": no strip produced a cross section" << std::endl;
    return kFALSE;
  }
  CompareReference(out);
  std::vector<const ChannelResult *> one(1, &out);
  Draw(one, "cross_section_" + ch.name + TagSuffix());
  if (ch.fit_model >= 0)
    DrawFit(out, "cross_section_" + ch.name + "_fit" + TagSuffix());
  return kTRUE;
}

Bool_t CrossSection::Run() {
  const CrossSectionConfig &X = Constants::cfg.CROSS_SECTION_CONFIG;
  const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;
  if (!(X.GAS_PRESSURE_TORR > 0.0) || X.BEAM_A <= 0) {
    std::cerr << "cross-section: this dataset has no CROSS_SECTION_CONFIG "
                 "(gas pressure and beam are both required)"
              << std::endl;
    return kFALSE;
  }
  if (X.CHANNELS.empty()) {
    std::cerr << "cross-section: this dataset declares no "
                 "CROSS_SECTION_CONFIG.CHANNELS"
              << std::endl;
    return kFALSE;
  }
  if (!LoadCache() || !LoadBeam())
    return kFALSE;

  std::cout << "cross-section: " << n_beam_
            << " beam events past every pre-tag "
            << "cut, of " << n_seen_ << " seen" << std::endl;
  std::cout << Form("  gas %.0f Torr at %.0f K -> %.4e atoms/cm^3; strip "
                    "%.3f cm -> %.4e atoms/cm^2",
                    X.GAS_PRESSURE_TORR, kGasTemperatureK, n_gas_,
                    kStripLengthCm, areal_)
            << std::endl;
  if (C.AN_REGION_MODE == StripSumScatterConfig::AN_REGION_ALL_TAGGED)
    std::cout << Form("  every tagged event counted; the systematic is the cut "
                      "variation with the gas pressure (+- %.1f Torr)",
                      X.GAS_PRESSURE_TORR_ERR)
              << std::endl;
  else
    std::cout << Form("  regions at %.1f sigma enclose %.1f%% of the fitted "
                      "component and counts are corrected for that only",
                      C.AN_REGION_NSIGMA, 100.0 * Enclosed(C.AN_REGION_NSIGMA))
              << std::endl;
  LoadTalys();
  if (talys_raw_.empty())
    std::cout << "cross-section: no TALYS graphs in root_files (run talys-xs); "
                 "midpoint energies, no overlay"
              << std::endl;

  std::vector<ChannelResult> results(X.CHANNELS.size());
  std::vector<const ChannelResult *> done;
  for (Int_t c = 0; c < Int_t(X.CHANNELS.size()); c++)
    if (RunChannel(X.CHANNELS[c], results[c]))
      done.push_back(&results[c]);
  if (done.empty())
    return kFALSE;
  if (done.size() > 1)
    Draw(done, "cross_section" + TagSuffix());
  cache_->Close();
  return kTRUE;
}
