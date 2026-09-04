#include "CrossSection.hpp"
#include "BeamEnergies.hpp"
#include "Constants.hpp"
#include "InitUtils.hpp"
#include "Paths.hpp"
#include "PlottingUtils.hpp"
#include "RegionCuts.hpp"
#include "StripSumScatter.hpp"
#include "TagEfficiency.hpp"
#include <TAxis.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TKey.h>
#include <TLegend.h>
#include <TMath.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TSystem.h>
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

// Marker styles and colours for the channels on a combined figure.
const Int_t kChannelMarker[4] = {20, 21, 22, 23};
const Int_t kChannelColor[4] = {kBlack, kGreen + 2, kMagenta + 2, kOrange + 7};

} // namespace

Double_t CrossSection::Point::Err() const {
  return std::sqrt(stat * stat + sys * sys);
}

// "n", "2n", "pn", "a", "g": light particles leaving the compound nucleus,
// each with an optional multiplicity digit in front. The residue is the
// compound (beam + alpha) minus what left.
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

Long64_t CrossSection::ReadCount(TFile &f, const char *name, Bool_t &ok) {
  TParameter<Long64_t> *p = static_cast<TParameter<Long64_t> *>(f.Get(name));
  if (!p) {
    ok = kFALSE;
    return 0;
  }
  return p->GetVal();
}

// The curve clipped to an energy window, for drawing.
TGraph *CrossSection::Clipped(TGraph *g, Double_t e_lo, Double_t e_hi) {
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
Double_t CrossSection::EffectiveEnergy(TGraph *axn, Double_t e_out,
                                       Double_t e_in) {
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

// Events in the scatter whose bin centre falls inside the cut scaled by
// `scale`. Counting at a few scales and correcting each by the fraction of a
// Gaussian it should enclose measures how much of what the region holds is
// not the component: a clean peak gives the same number at every scale,
// contamination grows with area.
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

// The simulated beam's loss per strip and the energy entering strip 0, and
// from them the centre-of-mass energy at the middle of every strip, so the
// alignment to a published table can be found over the whole detector rather
// than only the strips a cross section is reported for.
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

// Every model in root_files/talys/talys_xs.root, in the config's order, with
// every residual-production graph talys-xs wrote for it. Empty when there is
// no file; the plot goes on without them. The label is presentation, so the
// config's current one wins over the one stamped in the file; relabelling
// never needs a TALYS rerun.
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
CrossSection::ChannelCurves(const CrossSectionChannel &ch) const {
  const CrossSectionConfig &X = Constants::cfg.CROSS_SECTION_CONFIG;
  std::vector<TalysCurve> curves;
  if (talys_raw_.empty())
    return curves;
  for (Int_t m = 0; m < Int_t(talys_raw_.size()); m++) {
    std::map<Double_t, Double_t> sum;
    TString used, missing;
    for (Int_t k = 0; k < Int_t(ch.talys_exits.size()); k++) {
      Int_t z = 0, a = 0;
      ExitResidue(ch.talys_exits[k], X.BEAM_Z, X.BEAM_A, z, a);
      std::map<std::pair<Int_t, Int_t>, TGraph *>::const_iterator it =
          talys_raw_[m].find(std::make_pair(z, a));
      if (it == talys_raw_[m].end()) {
        missing += Form(" %s(Z=%d,A=%d)", ch.talys_exits[k].Data(), z, a);
        continue;
      }
      for (Int_t p = 0; p < it->second->GetN(); p++)
        sum[it->second->GetX()[p]] += it->second->GetY()[p];
      used += Form(" %s(Z=%d,A=%d)", ch.talys_exits[k].Data(), z, a);
    }
    std::vector<Double_t> tx, ty;
    for (std::map<Double_t, Double_t>::const_iterator it = sum.begin();
         it != sum.end(); ++it)
      if (it->second > 0.0) {
        tx.push_back(it->first);
        ty.push_back(it->second);
      }
    std::cout << "cross-section: " << ch.name << ": " << talys_labels_[m]
              << ": summed" << (used.Length() ? used : " nothing")
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

Bool_t CrossSection::Strip(const CrossSectionChannel &ch,
                           const std::vector<TalysCurve> &talys, Int_t reac,
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
  // Beam particles that reached this strip under the same conditions a
  // reaction here had to satisfy. Falls back to the flat count only if the
  // cache predates the per-strip one.
  Bool_t at_ok = kTRUE;
  const Long64_t n_at = ReadCount(*cache_, Form("n_normed_r%d", reac), at_ok);
  pt.n_denom = Double_t(at_ok && n_at > 0 ? n_at : n_beam_);

  // Beam energy at the strip's entrance and exit: what entered strip 0, less
  // each earlier strip's mean loss. Effective energy from the first model's
  // shape; the others say how much it depends on the shape.
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
  TagEfficiencyRecord eff;
  if (TagEfficiencyStore::Load(ch.name, reac, eff) && eff.eff > 0.0) {
    // The efficiency side's own count, unfolded: what the strip before fed
    // into it by migration is removed first, then the loss is undone. The
    // feed is only known when the strip before was itself unfolded.
    const Double_t feed =
        prev_reac_ == reac - 1 ? prev_migrate_ * prev_true_ : 0.0;
    pt.n_reac = (eff.n_counted - feed) / eff.eff;
    pt.sigma = pt.n_reac / norm;
    pt.stat = eff.n_counted > 0.0 ? pt.sigma / std::sqrt(eff.n_counted) : 0.0;
    pt.sys = pt.sigma * eff.eff_err / eff.eff;
    prev_reac_ = reac;
    prev_true_ = pt.n_reac;
    prev_migrate_ = eff.migrate;
    std::cout << Form("   %2d    [%5.2f, %5.2f]      %6.2f   %6.0f  %9.0f   "
                      "%7.1f +- %.1f (%.1f stat, %.1f eff; counted %.0f, "
                      "eff %.3f, fed %.0f from strip %d)",
                      reac, pt.e_in, pt.e_out, pt.e_eff, pt.n_reac, pt.n_denom,
                      pt.sigma, pt.Err(), pt.stat, pt.sys, eff.n_counted,
                      eff.eff, feed, reac - 1)
              << std::endl;
    delete cut;
    return kTRUE;
  }
  prev_reac_ = -1;

  // No efficiency record: two estimates of the reaction count. The geometric
  // one is everything inside the region, corrected for the fraction of a
  // Gaussian it should enclose; the attributed one is what the mixture fit
  // assigned to the reaction component (trimmed at 3 sigma, so corrected for
  // that). Where the island is well off the beam ridge they agree; where the
  // region also covers beam tail the geometric count runs away, since the
  // beam's real tail is far heavier than any Gaussian and cannot be
  // subtracted by model. The attributed count is the value; the two
  // estimates' disagreement is the region systematic, because that
  // disagreement is exactly the overlap ambiguity.
  const Double_t n_raw = CountInCut(h, cut, 1.0);
  const Double_t est_geo = n_raw / Enclosed(C.AN_REGION_NSIGMA);
  const Double_t n_assigned = RegionCutStore::LoadAssigned(region, reac);
  const Bool_t have_fit = n_assigned >= 0.0;
  const Double_t est_fit = have_fit ? n_assigned / Enclosed(3.0) : est_geo;
  pt.n_reac = est_fit;
  pt.sigma = pt.n_reac / norm;
  pt.stat = pt.n_reac > 0.0 ? pt.sigma / std::sqrt(pt.n_reac) : 0.0;
  if (have_fit) {
    // Against the region's core: the 1-sigma ellipse (half the drawn region),
    // corrected for the 39% it encloses. The full region's count is no check
    // where it also covers tail, but the core is where the island dominates,
    // so core and attributed counts should agree, and their disagreement is
    // the honest size of the overlap ambiguity.
    const Double_t est_core =
        CountInCut(h, cut, 1.0 / C.AN_REGION_NSIGMA) / Enclosed(1.0);
    pt.sys = std::fabs(est_core - est_fit) / norm;
  } else {
    // A drawn region without an attributed count: scale it by 0.75x and
    // 1.25x, each corrected for its enclosed fraction.
    Double_t lo = pt.sigma, hi = pt.sigma;
    const Double_t kScale[2] = {0.75, 1.25};
    for (Int_t k = 0; k < 2; k++) {
      const Double_t s = CountInCut(h, cut, kScale[k]) /
                         Enclosed(kScale[k] * C.AN_REGION_NSIGMA) / norm;
      lo = TMath::Min(lo, s);
      hi = TMath::Max(hi, s);
    }
    pt.sys = 0.5 * (hi - lo);
  }
  std::cout << Form("   %2d    [%5.2f, %5.2f]      %6.2f   %6.0f  %9.0f   "
                    "%7.1f +- %.1f (%.1f stat, %.1f region; in-region %.0f, "
                    "attributed %s)",
                    reac, pt.e_in, pt.e_out, pt.e_eff, pt.n_reac, pt.n_denom,
                    pt.sigma, pt.Err(), pt.stat, pt.sys, est_geo,
                    have_fit ? Form("%.0f", est_fit) : "none")
            << std::endl;
  delete cut;
  return kTRUE;
}

// Against the published values. A published table is a list of strips, so
// the two are aligned strip to row by the single integer offset that best
// matches the energies over the whole table, not row by row: pairing each
// strip with its nearest published energy silently pairs two strips with one
// row when a dataset reaches energies the reference never reported.
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
    std::cout << Form("   %2d     %6.2f    %6.2f        %7.1f +- %-5.1f  "
                      "%7.1f (%.1f)   %5.2f",
                      pt.reac, pt.e_eff, ref[row][0], pt.sigma, pt.Err(), v,
                      ref[row][4], v > 0.0 ? pt.sigma / v : 0.0)
              << std::endl;
  }
}

// Excitation function(s), this work against the published points. The frame
// spans every set drawn so the measured strips are seen in the context of the
// whole published curve, not just the rows they happen to sit beside. One
// channel: its own title. Several: the dataset alone, every channel labelled.
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
    std::vector<Double_t> vx, vy, vexl, vexh, vey;
    for (Int_t i = 0; i < Int_t(r.points.size()); i++) {
      const Point &pt = r.points[i];
      // The energy error is the strip's extent about the effective energy,
      // asymmetric since that energy sits above the midpoint on a rising
      // curve (the reference table's convention), with the effective
      // energy's dependence on the model shape added in quadrature on each
      // side.
      vx.push_back(pt.e_eff);
      vy.push_back(pt.sigma);
      vexl.push_back(std::hypot(pt.e_eff - TMath::Min(pt.e_in, pt.e_out),
                                pt.e_eff - pt.e_eff_lo));
      vexh.push_back(std::hypot(TMath::Max(pt.e_in, pt.e_out) - pt.e_eff,
                                pt.e_eff_hi - pt.e_eff));
      vey.push_back(pt.Err());
      fx_lo = TMath::Min(fx_lo, vx.back() - vexl.back());
      fx_hi = TMath::Max(fx_hi, vx.back() + vexh.back());
      fy_lo = TMath::Min(fy_lo, vy.back());
      fy_hi = TMath::Max(fy_hi, vy.back());
    }
    if (vx.empty())
      continue;
    Series s;
    s.g = new TGraphAsymmErrors(Int_t(vx.size()), &vx[0], &vy[0], &vexl[0],
                                &vexh[0], &vey[0], &vey[0]);
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
    // Hauser-Feshbach prediction, if talys-xs has written one. The first
    // model in full, the shape checks dashed; on a combined figure only the
    // first model per channel, in the channel's colour.
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
    const Int_t ch = rs.size() > 1 ? k : 0;
    g->SetLineColor(rs.size() > 1 ? kChannelColor[ch % 4] : kAzure + 1);
    g->SetLineWidth(k == 0 || rs.size() > 1 ? 2 : 1);
    g->SetLineStyle(k == 0 || rs.size() > 1 ? 1 : 2);
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
      PlottingUtils::AddLegend(0.42, 0.89, 0.16, 0.16 + 0.07 * n_entries);
  for (Int_t k = 0; k < Int_t(measured.size()); k++)
    leg->AddEntry(measured[k].g, measured[k].label, "pe");
  for (Int_t k = 0; k < Int_t(published.size()); k++)
    leg->AddEntry(published[k].g, published[k].label, "pe");
  for (Int_t k = 0; k < Int_t(drawn.size()); k++)
    leg->AddEntry(drawn[k].first, drawn[k].second, "l");
  leg->Draw();
  PlottingUtils::SaveFigure(c, name, "cross_section", PlotSaveOptions::kLOG);
  delete c;
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
  if (!talys_raw_.empty() && ch.talys_exits.empty()) {
    std::cerr << "cross-section: channel " << ch.name
              << ": TALYS models are declared but the channel names no exits"
              << std::endl;
    return kFALSE;
  }
  out.talys = ChannelCurves(ch);
  std::cout << std::endl
            << "  " << ch.name << " " << Label(ch) << ", region_" << ch.name
            << (out.talys.empty() ? " (no TALYS curve: midpoint energies)" : "")
            << std::endl;
  std::cout << "  strip   E_cm range [MeV]    "
            << (X.EFFECTIVE_ENERGY ? "E_cm,eff" : "E_cm,mid")
            << "   N_reac    N_beam    sigma [mb]" << std::endl;
  prev_reac_ = -1;
  for (Int_t reac = X.XS_STRIP_MIN; reac <= X.XS_STRIP_MAX; reac++) {
    Point pt;
    if (Strip(ch, out.talys, reac, pt))
      out.points.push_back(pt);
  }
  if (out.points.empty()) {
    std::cerr << "cross-section: channel " << ch.name
              << ": no strip produced a cross section" << std::endl;
    return kFALSE;
  }
  CompareReference(out);
  std::vector<const ChannelResult *> one(1, &out);
  Draw(one, "cross_section_" + ch.name);
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
  const TString method = TagEfficiencyStore::Method();
  if (method.Length() > 0)
    std::cout << "  tag efficiency: " << method
              << " (strips with a record are unfolded by it)" << std::endl;
  else
    std::cout << Form("  no tag-efficiency store; regions at %.1f sigma "
                      "enclose %.1f%% of the fitted component and counts "
                      "are corrected for that only",
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
    Draw(done, "cross_section");
  cache_->Close();
  return kTRUE;
}
