#include "PulseHistory.hpp"
#include "Constants.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TDecompSVD.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResultPtr.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TProfile.h>
#include <TTree.h>
#include <TVectorD.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <deque>
#include <functional>
#include <iostream>
#include <map>

namespace PulseHistory {

const char *GroupName(Int_t g) {
  switch (g) {
  case kLongLeft:
    return "long L (odd strips)";
  case kLongRight:
    return "long R (even strips)";
  case kShortLeft:
    return "short L (even strips)";
  case kShortRight:
    return "short R (odd strips)";
  case kStrip0:
    return "strip 0 (unsegmented)";
  case kStrip17:
    return "strip 17 (unsegmented)";
  }
  return "none";
}

const char *GroupTag(Int_t g) {
  switch (g) {
  case kLongLeft:
    return "L";
  case kLongRight:
    return "R";
  case kShortLeft:
    return "Ls";
  case kShortRight:
    return "Rs";
  case kStrip0:
    return "S0";
  case kStrip17:
    return "S17";
  }
  return "none";
}

Int_t ChainOf(Int_t g) {
  switch (g) {
  case kLongLeft:
  case kShortLeft:
    return 0;
  case kLongRight:
  case kShortRight:
    return 1;
  }
  return -1;
}

Bool_t IsLongGroup(Int_t g) { return g == kLongLeft || g == kLongRight; }

namespace {
const PulseHistoryGroupOption *GroupOption(Int_t g) {
  const PulseHistoryGroups &m = Constants::cfg.PULSE_HISTORY_GROUPS;
  switch (g) {
  case kLongLeft:
    return &m.long_left;
  case kLongRight:
    return &m.long_right;
  case kShortLeft:
    return &m.short_left;
  case kShortRight:
    return &m.short_right;
  case kStrip0:
    return &m.strip0;
  case kStrip17:
    return &m.strip17;
  default:
    return nullptr;
  }
}
} // namespace

Bool_t GroupEnabled(Int_t g) {
  const PulseHistoryGroupOption *o = GroupOption(g);
  return o ? o->enabled : kFALSE;
}

Bool_t GroupWantsForm(Int_t g) {
  const PulseHistoryGroupOption *o = GroupOption(g);
  return o ? (o->enabled && o->kernel == kPulseHistoryForm) : kFALSE;
}

Double_t GroupTauUs(Int_t g) {
  const PulseHistoryGroupOption *o = GroupOption(g);
  return (o && o->tau_us > 0.0) ? o->tau_us : 0.0;
}

Kernel::Kernel() {
  for (Int_t a = 0; a < kMaxAmpBins; a++) {
    c[a] = 0.0;
    cb[a] = 0.0;
    for (Int_t b = 0; b < kNBins; b++) {
      k[a][b] = 0.0;
      k_binned[a][b] = 0.0;
      k_low[a][b] = 0.0;
    }
  }
}

Double_t Kernel::FormAt(Int_t a, Double_t dt_us) const {
  if (a < 0 || a >= kMaxAmpBins || !(tau_us > 0.0))
    return 0.0;
  if (dt_us < t_lo_us) {
    // The previous trapezoid is still under the read point: free bins.
    const Int_t b = BinOf(dt_us * 1.0e-6);
    return (b >= 0 && b < n_free) ? k_low[a][b] : 0.0;
  }
  const Double_t ex = TMath::Exp(-dt_us / tau_us);
  Double_t v = -c[a] * ex * TMath::Exp(-t_m_us / tau_us);
  if (dt_us > t_f_us)
    v += cb[a] * (TMath::Exp(-t_f_us / tau_us) - ex);
  return v;
}

Double_t Kernel::EffectiveWindowUs(Int_t a) const {
  if (a < 0 || a >= kMaxAmpBins || cb[a] == 0.0)
    return 0.0;
  return c[a] * tau_us / cb[a];
}

Result::Result() {
  for (Int_t g = 0; g < kNGroups; g++) {
    enabled[g] = kTRUE;
    mean_shift[g] = 0.0;
    n_clamped_group[g] = 0;
    dev_vs_pred[g] = nullptr;
    dev_before[g] = dev_after[g] = shift[g] = nullptr;
    dtprev_before[g] = dtprev_after[g] = nullptr;
  }
}

Result::~Result() { FreeDiagnostics(); }

void Result::FreeDiagnostics() {
  for (Int_t g = 0; g < kNGroups; g++) {
    delete dev_vs_pred[g];
    delete dev_before[g];
    delete dev_after[g];
    delete shift[g];
    delete dtprev_before[g];
    delete dtprev_after[g];
    dev_vs_pred[g] = nullptr;
    dev_before[g] = dev_after[g] = shift[g] = nullptr;
    dtprev_before[g] = dtprev_after[g] = nullptr;
  }
  for (size_t c = 0; c < dtprev_before_ch.size(); c++)
    delete dtprev_before_ch[c];
  for (size_t c = 0; c < dtprev_after_ch.size(); c++)
    delete dtprev_after_ch[c];
  dtprev_before_ch.clear();
  dtprev_after_ch.clear();
}

Int_t BinOf(Double_t dt_s) {
  if (!(dt_s > 0.0))
    return -1;
  const Double_t l = TMath::Log10(dt_s);
  if (l < kLogLo || l >= kLogHi)
    return -1;
  return Int_t((l - kLogLo) / (kLogHi - kLogLo) * kNBins);
}

Double_t BinCentreUs(Int_t b) {
  return TMath::Power(10.0, kLogLo + (b + 0.5) * (kLogHi - kLogLo) / kNBins) *
         1.0e6;
}

/// p0 + p1 exp(-dt / p2) + p3 dt through (dt us, deviation ADC) pairs, the
/// form the pole-zero residual plus the slow baseline lift takes. A p1 within
/// 2 sigma of zero, or a p2 pinned at a limit, means the profile has no decay
/// to measure: the group's pole-zero is matched, and a decay time off such a
/// fit is noise.
void FitTauDecay(const std::vector<Double_t> &dt_us,
                 const std::vector<Double_t> &dev, Kernel::TauFit &out) {
  out = Kernel::TauFit{};
  const Int_t n = Int_t(dt_us.size());
  if (n < kTauFitMinEntries)
    return;
  // Binned into a TProfile: its bin errors weight populated early bins, where
  // the tail lives, against the sparse long tail, ROOT skips empty bins.
  const Int_t nb = 60;
  TProfile *prof = new TProfile("ph_tau_fit", "", nb, kTauFitLoUs, kTauFitHiUs);
  for (Int_t i = 0; i < n; i++)
    prof->Fill(dt_us[i], dev[i]);
  TF1 *f = new TF1("ph_tau_fit_fn", "[0] + [1] * exp(-x / [2]) + [3] * x",
                   kTauFitLoUs, kTauFitHiUs);
  const Double_t head = prof->GetBinContent(prof->FindBin(kTauFitLoUs + 1.0));
  const Double_t tail = prof->GetBinContent(prof->FindBin(kTauFitHiUs - 5.0));
  f->SetParameters(tail, head - tail, 10.0, 0.0);
  f->SetParLimits(2, 0.5, 200.0);
  TFitResultPtr r = prof->Fit(f, "RQS0");
  if (Int_t(r) == 0) {
    out.ok = kTRUE;
    out.p0 = f->GetParameter(0);
    out.p1 = f->GetParameter(1);
    out.p1err = f->GetParError(1);
    out.p2 = f->GetParameter(2);
    out.p2err = f->GetParError(2);
    out.p3 = f->GetParameter(3);
    out.chi2 = f->GetChisquare();
    out.ndf = f->GetNDF();
    out.flat = TMath::Abs(out.p1) < 2.0 * out.p1err || out.p2 <= 0.51 ||
               out.p2 >= 199.0;
  }
  // "S" gives the histogram ownership of the fit function; mark it
  // referenced so deleting the profile does not delete it under us.
  f->SetBit(TF1::kNotDraw);
  delete f;
  delete prof;
}

Int_t AmpBinOf(Double_t e_prev, Double_t mode, Int_t n_amp) {
  if (n_amp <= 1 || !(mode > 0.0))
    return 0;
  // Bands centred on multiples of the beam pulse: [0, 0.5), [0.5, 1.5), ...
  const Int_t a = Int_t(e_prev / mode + 0.5);
  return TMath::Min(a, n_amp - 1);
}

std::vector<Int_t> BuildGroupMap() {
  const Int_t nch = Constants::ActiveNChannels();
  std::vector<Int_t> gm(Constants::ActiveNBoards() * nch, kNone);
  const std::map<std::pair<Int_t, Int_t>, TString> &cm =
      Constants::ActiveChannelMap();
  for (std::map<std::pair<Int_t, Int_t>, TString>::const_iterator it =
           cm.begin();
       it != cm.end(); ++it) {
    const TString &name = it->second;
    const Int_t idx = it->first.first * nch + it->first.second;
    if (idx < 0 || idx >= Int_t(gm.size()))
      continue;
    if (name == "Strip0") {
      gm[idx] = kStrip0;
      continue;
    }
    if (name == "Strip17") {
      gm[idx] = kStrip17;
      continue;
    }
    if (name.Length() < 2 || (name[0] != 'L' && name[0] != 'R'))
      continue;
    TString num = name(1, name.Length() - 1);
    if (!num.IsDigit())
      continue;
    const Int_t s = num.Atoi();
    if (s < 1 || s > 16)
      continue;
    const Bool_t odd = (s % 2) != 0;
    const Bool_t is_long = (name[0] == 'L') ? odd : !odd;
    if (name[0] == 'L')
      gm[idx] = is_long ? kLongLeft : kShortLeft;
    else
      gm[idx] = is_long ? kLongRight : kShortRight;
  }
  return gm;
}

namespace {

// The channel map's name for a (board, channel) index; empty if unmapped.
TString NameOfIndex(Int_t idx) {
  const Int_t nch = Constants::ActiveNChannels();
  const std::map<std::pair<Int_t, Int_t>, TString> &cm =
      Constants::ActiveChannelMap();
  std::map<std::pair<Int_t, Int_t>, TString>::const_iterator it =
      cm.find(std::make_pair(idx / nch, idx % nch));
  return it == cm.end() ? TString("") : it->second;
}

Int_t IndexOfName(const TString &want) {
  const Int_t nch = Constants::ActiveNChannels();
  const std::map<std::pair<Int_t, Int_t>, TString> &cm =
      Constants::ActiveChannelMap();
  for (std::map<std::pair<Int_t, Int_t>, TString>::const_iterator it =
           cm.begin();
       it != cm.end(); ++it)
    if (it->second == want)
      return it->first.first * nch + it->first.second;
  return -1;
}

inline Int_t HitIndex(const RawHit &h) {
  return Int_t(h.board) * Constants::ActiveNChannels() + Int_t(h.channel);
}

struct Past {
  Double_t t; // ps
  Double_t e; // raw ADC
};

// Feature layout of the fit: 0 is the intercept; then the binned kernel's
// (band, dt bin); then, when the form is fitted too, the form's free low-dt
// bins (band, bin < n_free) and, per tau of the grid and band, the
// direct-tail and the baseline feature. One scan fills all of them and each
// model is solved from its own block of the normal matrix.
struct FeatLayout {
  Int_t n_amp = 1;
  Int_t n_free = 0;
  Int_t n_tau = 0;
  Int_t Bin(Int_t a, Int_t b) const { return 1 + a * kNBins + b; }
  Int_t Low(Int_t a, Int_t b) const {
    return 1 + n_amp * kNBins + a * n_free + b;
  }
  Int_t D(Int_t j, Int_t a) const {
    return 1 + n_amp * kNBins + n_amp * n_free + 2 * (j * n_amp + a);
  }
  Int_t B(Int_t j, Int_t a) const { return D(j, a) + 1; }
  Int_t Size() const {
    return 1 + n_amp * kNBins + n_amp * n_free + 2 * n_tau * n_amp;
  }
  Bool_t HasForm() const { return n_tau > 0; }
};

// The tau grid of the form and the timings from the trapezoid settings, with
// the two exponentials that do not depend on the gap precomputed per tau.
struct FormGrid {
  Double_t t_m_us = 0.0, t_f_us = 0.0, t_lo_us = 0.0;
  std::vector<Double_t> tau_us, em, ef; // tau, exp(-t_m/tau), exp(-t_f/tau)
  // exp(-dt/tau_j) tabulated against log10(dt [us]) over the kernel reach,
  // linearly interpolated: one table read per tau and past pulse in the scan
  // instead of one exp(). 4096 points over 2.5 decades put the interpolation
  // error below 1e-5, far inside the fit's own precision.
  static const Int_t kTabN = 4096;
  std::vector<Double_t> tab; // [j * kTabN + i]
  Double_t tab_lo = 0.0, tab_step = 0.0;
  void Tabulate() {
    tab_lo = kLogLo + 6.0;
    tab_step = (kLogHi - kLogLo) / Double_t(kTabN - 1);
    tab.assign(tau_us.size() * kTabN, 0.0);
    for (size_t j = 0; j < tau_us.size(); j++)
      for (Int_t i = 0; i < kTabN; i++)
        tab[j * kTabN + i] =
            TMath::Exp(-TMath::Power(10.0, tab_lo + i * tab_step) / tau_us[j]);
  }
  // exp(-dt/tau_j) at log10(dt [us]) = l.
  inline Double_t Exp(Int_t j, Double_t l) const {
    Double_t u = (l - tab_lo) / tab_step;
    if (u <= 0.0)
      return tab[j * kTabN];
    if (u >= kTabN - 1)
      return tab[j * kTabN + kTabN - 1];
    const Int_t i = Int_t(u);
    const Double_t f = u - i;
    return tab[j * kTabN + i] * (1.0 - f) + tab[j * kTabN + i + 1] * f;
  }
};

// Which features a scan has to build. The fit needs every one; the
// diagnostics only the applied model's, and the tau summary none.
struct ScanNeeds {
  Bool_t bins = kTRUE;
  Bool_t low = kTRUE;
  std::vector<Int_t> taus; // grid indices whose form features are built
};

// One accumulated least-squares problem: the right-hand side, the sums that
// give the residual, and the normal matrix in a plain array, of which only
// the blocks a fit reads are accumulated (see the scan in Measure): the
// binned block, the free low bins, and per tau its two features against the
// low bins and each other. Cross terms between the bins and the form, or
// between two taus, are never read, and the upper triangle is enough.
struct Normal {
  Int_t np = 0;
  std::vector<Double_t> A; // np x np, upper triangle of the read blocks
  std::vector<Double_t> b;
  Double_t syy = 0.0, sy = 0.0;
  Long64_t n = 0;
  void Init(Int_t np_) {
    np = np_;
    A.assign(size_t(np) * np, 0.0);
    b.assign(np, 0.0);
    syy = sy = 0.0;
    n = 0;
  }
  void Add(const Normal &o) {
    for (size_t i = 0; i < A.size(); i++)
      A[i] += o.A[i];
    for (size_t i = 0; i < b.size(); i++)
      b[i] += o.b[i];
    syy += o.syy;
    sy += o.sy;
    n += o.n;
  }
  inline Double_t &At(Int_t a, Int_t c) { return A[size_t(a) * np + c]; }
  inline Double_t Get(Int_t a, Int_t c) const {
    return a <= c ? A[size_t(a) * np + c] : A[size_t(c) * np + a];
  }
};

// Predicted shift from a channel's past, with the group's applied kernel.
Double_t Shift(const std::deque<Past> &d, Double_t t, const Kernel &K,
               Double_t mode) {
  Double_t s = 0.0;
  for (size_t k = 0; k < d.size(); k++) {
    const Past &p = d[k];
    const Double_t dt_s = (t - p.t) * 1.0e-12;
    const Int_t b = BinOf(dt_s);
    if (b < 0)
      continue;
    const Int_t a = AmpBinOf(p.e, mode, K.n_amp);
    s += (K.form ? K.FormAt(a, dt_s * 1.0e6) : K.k[a][b]) * p.e;
  }
  return s;
}

// The applied model's prediction from one feature vector (the diagnostics).
Double_t Predict(const Kernel &K, const FeatLayout &lay,
                 const std::vector<Double_t> &x) {
  Double_t pred = K.intercept;
  for (Int_t a = 0; a < K.n_amp; a++) {
    if (K.form) {
      for (Int_t b = 0; b < lay.n_free; b++)
        pred += K.k_low[a][b] * x[lay.Low(a, b)];
      if (K.j_tau >= 0)
        pred += -K.c[a] * x[lay.D(K.j_tau, a)] + K.cb[a] * x[lay.B(K.j_tau, a)];
    } else {
      for (Int_t b = 0; b < kNBins; b++)
        pred += K.k_binned[a][b] * x[lay.Bin(a, b)];
    }
  }
  return pred;
}

// Least squares on a subset of the features from the full normal matrix:
// scaled SVD (far-dt bins of rare bands are near-collinear with the
// intercept, and the tau grid's neighbours with each other). Returns the
// coefficients in the order of `keep` and the drop in the residual sum of
// squares, sum p b, so ss_res = syy - drop.
Bool_t SolveSubset(const Normal &N, const std::vector<Int_t> &keep, TVectorD &p,
                   Double_t &drop) {
  const Int_t nk = Int_t(keep.size());
  if (nk < 1)
    return kFALSE;
  TMatrixD Ak(nk, nk);
  TVectorD bk(nk);
  for (Int_t a = 0; a < nk; a++) {
    bk[a] = N.b[keep[a]];
    for (Int_t b = 0; b < nk; b++)
      Ak(a, b) = N.Get(keep[a], keep[b]);
  }
  TVectorD scale(nk);
  for (Int_t a = 0; a < nk; a++)
    scale[a] = Ak(a, a) > 0.0 ? 1.0 / TMath::Sqrt(Ak(a, a)) : 1.0;
  TMatrixD As(nk, nk);
  TVectorD bs(nk);
  for (Int_t a = 0; a < nk; a++) {
    bs[a] = bk[a] * scale[a];
    for (Int_t b = 0; b < nk; b++)
      As(a, b) = Ak(a, b) * scale[a] * scale[b];
  }
  TDecompSVD svd(As);
  svd.SetTol(1.0e-10);
  Bool_t solved = kFALSE;
  TVectorD ps = svd.Solve(bs, solved);
  if (!solved)
    return kFALSE;
  p.ResizeTo(nk);
  drop = 0.0;
  for (Int_t a = 0; a < nk; a++) {
    p[a] = ps[a] * scale[a];
    drop += p[a] * bk[a];
  }
  return std::isfinite(drop);
}

// The hit window and channel maps shared by the seed pass and the scan: the
// per-channel maxima of the current event window, refilled at every seed.
struct ScanState {
  const std::vector<RawHit> &hits;
  const std::vector<Int_t> &group_of;
  const std::vector<Int_t> &fit_ch;
  const std::vector<Int_t> &long_ch;
  const std::vector<Double_t> &mode;
  std::vector<Double_t> &ev_e;
  const std::vector<Double_t> &mean;
  Int_t nidx;
  Int_t strip0;
  Double_t window_ps;
  Double_t olo;
  Double_t ohi;
  Double_t blo;
  Double_t bhi;
  Int_t n_amp;
  Double_t keep_ps;
  const FeatLayout &lay;
  const FormGrid &grid;
};

// Max hit energy per channel within the event window from seed i0; false if
// strip 0 did not fire.
Bool_t Lookahead(const ScanState &st, size_t i0) {
  for (Int_t k = 0; k < Int_t(st.fit_ch.size()); k++)
    st.ev_e[st.fit_ch[k]] = 0.0;
  Bool_t strip0_fired = st.strip0 < 0;
  const ULong64_t t0 = st.hits[i0].timestamp;
  for (size_t j = i0 + 1;
       j < st.hits.size() && st.hits[j].timestamp - t0 < st.window_ps; j++) {
    const Int_t i = HitIndex(st.hits[j]);
    if (i == st.strip0 && st.hits[j].energy > 0)
      strip0_fired = kTRUE;
    if (i < 0 || i >= st.nidx || st.group_of[i] == kNone)
      continue;
    if (Double_t(st.hits[j].energy) > st.ev_e[i])
      st.ev_e[i] = Double_t(st.hits[j].energy);
  }
  return strip0_fired;
}

// Beam-like for group g: every long end inside its window around its mode,
// g's own chain loose, the other chain and strips 0/17 tight.
Bool_t BeamFor(const ScanState &st, Int_t g) {
  const Int_t chain = ChainOf(g);
  for (Int_t k = 0; k < Int_t(st.long_ch.size()); k++) {
    const Int_t c = st.long_ch[k];
    if (!(st.mode[c] > 0.0))
      return kFALSE;
    const Double_t r = st.ev_e[c] / st.mode[c];
    const Bool_t own = chain >= 0 && ChainOf(st.group_of[c]) == chain;
    if (own ? (r < st.olo || r > st.ohi) : (r < st.blo || r > st.bhi))
      return kFALSE;
  }
  return kTRUE;
}

/// One walk over the hits for group g, calling fn(channel, features x,
/// deviation y, dt to the previous pulse) at every beam-like seed of g. The
/// per-channel deque holds earlier pulses only, so the seed's own hit never
/// enters.
void Scan(const ScanState &st, Int_t g, const std::vector<size_t> &seeds_g,
          const ScanNeeds &needs, std::vector<Double_t> &x,
          const std::function<void(Int_t, const std::vector<Double_t> &,
                                   Double_t, Double_t)> &fn) {
  const Bool_t any_form =
      st.lay.HasForm() && (needs.low || !needs.taus.empty());
  std::vector<std::deque<Past>> past(st.nidx);
  size_t next_seed = 0;
  for (size_t j = 0; j < st.hits.size(); j++) {
    if (next_seed < seeds_g.size() && seeds_g[next_seed] == j) {
      next_seed++;
      const Double_t tg = Double_t(st.hits[j].timestamp);
      Lookahead(st, j);
      for (Int_t k = 0; k < Int_t(st.fit_ch.size()); k++) {
        const Int_t c = st.fit_ch[k];
        // A short end or unsegmented strip that did not fire has no height.
        if (st.group_of[c] != g || !(st.ev_e[c] > 0.0))
          continue;
        std::deque<Past> &d = past[c];
        while (!d.empty() && d.front().t < tg - st.keep_ps)
          d.pop_front();
        std::fill(x.begin(), x.end(), 0.0);
        x[0] = 1.0;
        for (size_t p = 0; p < d.size(); p++) {
          const Double_t dt_s = (tg - d[p].t) * 1.0e-12;
          if (!(dt_s > 0.0))
            continue;
          // One log per past pulse serves the bin and the exp table.
          const Double_t l_s = TMath::Log10(dt_s);
          if (l_s < kLogLo || l_s >= kLogHi)
            continue; // outside the kernel reach
          const Int_t b = Int_t((l_s - kLogLo) / (kLogHi - kLogLo) * kNBins);
          const Int_t a = AmpBinOf(d[p].e, st.mode[c], st.n_amp);
          const Double_t e = d[p].e;
          if (needs.bins)
            x[st.lay.Bin(a, b)] += e;
          if (!any_form)
            continue;
          const Double_t dt_us = dt_s * 1.0e6;
          if (dt_us < st.grid.t_lo_us) {
            if (needs.low && b < st.lay.n_free)
              x[st.lay.Low(a, b)] += e;
            continue;
          }
          if (needs.taus.empty())
            continue;
          const Double_t l = l_s + 6.0;
          const Bool_t past_freeze = dt_us > st.grid.t_f_us;
          for (size_t q = 0; q < needs.taus.size(); q++) {
            const Int_t j = needs.taus[q];
            const Double_t ex = st.grid.Exp(j, l);
            x[st.lay.D(j, a)] += e * ex * st.grid.em[j];
            if (past_freeze)
              x[st.lay.B(j, a)] += e * (st.grid.ef[j] - ex);
          }
        }
        const Double_t dt_prev = d.empty() ? -1.0 : (tg - d.back().t) * 1.0e-12;
        fn(c, x, st.ev_e[c] - st.mean[c], dt_prev);
      }
    }
    const Int_t i = HitIndex(st.hits[j]);
    if (i >= 0 && i < st.nidx && st.group_of[i] == g)
      past[i].push_back(
          {Double_t(st.hits[j].timestamp), Double_t(st.hits[j].energy)});
  }
}

// Fits one kernel from an accumulated problem: always the binned model and,
// when the layout carries the form's features, the form as well, at every tau
// of the grid (the profile) and at the given tau `given_j` when there is one
// (-1 when the tau is to be profiled). Then picks the applied one: the form
// where `want_form` and it fitted, the bins otherwise. K.j_tau records the
// grid index of the applied tau for the diagnostics. Returns K.ok; `why`
// receives a one-line reason when the given tau could not be used.
// R^2 that a kernel fitted elsewhere (the group's) reaches on this problem:
// the residual of its coefficient vector on the channel's own normal
// equations, so a channel's own fit can be judged against the fallback.
Double_t R2Of(const Kernel &K, const Normal &N, const FeatLayout &lay) {
  if (N.n < 1)
    return 0.0;
  const Double_t ybar = N.sy / N.n;
  const Double_t ss_tot = N.syy - N.n * ybar * ybar;
  if (!(ss_tot > 0.0))
    return 0.0;
  TVectorD p(lay.Size());
  p.Zero();
  p[0] = K.intercept;
  for (Int_t a = 0; a < K.n_amp; a++) {
    if (K.form) {
      for (Int_t b = 0; b < lay.n_free; b++)
        p[lay.Low(a, b)] = K.k_low[a][b];
      if (K.j_tau >= 0) {
        p[lay.D(K.j_tau, a)] = -K.c[a];
        p[lay.B(K.j_tau, a)] = K.cb[a];
      }
    } else {
      for (Int_t b = 0; b < kNBins; b++)
        p[lay.Bin(a, b)] = K.k_binned[a][b];
    }
  }
  // Only the applied model's features are nonzero in p, and only their
  // mutual blocks were accumulated, so sum over those.
  std::vector<Int_t> nz;
  for (Int_t i = 0; i < lay.Size(); i++)
    if (p[i] != 0.0)
      nz.push_back(i);
  Double_t ss_res = N.syy;
  for (size_t u = 0; u < nz.size(); u++) {
    ss_res -= 2.0 * p[nz[u]] * N.b[nz[u]];
    for (size_t v = 0; v < nz.size(); v++)
      ss_res += p[nz[u]] * p[nz[v]] * N.Get(nz[u], nz[v]);
  }
  return 1.0 - ss_res / ss_tot;
}

Bool_t FitKernel(const Normal &N, const FeatLayout &lay, const FormGrid &grid,
                 Int_t given_j, Bool_t want_form, Kernel &K, TString &why) {
  why = "";
  K.n = N.n;
  K.n_amp = lay.n_amp;
  K.n_free = lay.n_free;
  K.t_m_us = grid.t_m_us;
  K.t_f_us = grid.t_f_us;
  K.t_lo_us = grid.t_lo_us;
  K.ok = kFALSE;
  K.j_tau = -1;
  if (N.n < 100) {
    why = Form("only %lld pairs", N.n);
    return kFALSE;
  }
  const Double_t ybar = N.sy / N.n;
  const Double_t ss_tot = N.syy - N.n * ybar * ybar;
  if (!(ss_tot > 0.0)) {
    why = "no variance";
    return kFALSE;
  }
  const Int_t n_amp = lay.n_amp;
  K.rms_before = TMath::Sqrt(ss_tot / N.n);

  // The binned kernel: intercept plus every (band, bin) somebody populated.
  Bool_t binned_ok = kFALSE;
  Double_t intercept_binned = 0.0;
  {
    std::vector<Int_t> keep(1, 0);
    for (Int_t a = 0; a < n_amp; a++)
      for (Int_t b = 0; b < kNBins; b++)
        if (N.Get(lay.Bin(a, b), lay.Bin(a, b)) > 0.0)
          keep.push_back(lay.Bin(a, b));
    TVectorD p;
    Double_t drop = 0.0;
    if (SolveSubset(N, keep, p, drop)) {
      const Double_t ss_res = N.syy - drop;
      for (Int_t i = 0; i < Int_t(keep.size()); i++) {
        if (keep[i] == 0)
          intercept_binned = p[i];
        else
          K.k_binned[(keep[i] - 1) / kNBins][(keep[i] - 1) % kNBins] = p[i];
      }
      K.r2_binned = 1.0 - ss_res / ss_tot;
      K.rms_after_binned = TMath::Sqrt(TMath::Max(0.0, ss_res) / N.n);
      binned_ok = K.r2_binned > 0.0;
    }
  }

  // The form: for each tau of the grid, intercept plus the free low bins and
  // the two form features per band. The profile keeps the tau with the
  // smallest residual among the grid's own points; a given tau has its own
  // point and its solution is the one applied. c is minus the direct-tail
  // coefficient (undershoot > 0).
  K.form_ok = kFALSE;
  if (lay.HasForm()) {
    Double_t best_res = 0.0, given_res = 0.0;
    Int_t free_j = -1;
    TVectorD best_p, given_p;
    std::vector<Int_t> best_keep, given_keep;
    for (Int_t j = 0; j < lay.n_tau; j++) {
      if (j >= kNTauGrid && j != given_j)
        continue; // another group's given tau
      std::vector<Int_t> keep(1, 0);
      for (Int_t a = 0; a < n_amp; a++) {
        for (Int_t b = 0; b < lay.n_free; b++)
          if (N.Get(lay.Low(a, b), lay.Low(a, b)) > 0.0)
            keep.push_back(lay.Low(a, b));
        if (N.Get(lay.D(j, a), lay.D(j, a)) > 0.0)
          keep.push_back(lay.D(j, a));
        if (N.Get(lay.B(j, a), lay.B(j, a)) > 0.0)
          keep.push_back(lay.B(j, a));
      }
      TVectorD p;
      Double_t drop = 0.0;
      if (!SolveSubset(N, keep, p, drop))
        continue;
      const Double_t ss_res = N.syy - drop;
      if (j == given_j) {
        given_res = ss_res;
        given_p.ResizeTo(p.GetNrows());
        given_p = p;
        given_keep = keep;
        continue;
      }
      if (free_j < 0 || ss_res < best_res) {
        best_res = ss_res;
        free_j = j;
        best_p.ResizeTo(p.GetNrows());
        best_p = p;
        best_keep = keep;
      }
    }
    if (free_j >= 0) {
      K.tau_free_us = grid.tau_us[free_j];
      K.r2_form_free = 1.0 - best_res / ss_tot;
    }
    // Applied: the given tau's solution when there is one, else the
    // profile's; the free result stays alongside as the check.
    if (given_j >= 0 && given_keep.empty())
      why = Form("form at the given tau %.1f us did not solve; profiled tau "
                 "used",
                 grid.tau_us[given_j]);
    const Bool_t use_given = given_j >= 0 && !given_keep.empty();
    const Int_t j = use_given ? given_j : free_j;
    if (use_given) {
      best_res = given_res;
      best_p.ResizeTo(given_p.GetNrows());
      best_p = given_p;
      best_keep = given_keep;
    }
    if (j >= 0) {
      K.form_ok = kTRUE;
      K.j_tau = j;
      K.tau_us = grid.tau_us[j];
      K.tau_given = use_given;
      if (free_j < 0) {
        K.tau_free_us = K.tau_us;
        K.r2_form_free = 1.0 - best_res / ss_tot;
      }
      for (Int_t i = 0; i < Int_t(best_keep.size()); i++) {
        const Int_t f = best_keep[i];
        if (f == 0) {
          K.intercept_form = best_p[i];
          continue;
        }
        Bool_t placed = kFALSE;
        for (Int_t a = 0; a < n_amp && !placed; a++) {
          for (Int_t b = 0; b < lay.n_free && !placed; b++)
            if (f == lay.Low(a, b)) {
              K.k_low[a][b] = best_p[i];
              placed = kTRUE;
            }
          if (!placed && f == lay.D(j, a)) {
            K.c[a] = -best_p[i];
            placed = kTRUE;
          } else if (!placed && f == lay.B(j, a)) {
            K.cb[a] = best_p[i];
            placed = kTRUE;
          }
        }
      }
      K.r2_form = 1.0 - best_res / ss_tot;
      K.rms_after_form = TMath::Sqrt(TMath::Max(0.0, best_res) / N.n);
      if (!(K.r2_form > 0.0))
        K.form_ok = kFALSE;
    }
  }

  // Which one is applied: the form where asked for and fitted, the bins
  // otherwise (a form that failed falls back). K.k holds the applied kernel
  // at the bin centres.
  K.form = want_form && K.form_ok;
  if (K.form) {
    K.intercept = K.intercept_form;
    K.r2 = K.r2_form;
    K.rms_after = K.rms_after_form;
    for (Int_t a = 0; a < n_amp; a++)
      for (Int_t b = 0; b < kNBins; b++)
        K.k[a][b] = K.FormAt(a, BinCentreUs(b));
    K.ok = kTRUE;
  } else {
    K.intercept = intercept_binned;
    K.r2 = K.r2_binned;
    K.rms_after = K.rms_after_binned;
    for (Int_t a = 0; a < n_amp; a++)
      for (Int_t b = 0; b < kNBins; b++)
        K.k[a][b] = K.k_binned[a][b];
    K.ok = binned_ok;
    if (!K.ok && why.IsNull())
      why = "binned fit did not solve";
  }
  return K.ok;
}

} // namespace

Bool_t Measure(std::vector<RawHit> &hits, const std::vector<Int_t> &group_of,
               Result &res, const TString &file_label) {
  res.n_hits = Long64_t(hits.size());
  if (hits.empty())
    return kFALSE;
  Bool_t sorted = kTRUE;
  for (size_t j = 1; j < hits.size() && sorted; j++)
    sorted = hits[j].timestamp >= hits[j - 1].timestamp;
  res.sorted_input = sorted;
  if (!sorted)
    std::stable_sort(hits.begin(), hits.end(),
                     [](const RawHit &a, const RawHit &b) {
                       return a.timestamp < b.timestamp;
                     });

  const Int_t ref = IndexOfName(Constants::ActiveReferenceChannel());
  // Strip 0 must have fired too, else it is not a beam particle through the
  // whole chamber (as in the downstream beam-only selection).
  const Int_t strip0 = IndexOfName("Strip0");
  if (ref < 0) {
    std::cerr << "  " << file_label
              << ": pulse history: no reference channel in the map"
              << std::endl;
    return kFALSE;
  }
  const Double_t ref_lo = Constants::ActiveReferenceChannelMinAdc();
  const Double_t ref_hi = Constants::ActiveReferenceChannelMaxAdc();
  const Double_t window_ps = Constants::ActiveEventTimeWindowUs() * 1.0e6;
  const Int_t nidx = Int_t(group_of.size());
  const Int_t n_amp = TMath::Max(
      1, TMath::Min(kMaxAmpBins, Constants::cfg.PULSE_HISTORY_AMP_BINS));
  // Every channel with a kernel, and the long ends alone, which are what the
  // beam-like selection is made of.
  std::vector<Int_t> fit_ch, long_ch;
  for (Int_t i = 0; i < nidx; i++) {
    if (group_of[i] == kNone)
      continue;
    fit_ch.push_back(i);
    if (IsLongGroup(group_of[i]))
      long_ch.push_back(i);
  }
  if (long_ch.empty())
    return kFALSE;

  // Beam peak per channel: mode of the raw spectrum above the noise pile.
  // Beam-like: every long end in [lo, hi] x mode; short: rough amp bands.
  const Double_t emax = Constants::ActiveStripEMaxAdc();
  const Int_t nb = 128;
  std::vector<std::vector<Long64_t>> spec(nidx, std::vector<Long64_t>(nb, 0));
  for (size_t j = 0; j < hits.size(); j++) {
    const RawHit &h = hits[j];
    const Int_t i = HitIndex(h);
    if (i < 0 || i >= nidx || group_of[i] == kNone)
      continue;
    const Int_t b = Int_t(Double_t(h.energy) / emax * nb);
    if (b >= 0 && b < nb)
      spec[i][b]++;
  }
  res.mode.assign(nidx, 0.0);
  for (Int_t k = 0; k < Int_t(fit_ch.size()); k++) {
    const Int_t i = fit_ch[k];
    Int_t best = -1;
    for (Int_t b = nb / 16; b < nb; b++) // skip the lowest 1/16 of the range
      if (best < 0 || spec[i][b] > spec[i][best])
        best = b;
    res.mode[i] = best >= 0 ? (best + 0.5) * emax / nb : 0.0;
  }
  const std::vector<Double_t> &mode = res.mode;
  const Double_t blo = Constants::cfg.PULSE_HISTORY_BEAM_LO;
  const Double_t bhi = Constants::cfg.PULSE_HISTORY_BEAM_HI;

  // Pass A: seeds making a beam-like event for each group, and the channel
  // means; short ends and strips 0/17 enter the fit only where they fired.
  std::vector<size_t> seeds[kNGroups];
  std::vector<Double_t> mean(nidx, 0.0);
  std::vector<Double_t> ev_e(nidx, 0.0);
  const Double_t olo = Constants::cfg.PULSE_HISTORY_OWN_LO;
  const Double_t ohi = Constants::cfg.PULSE_HISTORY_OWN_HI;
  // The form is fitted for every group as soon as one asks for it, so the
  // report can compare it with the bins everywhere; nothing is fitted for it
  // otherwise, which keeps the scan at its old cost.
  Bool_t want_form = kFALSE;
  for (Int_t g = 1; g < kNGroups; g++)
    want_form = want_form || (GroupEnabled(g) && GroupWantsForm(g));
  FormGrid grid;
  {
    const Double_t rise = Constants::cfg.PULSE_HISTORY_TRAP_RISE_US;
    const Double_t flat = Constants::cfg.PULSE_HISTORY_TRAP_FLAT_US;
    grid.t_m_us = rise + Constants::cfg.PULSE_HISTORY_PEAKING_US;
    grid.t_f_us = 2.0 * rise + flat;
    grid.t_lo_us = grid.t_f_us - grid.t_m_us;
    if (want_form)
      for (Int_t j = 0; j < kNTauGrid; j++) {
        const Double_t tau =
            kTauGridLoUs * TMath::Power(kTauGridHiUs / kTauGridLoUs,
                                        Double_t(j) / Double_t(kNTauGrid - 1));
        grid.tau_us.push_back(tau);
        grid.em.push_back(TMath::Exp(-grid.t_m_us / tau));
        grid.ef.push_back(TMath::Exp(-grid.t_f_us / tau));
      }
  }
  // A tau given for a group is appended to the grid so the applied form is
  // solved at exactly that value; the profile stays over the first kNTauGrid
  // points. given_j[g] is its index, -1 when the group's tau is profiled.
  Int_t given_j[kNGroups];
  for (Int_t g = 0; g < kNGroups; g++)
    given_j[g] = -1;
  if (want_form)
    for (Int_t g = 1; g < kNGroups; g++) {
      const Double_t tau = GroupTauUs(g);
      if (!(tau > 0.0) || !GroupEnabled(g))
        continue;
      for (Int_t j = kNTauGrid; j < Int_t(grid.tau_us.size()); j++)
        if (TMath::Abs(grid.tau_us[j] - tau) < 1.0e-6)
          given_j[g] = j;
      if (given_j[g] < 0) {
        given_j[g] = Int_t(grid.tau_us.size());
        grid.tau_us.push_back(tau);
        grid.em.push_back(TMath::Exp(-grid.t_m_us / tau));
        grid.ef.push_back(TMath::Exp(-grid.t_f_us / tau));
      }
    }
  grid.Tabulate();
  FeatLayout lay;
  lay.n_amp = n_amp;
  lay.n_tau = Int_t(grid.tau_us.size());
  // Bins that start below t_lo stay free under the form; the last of them is
  // only partly below and gets the part that is.
  lay.n_free = want_form ? TMath::Max(0, BinOf(grid.t_lo_us * 1.0e-6) + 1) : 0;
  const Int_t np = lay.Size();
  const Double_t keep_ps = TMath::Power(10.0, kLogHi) * 1.0e12;
  std::vector<Double_t> x(np, 0.0);
  const ScanState st = {hits, group_of, fit_ch, long_ch,   mode, ev_e,
                        mean, nidx,     strip0, window_ps, olo,  ohi,
                        blo,  bhi,      n_amp,  keep_ps,   lay,  grid};
  std::vector<Long64_t> nmean(nidx, 0);
  for (size_t j = 0; j < hits.size(); j++) {
    const RawHit &h = hits[j];
    if (HitIndex(h) != ref)
      continue;
    if (Double_t(h.energy) < ref_lo || Double_t(h.energy) > ref_hi)
      continue;
    res.n_seeds++;
    if (!Lookahead(st, j))
      continue;
    for (Int_t g = 1; g < kNGroups; g++) {
      if (!BeamFor(st, g))
        continue;
      seeds[g].push_back(j);
      for (Int_t k = 0; k < Int_t(fit_ch.size()); k++) {
        const Int_t c = fit_ch[k];
        if (group_of[c] == g && ev_e[c] > 0.0) {
          mean[c] += ev_e[c];
          nmean[c]++;
        }
      }
    }
  }
  res.n_beam_events = 0;
  for (Int_t g = 1; g < kNGroups; g++)
    res.n_beam_events =
        TMath::Max(res.n_beam_events, Long64_t(seeds[g].size()));
  if (res.n_beam_events < Constants::cfg.PULSE_HISTORY_MIN_EVENTS) {
    std::cerr << "  " << file_label << ": pulse history: only "
              << res.n_beam_events << " beam-like events; not measured"
              << std::endl;
    return kFALSE;
  }
  for (Int_t k = 0; k < Int_t(fit_ch.size()); k++) {
    const Int_t c = fit_ch[k];
    mean[c] = nmean[c] ? mean[c] / Double_t(nmean[c]) : 0.0;
  }

  // Pass B: normal equations, one per channel; a group's is their sum. Every
  // channel gets its own kernel where it has the pairs for one, since the
  // preamps differ within a chain; the group kernel is the fallback.
  std::vector<Normal> N(nidx);
  for (Int_t k = 0; k < Int_t(fit_ch.size()); k++)
    N[fit_ch[k]].Init(np);
  // The fit scan builds every feature. Only the blocks a fit reads are
  // accumulated (see Normal), upper triangle: the nonzero features are
  // sorted into their classes first, then each class against itself and the
  // low bins against each tau.
  ScanNeeds all_needs;
  for (Int_t j = 0; j < lay.n_tau; j++)
    all_needs.taus.push_back(j);
  std::vector<Int_t> nz_bin, nz_low;
  std::vector<std::vector<Int_t>> nz_tau(lay.n_tau);
  const Int_t bin_lo = 1, bin_hi = 1 + n_amp * kNBins;
  const Int_t low_lo = bin_hi, low_hi = bin_hi + n_amp * lay.n_free;
  for (Int_t g = 1; g < kNGroups; g++) {
    res.enabled[g] = GroupEnabled(g);
    if (!res.enabled[g])
      continue; // configured off: no kernel, no correction
    Scan(st, g, seeds[g], all_needs, x,
         [&](Int_t c, const std::vector<Double_t> &xx, Double_t y, Double_t) {
           Normal &Nc = N[c];
           Nc.n++;
           Nc.sy += y;
           Nc.syy += y * y;
           Nc.At(0, 0) += 1.0;
           Nc.b[0] += y;
           nz_bin.clear();
           for (Int_t i = bin_lo; i < bin_hi; i++)
             if (xx[i] != 0.0)
               nz_bin.push_back(i);
           nz_low.clear();
           for (Int_t i = low_lo; i < low_hi; i++)
             if (xx[i] != 0.0)
               nz_low.push_back(i);
           for (Int_t j = 0; j < lay.n_tau; j++) {
             nz_tau[j].clear();
             for (Int_t a = 0; a < n_amp; a++) {
               if (xx[lay.D(j, a)] != 0.0)
                 nz_tau[j].push_back(lay.D(j, a));
               if (xx[lay.B(j, a)] != 0.0)
                 nz_tau[j].push_back(lay.B(j, a));
             }
           }
           // A class against the intercept, the right-hand side and itself.
           auto block = [&](const std::vector<Int_t> &v) {
             for (size_t ia = 0; ia < v.size(); ia++) {
               const Int_t i = v[ia];
               const Double_t xi = xx[i];
               Nc.At(0, i) += xi;
               Nc.b[i] += xi * y;
               for (size_t ib = ia; ib < v.size(); ib++)
                 Nc.At(i, v[ib]) += xi * xx[v[ib]];
             }
           };
           block(nz_bin);
           block(nz_low);
           for (Int_t j = 0; j < lay.n_tau; j++) {
             const std::vector<Int_t> &vt = nz_tau[j];
             if (vt.empty())
               continue;
             block(vt);
             for (size_t ia = 0; ia < nz_low.size(); ia++)
               for (size_t ib = 0; ib < vt.size(); ib++)
                 Nc.At(nz_low[ia], vt[ib]) += xx[nz_low[ia]] * xx[vt[ib]];
           }
         });
  }

  res.group_of = group_of;
  res.name_ch.assign(nidx, "");
  for (Int_t k = 0; k < Int_t(fit_ch.size()); k++)
    res.name_ch[fit_ch[k]] = NameOfIndex(fit_ch[k]);
  res.kernel_ch.assign(nidx, Kernel());
  Bool_t any = kFALSE;
  for (Int_t g = 1; g < kNGroups; g++) {
    Kernel &G = res.kernel[g];
    if (!res.enabled[g]) {
      G.ok = kFALSE;
      continue;
    }
    // The group problem is the sum over its channels.
    Normal Ng;
    Ng.Init(np);
    for (Int_t k = 0; k < Int_t(fit_ch.size()); k++)
      if (group_of[fit_ch[k]] == g)
        Ng.Add(N[fit_ch[k]]);
    TString why;
    FitKernel(Ng, lay, grid, given_j[g], GroupWantsForm(g), G, why);
    G.n_beam = Long64_t(seeds[g].size());
    if (!why.IsNull())
      std::cerr << "  " << file_label << ": pulse history " << GroupName(g)
                << ": " << why << std::endl;
    // Each channel: its own kernel when it has the pairs for one and the fit
    // succeeded, else a copy of the group's, flagged.
    for (Int_t k = 0; k < Int_t(fit_ch.size()); k++) {
      const Int_t c = fit_ch[k];
      if (group_of[c] != g)
        continue;
      Kernel &K = res.kernel_ch[c];
      TString why_c;
      Bool_t own =
          N[c].n >= kMinPairsPerChannel &&
          FitKernel(N[c], lay, grid, given_j[g], GroupWantsForm(g), K, why_c);
      if (!own) {
        if (why_c.IsNull())
          why_c = Form("only %lld pairs, fewer than %lld", N[c].n,
                       kMinPairsPerChannel);
        const Long64_t n_own = N[c].n;
        K = G;
        K.n = n_own;
        K.from_group = kTRUE;
        K.why_group = why_c;
      } else if (!why_c.IsNull()) {
        std::cerr << "  " << file_label << ": pulse history " << res.name_ch[c]
                  << ": " << why_c << std::endl;
      }
      // What the group's kernel would have done on this channel: the gain
      // from a kernel of its own is the difference to K.r2.
      K.r2_group = G.ok ? R2Of(G, N[c], lay) : 0.0;
      K.n_beam = G.n_beam;
      any = any || K.ok;
    }
  }
  if (!any)
    return kFALSE;

  // Decay summary per fitted group: pre-correction deviation vs previous-pulse
  // time. p2 = the preamp decay the pole-zero missed, hardware-study number.
  for (Int_t g = 1; g < kNGroups; g++) {
    Kernel &K = res.kernel[g];
    if (!K.ok)
      continue;
    std::vector<Double_t> dt_us, dev;
    // Only the deviation and the gap are needed: no features built.
    ScanNeeds none;
    none.bins = kFALSE;
    none.low = kFALSE;
    Scan(st, g, seeds[g], none, x,
         [&](Int_t, const std::vector<Double_t> &, Double_t y,
             Double_t dt_prev) {
           if (dt_prev >= kTauFitLoUs && dt_prev <= kTauFitHiUs) {
             dt_us.push_back(dt_prev);
             dev.push_back(y);
           }
         });
    FitTauDecay(dt_us, dev, K.tau);
  }

  // Pass C: diagnostics, with the applied (per-channel) kernels. Log-time
  // axes are log10(dt) where dt is the time difference in microseconds.
  if (Constants::cfg.SAVE_PLOTS) {
    const TString tag = file_label;
    const Double_t xlo = kLogLo + 6.0, xhi = kLogHi + 6.0;
    for (Int_t g = 1; g < kNGroups; g++) {
      const char *gn = GroupTag(g);
      res.dev_vs_pred[g] =
          new TH2D(Form("h_ph_dev_vs_pred_%s_%s", gn, tag.Data()),
                   ";Predicted Deviation from Mean "
                   "[ADC];#splitline{Measured}{Deviation from Mean [ADC]}",
                   120, -600.0, 600.0, 120, -600.0, 600.0);
      res.dev_before[g] = new TH1D(
          Form("h_ph_dev_before_%s_%s", gn, tag.Data()),
          Form(";%s #minus Mean %s [ADC];Events", gn, gn), 240, -600.0, 600.0);
      res.dev_after[g] = new TH1D(
          Form("h_ph_dev_after_%s_%s", gn, tag.Data()),
          Form(";%s #minus Mean %s [ADC];Events", gn, gn), 240, -600.0, 600.0);
      res.dtprev_before[g] = new TProfile(
          Form("p_ph_dtprev_before_%s_%s", gn, tag.Data()),
          Form(";log_{10}(#Deltat [#mus]);%s #minus Mean %s [ADC]", gn, gn), 48,
          xlo, xhi);
      res.dtprev_after[g] = new TProfile(
          Form("p_ph_dtprev_after_%s_%s", gn, tag.Data()),
          Form(";log_{10}(#Deltat [#mus]);%s #minus Mean %s [ADC]", gn, gn), 48,
          xlo, xhi);
      res.shift[g] = new TH1D(Form("h_ph_shift_%s_%s", gn, tag.Data()),
                              ";Applied Shift [ADC];Hits", 240, -600.0, 600.0);
      TH1 *hists[6] = {static_cast<TH1 *>(res.dev_vs_pred[g]),
                       static_cast<TH1 *>(res.dev_before[g]),
                       static_cast<TH1 *>(res.dev_after[g]),
                       static_cast<TH1 *>(res.dtprev_before[g]),
                       static_cast<TH1 *>(res.dtprev_after[g]),
                       static_cast<TH1 *>(res.shift[g])};
      for (Int_t h = 0; h < 6; h++)
        hists[h]->SetDirectory(nullptr);
    }
    // Per channel: the deviation against dt before and after its own kernel,
    // which is where a channel the group kernel did not serve shows up.
    res.dtprev_before_ch.assign(nidx, nullptr);
    res.dtprev_after_ch.assign(nidx, nullptr);
    for (Int_t k = 0; k < Int_t(fit_ch.size()); k++) {
      const Int_t c = fit_ch[k];
      const char *cn = res.name_ch[c].Data();
      res.dtprev_before_ch[c] = new TProfile(
          Form("p_ph_dtprev_before_ch_%s_%s", cn, tag.Data()),
          Form(";log_{10}(#Deltat [#mus]);%s #minus Mean %s [ADC]", cn, cn), 48,
          xlo, xhi);
      res.dtprev_after_ch[c] = new TProfile(
          Form("p_ph_dtprev_after_ch_%s_%s", cn, tag.Data()),
          Form(";log_{10}(#Deltat [#mus]);%s #minus Mean %s [ADC]", cn, cn), 48,
          xlo, xhi);
      res.dtprev_before_ch[c]->SetDirectory(nullptr);
      res.dtprev_after_ch[c]->SetDirectory(nullptr);
    }
    for (Int_t g = 1; g < kNGroups; g++) {
      if (!res.kernel[g].ok)
        continue;
      // Only the applied models' features: bins where a channel is binned,
      // the low bins and the applied taus where one is on the form.
      ScanNeeds needs;
      needs.bins = kFALSE;
      needs.low = kFALSE;
      for (Int_t k = 0; k < Int_t(fit_ch.size()); k++) {
        const Kernel &K = res.kernel_ch[fit_ch[k]];
        if (group_of[fit_ch[k]] != g || !K.ok)
          continue;
        if (!K.form) {
          needs.bins = kTRUE;
          continue;
        }
        needs.low = kTRUE;
        if (K.j_tau >= 0 && std::find(needs.taus.begin(), needs.taus.end(),
                                      K.j_tau) == needs.taus.end())
          needs.taus.push_back(K.j_tau);
      }
      Scan(st, g, seeds[g], needs, x,
           [&](Int_t c, const std::vector<Double_t> &xx, Double_t y,
               Double_t dt_prev) {
             const Kernel &K = res.kernel_ch[c];
             if (!K.ok)
               return;
             const Double_t pred = Predict(K, lay, xx);
             const Double_t after = y - pred;
             res.dev_vs_pred[g]->Fill(pred, y);
             res.dev_before[g]->Fill(y);
             res.dev_after[g]->Fill(after);
             if (dt_prev > 0.0) {
               const Double_t l = TMath::Log10(dt_prev * 1.0e6);
               res.dtprev_before[g]->Fill(l, y);
               res.dtprev_after[g]->Fill(l, after);
               res.dtprev_before_ch[c]->Fill(l, y);
               res.dtprev_after_ch[c]->Fill(l, after);
             }
           });
    }
  }
  return kTRUE;
}
void Apply(std::vector<RawHit> &hits, const std::vector<Int_t> &group_of,
           Result &res) {
  const Int_t nidx = Int_t(group_of.size());
  const Double_t keep_ps = Constants::cfg.PULSE_HISTORY_APPLY_MAX_US * 1.0e6;
  std::vector<std::deque<Past>> past(nidx);
  Double_t shift_sum[kNGroups] = {0};
  Long64_t shift_n[kNGroups] = {0};
  for (size_t j = 0; j < hits.size(); j++) {
    RawHit &h = hits[j];
    const Int_t i = HitIndex(h);
    if (i < 0 || i >= nidx || group_of[i] == kNone)
      continue;
    const Int_t g = group_of[i];
    const Double_t t = Double_t(h.timestamp);
    const Double_t e_raw = Double_t(h.energy);
    std::deque<Past> &d = past[i];
    while (!d.empty() && d.front().t < t - keep_ps)
      d.pop_front();
    // The channel's own kernel (or the group's copy it was handed).
    const Kernel *K =
        i < Int_t(res.kernel_ch.size()) ? &res.kernel_ch[i] : nullptr;
    if (K && K->ok) {
      const Double_t s =
          Shift(d, t, *K, i < Int_t(res.mode.size()) ? res.mode[i] : 0.0);
      Double_t e = e_raw - s;
      if (e < 0.0) {
        e = 0.0;
        h.flags |= kFlagClamped;
        res.n_clamped++;
        res.n_clamped_group[g]++;
      } else if (e > 65535.0) {
        e = 65535.0;
      }
      h.energy = UShort_t(e + 0.5);
      res.n_corrected++;
      shift_sum[g] += s;
      shift_n[g]++;
      if (res.shift[g])
        res.shift[g]->Fill(-s);
    }
    // The kernel was fitted on raw amplitudes of the previous pulses.
    d.push_back({t, e_raw});
  }
  for (Int_t g = 0; g < kNGroups; g++)
    res.mean_shift[g] = shift_n[g] ? shift_sum[g] / shift_n[g] : 0.0;
}

namespace {

// The form's line of a report entry: which tau, how the free profile compares
// with a given one, and the per-band amplitudes.
TString FormLines(const Kernel &K, const char *indent) {
  TString s;
  s += Form("%sform: tau %.1f us (%s)  R^2 %.3f vs binned %.3f  rms %.1f vs "
            "%.1f ADC%s\n",
            indent, K.tau_us, K.tau_given ? "given" : "profiled", K.r2_form,
            K.r2_binned, K.rms_after_form, K.rms_after_binned,
            (!K.form && K.form_ok) ? "  [bins applied]" : "");
  if (K.tau_given) {
    // The check: where the free profile lands and what it costs to hold the
    // given value. The grid step is about 15 percent, so a free tau within
    // one step of the given one is a match.
    const Double_t ratio = K.tau_us > 0.0 ? K.tau_free_us / K.tau_us : 0.0;
    s += Form("%s  profile prefers tau %.1f us (%.2f x given, grid step "
              "1.15)  R^2 %.4f there vs %.4f at given%s\n",
              indent, K.tau_free_us, ratio, K.r2_form_free, K.r2_form,
              (ratio > 0.0 && ratio < 1.0 / 1.15) || ratio > 1.15
                  ? "  <-- given tau not preferred"
                  : "");
  }
  for (Int_t a = 0; a < K.n_amp; a++)
    s += Form("%s  band %d: c %+.4f (undershoot > 0)  c_B %+.5f  W_eff %.0f "
              "us\n",
              indent, a, K.c[a], K.cb[a], K.EffectiveWindowUs(a));
  return s;
}

const char *StatusOf(const Kernel &K, Bool_t enabled) {
  if (!enabled)
    return "DISABLED";
  if (!K.ok)
    return "NOT USED";
  return K.form ? "FORM    " : "BINNED  ";
}

} // namespace

TString Report(const Result &res, const TString &file_label) {
  TString s;
  s += Form("  pulse history %s: %lld hits%s, %lld reference seeds, %lld "
            "beam-like events; %lld hits corrected, %lld clamped to 0\n",
            file_label.Data(), res.n_hits,
            res.sorted_input ? "" : " (SORTED: input was not time ordered)",
            res.n_seeds, res.n_beam_events, res.n_corrected, res.n_clamped);
  for (Int_t g = 1; g < kNGroups; g++) {
    const Kernel &K = res.kernel[g];
    // The group's own fit (the sum of its channels) is the fallback kernel
    // and the reference the per-channel lines below are read against.
    s += Form("    %-22s %s  beam events %lld  n %lld  R^2 %.3f  rms %.1f -> "
              "%.1f ADC  mean shift %.1f ADC  clamped %lld  [group fit]\n",
              GroupName(g), StatusOf(K, res.enabled[g]), K.n_beam, K.n, K.r2,
              K.rms_before, K.rms_after, res.mean_shift[g],
              res.n_clamped_group[g]);
    if (!res.enabled[g])
      continue;
    if (K.form_ok)
      s += FormLines(K, "      ");
    for (Int_t a = 0; a < K.n_amp; a++) {
      if (K.n_amp == 1)
        s += "      kernel:";
      else
        s += Form("      kernel, previous pulse ~%dx beam:", a);
      for (Int_t b = 0; b < kNBins; b++)
        s += Form(" %.0fus:%+.3f", BinCentreUs(b), K.k[a][b]);
      s += "\n";
    }
    // Decay-time summary of the raw dt profile: p2 = the preamp decay the
    // pole-zero missed; flat profiles (matched chains) report no meaningful p2.
    if (K.tau.ok) {
      if (K.tau.flat)
        s += Form(
            "      tau fit: flat (A %.1f +- %.1f ADC within 2 sigma of 0)\n",
            K.tau.p1, K.tau.p1err);
      else
        s += Form(
            "      tau fit: %.1f +- %.1f us  (A %+.1f +- %.1f ADC,"
            " level %+.1f, baseline lift %+.4f ADC/us, chi2/ndf %.1f/%d)\n",
            K.tau.p2, K.tau.p2err, K.tau.p1, K.tau.p1err, K.tau.p0, K.tau.p3,
            K.tau.chi2, K.tau.ndf);
    }
    // The channels of the group: the kernel each one is corrected with.
    for (Int_t c = 0; c < Int_t(res.kernel_ch.size()); c++) {
      if (res.group_of[c] != g)
        continue;
      const Kernel &C = res.kernel_ch[c];
      s += Form("      %-8s %s  n %-9lld R^2 %.3f  rms %.1f -> %.1f ADC",
                res.name_ch[c].Data(), StatusOf(C, kTRUE), C.n, C.r2,
                C.rms_before, C.rms_after);
      if (C.from_group) {
        s += Form("  [group kernel: %s]\n", C.why_group.Data());
        continue;
      }
      if (C.form_ok) {
        s += Form("  tau %.1f us%s", C.tau_us, C.tau_given ? " given" : "");
        if (C.tau_given)
          s += Form(", profile %.1f", C.tau_free_us);
        s += "  c";
        for (Int_t a = 0; a < C.n_amp; a++)
          s += Form(" %+.3f", C.c[a]);
        s += Form("  (binned R^2 %.3f, group kernel %.3f)", C.r2_binned,
                  C.r2_group);
      } else {
        // No form: the size of the binned kernel past the trapezoid, so a
        // channel that differs from its group can be seen in the log.
        Double_t kmax = 0.0;
        for (Int_t a = 0; a < C.n_amp; a++)
          for (Int_t b = 4; b < kNBins; b++)
            if (TMath::Abs(C.k[a][b]) > TMath::Abs(kmax))
              kmax = C.k[a][b];
        s += Form("  max|k| beyond 4.6 us %+.3f  (group kernel R^2 %.3f)", kmax,
                  C.r2_group);
      }
      s += "\n";
    }
  }
  return s;
}

namespace {

// One kernel figure: the binned coefficients as points, the form as a curve
// over them when it was fitted (solid where applied, dashed where the bins
// are), one colour per amplitude band.
void DrawKernelFigure(const Kernel &K, const TString &subdir,
                      const char *fname) {
  const Int_t colors[kMaxAmpBins] = {kBlack,     kRed + 1,    kAzure + 1,
                                     kGreen + 2, kOrange + 7, kMagenta + 1};
  const Double_t xlo = kLogLo + 6.0, xhi = kLogHi + 6.0;
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  Double_t ylo = 0.0, yhi = 0.0;
  for (Int_t a = 0; a < K.n_amp; a++)
    for (Int_t b = 0; b < kNBins; b++) {
      ylo = TMath::Min(ylo, K.k_binned[a][b]);
      yhi = TMath::Max(yhi, K.k_binned[a][b]);
    }
  TH1F *frame = c->DrawFrame(xlo, ylo - 0.02, xhi, yhi + 0.02);
  frame->SetTitle(";log_{10}(#Deltat [#mus]);"
                  "Relative Amplitude Shift");
  TLegend *leg = nullptr;
  if (K.n_amp > 1 || K.form_ok) {
    // The tail of every kernel runs to zero at the right; the legend goes
    // into whichever right-hand corner the kernel leaves empty, the top for
    // an undershoot (negative kernel) and the bottom for an overshoot.
    if (TMath::Abs(ylo) > yhi)
      leg = PlottingUtils::AddLegend(0.55, 0.89, 0.15, 0.40);
    else
      leg = PlottingUtils::AddLegend(0.55, 0.89, 0.62, 0.87);
  }
  for (Int_t a = 0; a < K.n_amp; a++) {
    TGraph *gr = new TGraph();
    for (Int_t b = 0; b < kNBins; b++)
      gr->SetPoint(b, xlo + (b + 0.5) * (xhi - xlo) / kNBins, K.k_binned[a][b]);
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(colors[a % kMaxAmpBins]);
    gr->SetLineColor(colors[a % kMaxAmpBins]);
    gr->SetLineWidth(2);
    gr->Draw(K.form_ok ? "P SAME" : "PL SAME");
    if (leg)
      leg->AddEntry(gr,
                    K.n_amp > 1 ? Form("Previous pulse ~%d#times beam, bins", a)
                                : "Bins",
                    "p");
    if (K.form_ok) {
      const Int_t nf = 200;
      TGraph *gf = new TGraph();
      for (Int_t i = 0; i < nf; i++) {
        const Double_t l = xlo + (i + 0.5) * (xhi - xlo) / nf;
        gf->SetPoint(i, l, K.FormAt(a, TMath::Power(10.0, l)));
      }
      gf->SetLineColor(colors[a % kMaxAmpBins]);
      gf->SetLineWidth(2);
      gf->SetLineStyle(K.form ? 1 : 2);
      gf->Draw("L SAME");
      leg->AddEntry(gf,
                    Form("%sform, #tau %.0f #mus%s%s",
                         K.n_amp > 1 ? Form("~%d#times beam, ", a) : "",
                         K.tau_us, K.tau_given ? " given" : "",
                         K.form ? " (applied)" : ""),
                    "l");
    }
  }
  if (leg)
    leg->Draw();
  TLine *zero = new TLine(xlo, 0.0, xhi, 0.0);
  zero->SetLineStyle(2);
  zero->Draw();
  PlottingUtils::SaveFigure(c, fname, subdir, PlotSaveOptions::kLINEAR);
  delete c;
}

// A palette for up to 16 channels on one canvas.
Int_t ChannelColor(Int_t k) {
  const Int_t colors[16] = {kBlack,      kRed + 1,     kAzure + 1, kGreen + 2,
                            kOrange + 7, kMagenta + 1, kCyan + 2,  kYellow + 2,
                            kGray + 2,   kRed - 7,     kAzure - 3, kGreen - 6,
                            kOrange - 3, kMagenta - 7, kCyan - 6,  kPink + 7};
  return colors[k % 16];
}

} // namespace

void SavePlots(Result &res, const TString &file_label) {
  if (!Constants::cfg.SAVE_PLOTS)
    return;
  std::lock_guard<std::mutex> lock(g_plot_mutex);
  const TString subdir = "pulse_history/" + file_label;
  const TString subdir_ch = subdir + "/channels";
  const Double_t xlo = kLogLo + 6.0, xhi = kLogHi + 6.0;
  // Kernels: the group fit, one canvas per group.
  for (Int_t g = 1; g < kNGroups; g++)
    if (res.kernel[g].ok)
      DrawKernelFigure(res.kernel[g], subdir, Form("kernel_%s", GroupTag(g)));
  // And every channel's own, under channels/.
  for (Int_t c = 0; c < Int_t(res.kernel_ch.size()); c++)
    if (res.kernel_ch[c].ok && !res.kernel_ch[c].from_group)
      DrawKernelFigure(res.kernel_ch[c], subdir_ch,
                       Form("kernel_%s", res.name_ch[c].Data()));
  // Per group, the applied kernels of its channels on one canvas (beam band),
  // the group's own dashed, so a channel that differs stands out.
  for (Int_t g = 1; g < kNGroups; g++) {
    const Kernel &G = res.kernel[g];
    if (!G.ok)
      continue;
    std::vector<Int_t> chans;
    for (Int_t c = 0; c < Int_t(res.kernel_ch.size()); c++)
      if (res.group_of[c] == g && res.kernel_ch[c].ok)
        chans.push_back(c);
    if (chans.empty())
      continue;
    const Int_t band = G.n_amp > 1 ? 1 : 0;
    Double_t ylo = 0.0, yhi = 0.0;
    for (size_t k = 0; k < chans.size(); k++)
      for (Int_t b = 0; b < kNBins; b++) {
        ylo = TMath::Min(ylo, res.kernel_ch[chans[k]].k[band][b]);
        yhi = TMath::Max(yhi, res.kernel_ch[chans[k]].k[band][b]);
      }
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    TH1F *frame = c->DrawFrame(xlo, ylo - 0.02, xhi, yhi + 0.02);
    frame->SetTitle(Form(";log_{10}(#Deltat [#mus]);Relative Amplitude "
                         "Shift%s",
                         G.n_amp > 1 ? " (previous pulse ~1#times beam)" : ""));
    TLegend *leg = TMath::Abs(ylo) > yhi
                       ? PlottingUtils::AddLegend(0.62, 0.89, 0.15, 0.50)
                       : PlottingUtils::AddLegend(0.62, 0.89, 0.52, 0.87);
    for (size_t k = 0; k < chans.size(); k++) {
      const Kernel &K = res.kernel_ch[chans[k]];
      TGraph *gr = new TGraph();
      for (Int_t b = 0; b < kNBins; b++)
        gr->SetPoint(b, xlo + (b + 0.5) * (xhi - xlo) / kNBins, K.k[band][b]);
      gr->SetLineColor(ChannelColor(Int_t(k)));
      gr->SetMarkerColor(ChannelColor(Int_t(k)));
      gr->SetMarkerStyle(20);
      gr->SetMarkerSize(0.7);
      gr->SetLineWidth(2);
      gr->SetLineStyle(K.from_group ? 3 : 1);
      gr->Draw("PL SAME");
      leg->AddEntry(gr,
                    Form("%s%s", res.name_ch[chans[k]].Data(),
                         K.from_group ? " (group)" : ""),
                    "l");
    }
    TGraph *gg = new TGraph();
    for (Int_t b = 0; b < kNBins; b++)
      gg->SetPoint(b, xlo + (b + 0.5) * (xhi - xlo) / kNBins, G.k[band][b]);
    gg->SetLineColor(kGray + 1);
    gg->SetLineWidth(3);
    gg->SetLineStyle(2);
    gg->Draw("L SAME");
    leg->AddEntry(gg, "Group fit", "l");
    leg->Draw();
    TLine *zero = new TLine(xlo, 0.0, xhi, 0.0);
    zero->SetLineStyle(2);
    zero->Draw();
    PlottingUtils::SaveFigure(c, Form("kernel_channels_%s", GroupTag(g)),
                              subdir, PlotSaveOptions::kLINEAR);
    delete c;
  }
  for (Int_t g = 1; g < kNGroups; g++) {
    const char *gn = GroupTag(g);
    if (res.dev_vs_pred[g] && res.dev_vs_pred[g]->GetEntries() > 0) {
      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
      c->SetLogz(kTRUE);
      PlottingUtils::ConfigureAndDraw2DHistogram(res.dev_vs_pred[g], c);
      TLine *diag = new TLine(-600.0, -600.0, 600.0, 600.0);
      diag->SetLineColor(kRed + 1);
      diag->SetLineStyle(2);
      diag->Draw();
      PlottingUtils::SaveFigure(c, Form("deviation_vs_predicted_%s", gn),
                                subdir, PlotSaveOptions::kLINEAR);
      delete c;
    }
    if (res.dev_before[g] && res.dev_before[g]->GetEntries() > 0) {
      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kTRUE);
      PlottingUtils::ConfigureAndDrawHistogram(res.dev_before[g], kBlack);
      res.dev_after[g]->SetLineColor(kRed + 1);
      res.dev_after[g]->SetLineWidth(2);
      res.dev_after[g]->Draw("HIST SAME");
      TLegend *leg = PlottingUtils::AddLegend(0.62, 0.89, 0.72, 0.88);
      leg->AddEntry(res.dev_before[g],
                    Form("Before, RMS %.1f", res.kernel[g].rms_before), "l");
      leg->AddEntry(res.dev_after[g],
                    Form("After, RMS %.1f", res.kernel[g].rms_after), "l");
      leg->Draw();
      PlottingUtils::SaveFigure(c, Form("deviation_before_after_%s", gn),
                                subdir, PlotSaveOptions::kLOG);
      delete c;
    }
    if (res.dtprev_before[g] && res.dtprev_before[g]->GetEntries() > 0) {
      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
      res.dtprev_before[g]->SetMinimum(-300.0);
      res.dtprev_before[g]->SetMaximum(200.0);
      PlottingUtils::ConfigureAndDrawHistogram(res.dtprev_before[g], kBlack);
      res.dtprev_after[g]->SetLineColor(kRed + 1);
      res.dtprev_after[g]->SetMarkerColor(kRed + 1);
      res.dtprev_after[g]->SetLineWidth(2);
      res.dtprev_after[g]->Draw("SAME");
      TLegend *leg = PlottingUtils::AddLegend(0.62, 0.89, 0.15, 0.30);
      leg->AddEntry(res.dtprev_before[g], "Before", "l");
      leg->AddEntry(res.dtprev_after[g], "After", "l");
      leg->Draw();
      PlottingUtils::SaveFigure(c, Form("deviation_vs_dt_previous_%s", gn),
                                subdir, PlotSaveOptions::kLINEAR);
      delete c;
    }
    // The channels of the group after their own kernels, on one canvas:
    // a channel its kernel did not serve keeps a dt dependence here.
    for (Int_t pass = 0; pass < 2; pass++) {
      std::vector<TProfile *> &prof =
          pass == 0 ? res.dtprev_before_ch : res.dtprev_after_ch;
      std::vector<Int_t> chans;
      for (Int_t c = 0; c < Int_t(prof.size()); c++)
        if (res.group_of[c] == g && prof[c] && prof[c]->GetEntries() > 0)
          chans.push_back(c);
      if (chans.empty())
        continue;
      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
      TProfile *first = prof[chans[0]];
      first->SetMinimum(pass == 0 ? -300.0 : -100.0);
      first->SetMaximum(pass == 0 ? 200.0 : 100.0);
      first->SetTitle(Form(";log_{10}(#Deltat [#mus]);Channel #minus Mean "
                           "[ADC], %s",
                           pass == 0 ? "before" : "after"));
      PlottingUtils::ConfigureAndDrawHistogram(first, ChannelColor(0));
      TLegend *leg = PlottingUtils::AddLegend(0.62, 0.89, 0.15, 0.50);
      leg->AddEntry(first, res.name_ch[chans[0]].Data(), "l");
      for (size_t k = 1; k < chans.size(); k++) {
        TProfile *p = prof[chans[k]];
        p->SetLineColor(ChannelColor(Int_t(k)));
        p->SetMarkerColor(ChannelColor(Int_t(k)));
        p->SetLineWidth(2);
        p->Draw("SAME");
        leg->AddEntry(p, res.name_ch[chans[k]].Data(), "l");
      }
      leg->Draw();
      PlottingUtils::SaveFigure(c,
                                Form("deviation_vs_dt_channels_%s_%s",
                                     pass == 0 ? "before" : "after", gn),
                                subdir, PlotSaveOptions::kLINEAR);
      delete c;
    }
    if (res.shift[g] && res.shift[g]->GetEntries() > 0) {
      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kTRUE);
      PlottingUtils::ConfigureAndDrawHistogram(res.shift[g], kBlack);
      PlottingUtils::SaveFigure(c, Form("applied_shift_%s", gn), subdir,
                                PlotSaveOptions::kLOG);
      delete c;
    }
  }
  res.FreeDiagnostics();
}

void WriteToEventsFile(const TString &events_subpath, const Result &res) {
  TFile *f = IO::OpenForWriting(events_subpath, "UPDATE");
  if (!f || f->IsZombie()) {
    if (f)
      delete f;
    return;
  }
  f->cd();
  if (TObject *old = f->Get("pulse_history"))
    old->Delete();
  // One entry per channel with a kernel (Board/Channel from the map, Name
  // the channel's), plus one per group (Board and Channel -1, Name the
  // group's tag): the group fit that is the channels' fallback.
  TTree *t = new TTree("pulse_history",
                       "Pulse-history kernel per channel and per group");
  Int_t group = 0, board = -1, channel = -1, n_amp = 1;
  Char_t name[32] = {0};
  Bool_t ok = kFALSE, enabled = kTRUE, from_group = kFALSE, form = kFALSE,
         form_ok = kFALSE;
  Bool_t form_tau_given = kFALSE;
  Double_t form_tau_us = 0.0, form_tau_free_us = 0.0, r2_form_free = 0.0,
           form_c[kMaxAmpBins], form_cb[kMaxAmpBins], r2_form = 0.0,
           r2_binned = 0.0, r2_group = 0.0, k_binned[kMaxAmpBins * kNBins];
  Double_t tau_us = 0.0, tau_err_us = 0.0, tau_p0 = 0.0, tau_p1 = 0.0,
           tau_p3 = 0.0;
  Bool_t tau_ok = kFALSE, tau_flat = kTRUE;
  Double_t k[kMaxAmpBins * kNBins], centre_us[kNBins],
      intercept = 0.0, r2 = 0.0, rms_before = 0.0, rms_after = 0.0,
      mean_shift = 0.0, apply_max_us = 0.0;
  Long64_t n = 0, n_beam = 0, n_corrected = 0, n_clamped = 0;
  t->Branch("Group", &group, "Group/I");
  t->Branch("Board", &board, "Board/I");
  t->Branch("Channel", &channel, "Channel/I");
  t->Branch("Name", name, "Name/C");
  t->Branch("Ok", &ok, "Ok/O");
  // Whether PULSE_HISTORY_GROUPS had the group on; off means the group was
  // never fitted, as opposed to fitted and rejected (Ok false, Enabled true).
  t->Branch("Enabled", &enabled, "Enabled/O");
  // Channel rows only: the channel is corrected with the group's kernel
  // because its own could not be fitted.
  t->Branch("FromGroup", &from_group, "FromGroup/O");
  t->Branch("NAmpBins", &n_amp, "NAmpBins/I");
  // Kernel is what was applied, at the bin centres (the form evaluated there
  // when Form is true); KernelBinned the binned fit, always.
  t->Branch("Kernel", k, Form("Kernel[%d]/D", kMaxAmpBins * kNBins));
  t->Branch("KernelBinned", k_binned,
            Form("KernelBinned[%d]/D", kMaxAmpBins * kNBins));
  t->Branch("Form", &form, "Form/O");
  t->Branch("FormOk", &form_ok, "FormOk/O");
  t->Branch("FormTauUs", &form_tau_us, "FormTauUs/D");
  t->Branch("FormTauGiven", &form_tau_given, "FormTauGiven/O");
  t->Branch("FormTauFreeUs", &form_tau_free_us, "FormTauFreeUs/D");
  t->Branch("R2FormFree", &r2_form_free, "R2FormFree/D");
  t->Branch("FormC", form_c, Form("FormC[%d]/D", kMaxAmpBins));
  t->Branch("FormCB", form_cb, Form("FormCB[%d]/D", kMaxAmpBins));
  t->Branch("R2Form", &r2_form, "R2Form/D");
  t->Branch("R2Binned", &r2_binned, "R2Binned/D");
  t->Branch("R2Group", &r2_group, "R2Group/D");
  t->Branch("BinCentreUs", centre_us, Form("BinCentreUs[%d]/D", kNBins));
  t->Branch("Intercept", &intercept, "Intercept/D");
  t->Branch("R2", &r2, "R2/D");
  t->Branch("RmsBefore", &rms_before, "RmsBefore/D");
  t->Branch("RmsAfter", &rms_after, "RmsAfter/D");
  t->Branch("MeanShift", &mean_shift, "MeanShift/D");
  t->Branch("ApplyMaxUs", &apply_max_us, "ApplyMaxUs/D");
  // Decay-time summary of the raw dt profile (p0 + p1 exp(-dt/p2) + p3 dt);
  // see Kernel::TauFit. Group rows only; flat groups have tau_ok but
  // tau_flat and p2 = 0.
  t->Branch("TauOk", &tau_ok, "TauOk/O");
  t->Branch("TauFlat", &tau_flat, "TauFlat/O");
  t->Branch("TauUs", &tau_us, "TauUs/D");
  t->Branch("TauErrUs", &tau_err_us, "TauErrUs/D");
  t->Branch("TauP0", &tau_p0, "TauP0/D");
  t->Branch("TauAmp", &tau_p1, "TauAmp/D");
  t->Branch("TauLift", &tau_p3, "TauLift/D");
  t->Branch("N", &n, "N/L");
  t->Branch("NBeamEvents", &n_beam, "NBeamEvents/L");
  t->Branch("NCorrected", &n_corrected, "NCorrected/L");
  t->Branch("NClamped", &n_clamped, "NClamped/L");
  const Int_t nch = Constants::ActiveNChannels();
  auto fill = [&](const Kernel &K, Int_t g, Int_t idx) {
    group = g;
    board = idx >= 0 ? idx / nch : -1;
    channel = idx >= 0 ? idx % nch : -1;
    strncpy(name, idx >= 0 ? res.name_ch[idx].Data() : GroupTag(g),
            sizeof(name) - 1);
    name[sizeof(name) - 1] = 0;
    ok = K.ok;
    enabled = res.enabled[g];
    from_group = K.from_group;
    form = K.form;
    form_ok = K.form_ok;
    form_tau_us = K.tau_us;
    form_tau_given = K.tau_given;
    form_tau_free_us = K.tau_free_us;
    r2_form_free = K.r2_form_free;
    r2_form = K.r2_form;
    r2_binned = K.r2_binned;
    r2_group = K.r2_group;
    n_amp = K.n_amp;
    for (Int_t a = 0; a < kMaxAmpBins; a++) {
      form_c[a] = K.c[a];
      form_cb[a] = K.cb[a];
      for (Int_t b = 0; b < kNBins; b++) {
        k[a * kNBins + b] = K.k[a][b];
        k_binned[a * kNBins + b] = K.k_binned[a][b];
      }
    }
    for (Int_t b = 0; b < kNBins; b++)
      centre_us[b] = BinCentreUs(b);
    intercept = K.intercept;
    r2 = K.r2;
    rms_before = K.rms_before;
    rms_after = K.rms_after;
    mean_shift = res.mean_shift[g];
    apply_max_us = Constants::cfg.PULSE_HISTORY_APPLY_MAX_US;
    tau_ok = K.tau.ok;
    tau_flat = K.tau.flat;
    tau_us = K.tau.ok ? K.tau.p2 : 0.0;
    tau_err_us = K.tau.ok ? K.tau.p2err : 0.0;
    tau_p0 = K.tau.p0;
    tau_p1 = K.tau.p1;
    tau_p3 = K.tau.p3;
    n = K.n;
    n_beam = res.n_beam_events;
    n_corrected = res.n_corrected;
    n_clamped = res.n_clamped;
    t->Fill();
  };
  for (Int_t g = 1; g < kNGroups; g++)
    fill(res.kernel[g], g, -1);
  for (Int_t c = 0; c < Int_t(res.kernel_ch.size()); c++)
    if (res.group_of[c] != kNone)
      fill(res.kernel_ch[c], res.group_of[c], c);
  t->Write("pulse_history", TObject::kOverwrite);
  f->Close();
  delete f;
}

} // namespace PulseHistory
