#include "PulseHistory.hpp"
#include "Constants.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TDecompSVD.h>
#include <TFile.h>
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
#include <deque>
#include <functional>
#include <iostream>

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
  case kGuard0:
    return "guard strip 0";
  case kGuard17:
    return "guard strip 17";
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
  case kGuard0:
    return "S0";
  case kGuard17:
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

Kernel::Kernel() {
  for (Int_t a = 0; a < kMaxAmpBins; a++)
    for (Int_t b = 0; b < kNBins; b++)
      k[a][b] = 0.0;
}

Result::Result() {
  for (Int_t g = 0; g < kNGroups; g++) {
    mean_shift[g] = 0.0;
    n_clamped_group[g] = 0;
    dev_vs_pred[g] = nullptr;
    dev_before[g] = dev_after[g] = shift[g] = nullptr;
    dtprev_before[g] = dtprev_after[g] = nullptr;
  }
}

Result::~Result() {
  for (Int_t g = 0; g < kNGroups; g++) {
    delete dev_vs_pred[g];
    delete dev_before[g];
    delete dev_after[g];
    delete shift[g];
    delete dtprev_before[g];
    delete dtprev_after[g];
  }
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
      gm[idx] = kGuard0;
      continue;
    }
    if (name == "Strip17") {
      gm[idx] = kGuard17;
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

// Feature index of (amplitude band a, dt bin b); 0 is the intercept.
inline Int_t Feat(Int_t a, Int_t b) { return 1 + a * kNBins + b; }

// Predicted shift from a channel's past, with the group's kernel.
Double_t Shift(const std::deque<Past> &d, Double_t t, const Kernel &K,
               Double_t mode) {
  Double_t s = 0.0;
  for (const Past &p : d) {
    const Int_t b = BinOf((t - p.t) * 1.0e-12);
    if (b < 0)
      continue;
    s += K.k[AmpBinOf(p.e, mode, K.n_amp)][b] * p.e;
  }
  return s;
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
  // The entrance guard has to have fired too, or the event is not a beam
  // particle through the whole chamber (same requirement as the beam-only
  // selection downstream).
  const Int_t guard0 = IndexOfName("Strip0");
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

  // Beam peak per channel: the mode of its raw spectrum above the noise pile.
  // Beam-like means every long end inside [lo, hi] x its mode. (For a short
  // end most of the spectrum sits under that pile, so its mode is only a
  // rough scale, used for the amplitude bands and nothing else.)
  const Double_t emax = Constants::ActiveStripEMaxAdc();
  const Int_t nb = 128;
  std::vector<std::vector<Long64_t>> spec(nidx, std::vector<Long64_t>(nb, 0));
  for (const RawHit &h : hits) {
    const Int_t i = HitIndex(h);
    if (i < 0 || i >= nidx || group_of[i] == kNone)
      continue;
    const Int_t b = Int_t(Double_t(h.energy) / emax * nb);
    if (b >= 0 && b < nb)
      spec[i][b]++;
  }
  res.mode.assign(nidx, 0.0);
  for (Int_t i : fit_ch) {
    Int_t best = -1;
    for (Int_t b = nb / 16; b < nb; b++) // skip the lowest 1/16 of the range
      if (best < 0 || spec[i][b] > spec[i][best])
        best = b;
    res.mode[i] = best >= 0 ? (best + 0.5) * emax / nb : 0.0;
  }
  const std::vector<Double_t> &mode = res.mode;
  const Double_t blo = Constants::cfg.PULSE_HISTORY_BEAM_LO;
  const Double_t bhi = Constants::cfg.PULSE_HISTORY_BEAM_HI;

  // Pass A: seeds that make a beam-like event for each group, and the channel
  // means. A beam-like event for group g has the guard fired, every long end
  // of the OTHER chain inside the tight window, and every long end of g's own
  // chain inside the loose one: the tight window on the fitted chain would
  // cut off exactly the large undershoots the kernel has to fit. The guards
  // belong to neither chain, so both are held tight for them. A short end or
  // guard enters the fit only in the events where it fired.
  std::vector<size_t> seeds[kNGroups];
  std::vector<Double_t> mean(nidx, 0.0);
  std::vector<Double_t> ev_e(nidx, 0.0);
  const Double_t olo = Constants::cfg.PULSE_HISTORY_OWN_LO;
  const Double_t ohi = Constants::cfg.PULSE_HISTORY_OWN_HI;
  auto lookahead = [&](size_t i0) {
    for (Int_t c : fit_ch)
      ev_e[c] = 0.0;
    Bool_t guard_fired = guard0 < 0;
    const ULong64_t t0 = hits[i0].timestamp;
    for (size_t j = i0 + 1;
         j < hits.size() && hits[j].timestamp - t0 < window_ps; j++) {
      const Int_t i = HitIndex(hits[j]);
      if (i == guard0 && hits[j].energy > 0)
        guard_fired = kTRUE;
      if (i < 0 || i >= nidx || group_of[i] == kNone)
        continue;
      if (Double_t(hits[j].energy) > ev_e[i])
        ev_e[i] = Double_t(hits[j].energy);
    }
    return guard_fired;
  };
  auto beam_for = [&](Int_t g) {
    const Int_t chain = ChainOf(g);
    for (Int_t c : long_ch) {
      if (!(mode[c] > 0.0))
        return kFALSE;
      const Double_t r = ev_e[c] / mode[c];
      const Bool_t own = chain >= 0 && ChainOf(group_of[c]) == chain;
      if (own ? (r < olo || r > ohi) : (r < blo || r > bhi))
        return kFALSE;
    }
    return kTRUE;
  };
  std::vector<Long64_t> nmean(nidx, 0);
  for (size_t j = 0; j < hits.size(); j++) {
    const RawHit &h = hits[j];
    if (HitIndex(h) != ref)
      continue;
    if (Double_t(h.energy) < ref_lo || Double_t(h.energy) > ref_hi)
      continue;
    res.n_seeds++;
    if (!lookahead(j))
      continue;
    for (Int_t g = 1; g < kNGroups; g++) {
      if (!beam_for(g))
        continue;
      seeds[g].push_back(j);
      for (Int_t c : fit_ch)
        if (group_of[c] == g && ev_e[c] > 0.0) {
          mean[c] += ev_e[c];
          nmean[c]++;
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
  for (Int_t c : fit_ch)
    mean[c] = nmean[c] ? mean[c] / Double_t(nmean[c]) : 0.0;

  // One walk over the hits for group g, calling fn(channel, features x,
  // deviation y, dt to the previous pulse) for every channel of g at every
  // beam-like seed of g. The per-channel deque holds that channel's earlier
  // pulses only, so the event's own hit (which comes after the seed) never
  // enters.
  const Int_t np = 1 + n_amp * kNBins;
  const Double_t keep_ps = TMath::Power(10.0, kLogHi) * 1.0e12;
  std::vector<Double_t> x(np, 0.0);
  auto scan = [&](Int_t g,
                  const std::function<void(Int_t, const std::vector<Double_t> &,
                                           Double_t, Double_t)> &fn) {
    std::vector<std::deque<Past>> past(nidx);
    size_t next_seed = 0;
    for (size_t j = 0; j < hits.size(); j++) {
      if (next_seed < seeds[g].size() && seeds[g][next_seed] == j) {
        next_seed++;
        const Double_t tg = Double_t(hits[j].timestamp);
        lookahead(j);
        for (Int_t c : fit_ch) {
          // A short end or guard that did not fire has no height to fit.
          if (group_of[c] != g || !(ev_e[c] > 0.0))
            continue;
          std::deque<Past> &d = past[c];
          while (!d.empty() && d.front().t < tg - keep_ps)
            d.pop_front();
          std::fill(x.begin(), x.end(), 0.0);
          x[0] = 1.0;
          for (const Past &p : d) {
            const Int_t b = BinOf((tg - p.t) * 1.0e-12);
            if (b >= 0)
              x[Feat(AmpBinOf(p.e, mode[c], n_amp), b)] += p.e;
          }
          const Double_t dt_prev =
              d.empty() ? -1.0 : (tg - d.back().t) * 1.0e-12;
          fn(c, x, ev_e[c] - mean[c], dt_prev);
        }
      }
      const Int_t i = HitIndex(hits[j]);
      if (i >= 0 && i < nidx && group_of[i] == g)
        past[i].push_back(
            {Double_t(hits[j].timestamp), Double_t(hits[j].energy)});
    }
  };

  // Pass B: normal equations per group.
  TMatrixD A[kNGroups];
  TVectorD bv[kNGroups];
  Double_t syy[kNGroups] = {0}, sy[kNGroups] = {0};
  Long64_t ng[kNGroups] = {0};
  for (Int_t g = 0; g < kNGroups; g++) {
    A[g].ResizeTo(np, np);
    bv[g].ResizeTo(np);
    A[g].Zero();
    bv[g].Zero();
  }
  for (Int_t g = 1; g < kNGroups; g++)
    scan(g, [&](Int_t, const std::vector<Double_t> &xx, Double_t y, Double_t) {
      for (Int_t a = 0; a < np; a++) {
        if (xx[a] == 0.0)
          continue;
        bv[g][a] += xx[a] * y;
        for (Int_t b = 0; b < np; b++)
          A[g](a, b) += xx[a] * xx[b];
      }
      syy[g] += y * y;
      sy[g] += y;
      ng[g]++;
    });
  Bool_t any = kFALSE;
  for (Int_t g = 1; g < kNGroups; g++) {
    Kernel &K = res.kernel[g];
    K.n = ng[g];
    K.n_amp = n_amp;
    K.n_beam = Long64_t(seeds[g].size());
    if (ng[g] < 100)
      continue;
    // Drop features nobody populated (an amplitude band with no pulses).
    std::vector<Int_t> keep;
    for (Int_t a = 0; a < np; a++)
      if (a == 0 || A[g](a, a) > 0.0)
        keep.push_back(a);
    const Int_t nk = Int_t(keep.size());
    TMatrixD Ak(nk, nk);
    TVectorD bk(nk);
    for (Int_t a = 0; a < nk; a++) {
      bk[a] = bv[g][keep[a]];
      for (Int_t b = 0; b < nk; b++)
        Ak(a, b) = A[g](keep[a], keep[b]);
    }
    // SVD rather than LU: the far-dt bins of the rarer amplitude bands are
    // nearly collinear with the intercept, and a plain inversion gives up on
    // them. Column scaling first, so the tolerance means the same for every
    // feature.
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
      continue;
    TVectorD p(nk);
    for (Int_t a = 0; a < nk; a++)
      p[a] = ps[a] * scale[a];
    Double_t ss_res = syy[g];
    for (Int_t a = 0; a < nk; a++)
      ss_res -= p[a] * bk[a];
    const Double_t ybar = sy[g] / ng[g];
    const Double_t ss_tot = syy[g] - ng[g] * ybar * ybar;
    if (!(ss_tot > 0.0) || !std::isfinite(ss_res))
      continue;
    for (Int_t a = 0; a < nk; a++) {
      if (keep[a] == 0)
        K.intercept = p[a];
      else
        K.k[(keep[a] - 1) / kNBins][(keep[a] - 1) % kNBins] = p[a];
    }
    K.r2 = 1.0 - ss_res / ss_tot;
    K.rms_before = TMath::Sqrt(ss_tot / ng[g]);
    K.rms_after = TMath::Sqrt(TMath::Max(0.0, ss_res) / ng[g]);
    K.ok = K.r2 > 0.0;
    any = any || K.ok;
  }
  if (!any)
    return kFALSE;

  // Pass C: diagnostics, with the fitted kernels. Log-time axes are
  // log10(dt) where dt is the time difference in microseconds
  if (Constants::cfg.SAVE_PLOTS) {
    const TString tag = file_label;
    const Double_t xlo = kLogLo + 6.0, xhi = kLogHi + 6.0;
    for (Int_t g = 1; g < kNGroups; g++) {
      const char *gn = GroupTag(g);
      res.dev_vs_pred[g] = new TH2D(
          Form("h_ph_dev_vs_pred_%s_%s", gn, tag.Data()),
          ";Predicted Deviation from Mean [ADC];#splitline{Measured}{Deviation from Mean [ADC]}",
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
      for (TH1 *h : {static_cast<TH1 *>(res.dev_vs_pred[g]),
                     static_cast<TH1 *>(res.dev_before[g]),
                     static_cast<TH1 *>(res.dev_after[g]),
                     static_cast<TH1 *>(res.dtprev_before[g]),
                     static_cast<TH1 *>(res.dtprev_after[g]),
                     static_cast<TH1 *>(res.shift[g])})
        h->SetDirectory(nullptr);
    }
    for (Int_t g = 1; g < kNGroups; g++) {
      const Kernel &K = res.kernel[g];
      if (!K.ok)
        continue;
      scan(g, [&](Int_t, const std::vector<Double_t> &xx, Double_t y,
                  Double_t dt_prev) {
        Double_t pred = K.intercept;
        for (Int_t a = 0; a < n_amp; a++)
          for (Int_t b = 0; b < kNBins; b++)
            pred += K.k[a][b] * xx[Feat(a, b)];
        const Double_t after = y - pred;
        res.dev_vs_pred[g]->Fill(pred, y);
        res.dev_before[g]->Fill(y);
        res.dev_after[g]->Fill(after);
        if (dt_prev > 0.0) {
          const Double_t l = TMath::Log10(dt_prev * 1.0e6);
          res.dtprev_before[g]->Fill(l, y);
          res.dtprev_after[g]->Fill(l, after);
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
    if (res.kernel[g].ok) {
      const Double_t s = Shift(d, t, res.kernel[g],
                               i < Int_t(res.mode.size()) ? res.mode[i] : 0.0);
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

TString Report(const Result &res, const TString &file_label) {
  TString s;
  s += Form("  pulse history %s: %lld hits%s, %lld reference seeds, %lld "
            "beam-like events; %lld hits corrected, %lld clamped to 0\n",
            file_label.Data(), res.n_hits,
            res.sorted_input ? "" : " (SORTED: input was not time ordered)",
            res.n_seeds, res.n_beam_events, res.n_corrected, res.n_clamped);
  for (Int_t g = 1; g < kNGroups; g++) {
    const Kernel &K = res.kernel[g];
    s += Form("    %-22s %s  beam events %lld  n %lld  R^2 %.3f  rms %.1f -> "
              "%.1f ADC  mean shift %.1f ADC  clamped %lld\n",
              GroupName(g), K.ok ? "fitted " : "NOT USED", K.n_beam, K.n, K.r2,
              K.rms_before, K.rms_after, res.mean_shift[g],
              res.n_clamped_group[g]);
    for (Int_t a = 0; a < K.n_amp; a++) {
      if (K.n_amp == 1)
        s += "      kernel:";
      else
        s += Form("      kernel, previous pulse ~%dx beam:", a);
      for (Int_t b = 0; b < kNBins; b++)
        s += Form(" %.0fus:%+.3f", BinCentreUs(b), K.k[a][b]);
      s += "\n";
    }
  }
  return s;
}

void SavePlots(Result &res, const TString &file_label) {
  if (!Constants::cfg.SAVE_PLOTS)
    return;
  std::lock_guard<std::mutex> lock(g_plot_mutex);
  const TString subdir = "pulse_history/" + file_label;
  const Int_t colors[kMaxAmpBins] = {kBlack,     kRed + 1,    kAzure + 1,
                                     kGreen + 2, kOrange + 7, kMagenta + 1};
  const Double_t xlo = kLogLo + 6.0, xhi = kLogHi + 6.0;
  // Kernels: one canvas per group, one curve per amplitude band.
  for (Int_t g = 1; g < kNGroups; g++) {
    const Kernel &K = res.kernel[g];
    if (!K.ok)
      continue;
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    Double_t ylo = 0.0, yhi = 0.0;
    for (Int_t a = 0; a < K.n_amp; a++)
      for (Int_t b = 0; b < kNBins; b++) {
        ylo = TMath::Min(ylo, K.k[a][b]);
        yhi = TMath::Max(yhi, K.k[a][b]);
      }
    TH1F *frame = c->DrawFrame(xlo, ylo - 0.02, xhi, yhi + 0.02);
    frame->SetTitle(";log_{10}(#Deltat [#mus]);"
                    "Relative Amplitude Shift");
    TLegend *leg = nullptr;
    if (K.n_amp > 1) {
      leg = PlottingUtils::AddLegend(0.55, 0.89, 0.15, 0.40);
    }
    for (Int_t a = 0; a < K.n_amp; a++) {
      TGraph *gr = new TGraph();
      for (Int_t b = 0; b < kNBins; b++)
        gr->SetPoint(b, xlo + (b + 0.5) * (xhi - xlo) / kNBins, K.k[a][b]);
      gr->SetMarkerStyle(20);
      gr->SetMarkerColor(colors[a % kMaxAmpBins]);
      gr->SetLineColor(colors[a % kMaxAmpBins]);
      gr->SetLineWidth(2);
      gr->Draw("PL SAME");
      if (K.n_amp > 1) {
        leg->AddEntry(gr, Form("Previous pulse ~%d#times beam", a), "pl");
      }
    }
    if (K.n_amp > 1) {
      leg->Draw();
    }
    TLine *zero = new TLine(xlo, 0.0, xhi, 0.0);
    zero->SetLineStyle(2);
    zero->Draw();
    PlottingUtils::SaveFigure(c, Form("kernel_%s", GroupTag(g)), subdir,
                              PlotSaveOptions::kLINEAR);
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
    if (res.shift[g] && res.shift[g]->GetEntries() > 0) {
      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kTRUE);
      PlottingUtils::ConfigureAndDrawHistogram(res.shift[g], kBlack);
      PlottingUtils::SaveFigure(c, Form("applied_shift_%s", gn), subdir,
                                PlotSaveOptions::kLOG);
      delete c;
    }
  }
  for (Int_t g = 0; g < kNGroups; g++) {
    delete res.dev_vs_pred[g];
    delete res.dev_before[g];
    delete res.dev_after[g];
    delete res.shift[g];
    delete res.dtprev_before[g];
    delete res.dtprev_after[g];
    res.dev_vs_pred[g] = nullptr;
    res.dev_before[g] = res.dev_after[g] = res.shift[g] = nullptr;
    res.dtprev_before[g] = res.dtprev_after[g] = nullptr;
  }
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
  TTree *t = new TTree("pulse_history", "Pulse-history kernel per group");
  Int_t group = 0, n_amp = 1;
  Bool_t ok = kFALSE;
  Double_t k[kMaxAmpBins * kNBins], centre_us[kNBins],
      intercept = 0.0, r2 = 0.0, rms_before = 0.0, rms_after = 0.0,
      mean_shift = 0.0, apply_max_us = 0.0;
  Long64_t n = 0, n_beam = 0, n_corrected = 0, n_clamped = 0;
  t->Branch("Group", &group, "Group/I");
  t->Branch("Ok", &ok, "Ok/O");
  t->Branch("NAmpBins", &n_amp, "NAmpBins/I");
  t->Branch("Kernel", k, Form("Kernel[%d]/D", kMaxAmpBins * kNBins));
  t->Branch("BinCentreUs", centre_us, Form("BinCentreUs[%d]/D", kNBins));
  t->Branch("Intercept", &intercept, "Intercept/D");
  t->Branch("R2", &r2, "R2/D");
  t->Branch("RmsBefore", &rms_before, "RmsBefore/D");
  t->Branch("RmsAfter", &rms_after, "RmsAfter/D");
  t->Branch("MeanShift", &mean_shift, "MeanShift/D");
  t->Branch("ApplyMaxUs", &apply_max_us, "ApplyMaxUs/D");
  t->Branch("N", &n, "N/L");
  t->Branch("NBeamEvents", &n_beam, "NBeamEvents/L");
  t->Branch("NCorrected", &n_corrected, "NCorrected/L");
  t->Branch("NClamped", &n_clamped, "NClamped/L");
  for (Int_t g = 1; g < kNGroups; g++) {
    const Kernel &K = res.kernel[g];
    group = g;
    ok = K.ok;
    n_amp = K.n_amp;
    for (Int_t a = 0; a < kMaxAmpBins; a++)
      for (Int_t b = 0; b < kNBins; b++)
        k[a * kNBins + b] = K.k[a][b];
    for (Int_t b = 0; b < kNBins; b++)
      centre_us[b] = BinCentreUs(b);
    intercept = K.intercept;
    r2 = K.r2;
    rms_before = K.rms_before;
    rms_after = K.rms_after;
    mean_shift = res.mean_shift[g];
    apply_max_us = Constants::cfg.PULSE_HISTORY_APPLY_MAX_US;
    n = K.n;
    n_beam = res.n_beam_events;
    n_corrected = res.n_corrected;
    n_clamped = res.n_clamped;
    t->Fill();
  }
  t->Write("pulse_history", TObject::kOverwrite);
  f->Close();
  delete f;
}

} // namespace PulseHistory
