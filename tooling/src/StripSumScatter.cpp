#include "StripSumScatter.hpp"
#include "RegionCuts.hpp"
#include <TParameter.h>
#include <algorithm>
#include <fstream>

const char *const kTagCutName[kNTagCuts] = {
    "tagged",     "all strips", "upstream beam", "jump",
    "reac level", "smoothness", "tail rise",     "re-rise",
    "post above", "end strip",  "cliff"};
const char *const kPreCutName[kNPreCuts] = {"reached tag", "beam gate",
                                            "pileup", "noise", "both mult"};

StripSumScatter::StripSumScatter() {
  for (Int_t i = 0; i < 64; i++) {
    m_yLo[i] = 0.0;
    m_yHi[i] = 0.0;
  }
  m_nSeen = 0;
  m_nNormed = 0;
  m_normedAt.assign(
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX -
          Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN + 1,
      0);
  m_tagged.assign(
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX -
          Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN + 1,
      0);
  m_preCounts.assign(kNPreCuts, 0);
  m_cutCounts.assign(m_tagged.size() * kNTagCuts, 0);
}

StripSumScatter::~StripSumScatter() {
  std::map<Int_t, TH2F *>::iterator it;
  for (it = m_scatter.begin(); it != m_scatter.end(); ++it)
    delete it->second;
  m_scatter.clear();
}

Int_t StripSumScatter::ReacIndex(Int_t reac) {
  return reac - Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
}

Int_t StripSumScatter::YLoOf(Int_t reac) { return reac + 1; }

Int_t StripSumScatter::YHiOf(Int_t reac) {
  const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;
  // The window shrinks where it would otherwise run past the last strip it
  // may reach, and is never shorter than one strip.
  const Int_t hi = TMath::Min(reac + C.POST_TRIGGER_SUM_STRIPS,
                              TMath::Min(C.POST_WINDOW_LAST_STRIP, 17));
  return TMath::Max(hi, reac + 1);
}

Double_t StripSumScatter::s_stripSigma[18] = {0.0};
Double_t StripSumScatter::s_stripMean[18] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                             1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                             1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

Double_t StripSumScatter::StripSigma(Int_t strip) {
  return (strip >= 0 && strip < 18) ? s_stripSigma[strip] : 0.0;
}

Double_t StripSumScatter::StripMean(Int_t strip) {
  return (strip >= 0 && strip < 18) ? s_stripMean[strip] : 1.0;
}

void StripSumScatter::SetStripNoise(const Double_t *mean,
                                    const Double_t *sigma) {
  for (Int_t s = 0; s < 18; s++) {
    s_stripMean[s] = mean ? mean[s] : 1.0;
    s_stripSigma[s] = sigma[s];
  }
}

TagThresholds StripSumScatter::NominalThresholds() {
  const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;
  TagThresholds T;
  T.require_upstream = C.REQUIRE_BEAM_UPSTREAM_OF_REAC;
  T.upstream_nsigma = C.BEAM_UPSTREAM_NSIGMA;
  T.jump_nsigma = C.REAC_JUMP_NSIGMA;
  T.smooth_nsigma = C.SMOOTHNESS_NSIGMA;
  T.tail_fall_from_strip = C.TAIL_FALL_FROM_STRIP;
  T.tail_rise_nsigma = C.TAIL_RISE_NSIGMA;
  T.tail_return_nsigma = C.TAIL_RETURN_NSIGMA;
  T.tail_rerise_nsigma = C.TAIL_RERISE_NSIGMA;
  T.post_above_nsigma = C.POST_ABOVE_NSIGMA;
  T.post_above_strips = C.POST_ABOVE_STRIPS;
  T.end_strip_nsigma = C.END_STRIP_NSIGMA;
  T.cliff_max = C.TAIL_CLIFF_MAX_FRACTION;
  return T;
}

std::vector<std::pair<TString, TagThresholds>>
StripSumScatter::ThresholdVariants() {
  const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;
  std::vector<std::pair<TString, TagThresholds>> out;
  if (!C.CUT_VARIATION)
    return out;
  const TagThresholds nom = NominalThresholds();
  const Double_t ds = C.CUT_VARIATION_NSIGMA_STEP;
  // One threshold at a time, up and down; only the conditions that are on.
  // The sign convention is the threshold's own: "+" is the larger value.
  struct Var {
    const char *name;
    Double_t TagThresholds::*field;
    Bool_t on;
    Double_t step;
    Bool_t any_sign; // the threshold stays meaningful at or below zero
  };
  const Var vars[] = {
      {"upstream", &TagThresholds::upstream_nsigma, nom.require_upstream, ds,
       kFALSE},
      {"jump", &TagThresholds::jump_nsigma, kTRUE, ds, kFALSE},
      {"smooth", &TagThresholds::smooth_nsigma, nom.smooth_nsigma > 0.0, ds,
       kFALSE},
      {"tailrise", &TagThresholds::tail_rise_nsigma,
       nom.tail_rise_nsigma > 0.0 && nom.tail_fall_from_strip > 0, ds, kFALSE},
      {"return", &TagThresholds::tail_return_nsigma,
       nom.tail_rerise_nsigma > 0.0, ds, kFALSE},
      {"rerise", &TagThresholds::tail_rerise_nsigma,
       nom.tail_rerise_nsigma > 0.0, ds, kFALSE},
      {"post", &TagThresholds::post_above_nsigma,
       nom.post_above_nsigma > 0.0 && nom.post_above_strips > 0, ds, kFALSE},
      // The end strip is always tested; its threshold, mean - n sigma, is
      // meaningful for any n, so the shift is never skipped.
      {"end", &TagThresholds::end_strip_nsigma, kTRUE, ds, kTRUE},
      {"cliff", &TagThresholds::cliff_max, nom.cliff_max > 0.0,
       C.CUT_VARIATION_CLIFF_STEP, kFALSE},
  };
  for (size_t v = 0; v < sizeof(vars) / sizeof(vars[0]); v++) {
    if (!vars[v].on || !(vars[v].step > 0.0))
      continue;
    for (Int_t dir = +1; dir >= -1; dir -= 2) {
      TagThresholds T = nom;
      T.*(vars[v].field) = nom.*(vars[v].field) + dir * vars[v].step;
      // A variation must not switch a condition off (or on): keep it above 0
      // where zero means off.
      if (!vars[v].any_sign && !(T.*(vars[v].field) > 0.0))
        continue;
      out.push_back(
          std::make_pair(TString(vars[v].name) + (dir > 0 ? "+" : "-"), T));
    }
  }
  return out;
}

Double_t StripSumScatter::JumpMin(Int_t reac) {
  return Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REAC_JUMP_NSIGMA *
         StripSigma(reac);
}

Bool_t StripSumScatter::BeamUpstreamOf(const EnergyView &ev, Int_t reac,
                                       const TagThresholds &T) {
  if (!T.require_upstream)
    return kTRUE;
  for (Int_t s = 1; s < reac; s++)
    if (TMath::Abs(ev.Total(s) - StripMean(s)) >
        T.upstream_nsigma * StripSigma(s))
      return kFALSE;
  return kTRUE;
}

Bool_t StripSumScatter::BeamUpstreamOf(const EnergyView &ev, Int_t reac) {
  return BeamUpstreamOf(ev, reac, NominalThresholds());
}

/// The tail-shape conditions of the published 87Rb per-strip macros
/// (TracesVisu2.C at A9-A11), each in sigma of the measured noise and off
/// unless its config value is set. All are numerator-only: they describe the
/// residue, not the beam. A strip whose noise was not measured is skipped.
TagCut StripSumScatter::TailReason(const EnergyView &ev, Int_t reac,
                                   const TagThresholds &T) {
  const Int_t kLast = Constants::cfg.IGNORE_STRIP_17 ? 16 : 17;
  // Falling tail: no rise beyond the noise from the start strip on.
  if (T.tail_rise_nsigma > 0.0 && (T.tail_fall_from_strip > 0)) {
    Int_t start = reac + 1;
    start = TMath::Max(reac + 1, T.tail_fall_from_strip + 1);

    for (Int_t s = start; s <= kLast; s++)
      if (StripSigma(s) > 0.0 &&
          ev.Total(s) - ev.Total(s - 1) > T.tail_rise_nsigma * StripSigma(s))
        return kCutTailRise;
  }
  // No return: back at the beam after the peak, then above it again.
  if (T.tail_rerise_nsigma > 0.0) {
    Int_t peak = reac;
    for (Int_t s = reac + 1; s <= TMath::Min(YHiOf(reac), kLast); s++)
      if (ev.Total(s) > ev.Total(peak))
        peak = s;
    Bool_t returned = kFALSE;
    for (Int_t s = peak + 1; s <= kLast; s++) {
      if (!(StripSigma(s) > 0.0))
        continue;
      if (returned &&
          ev.Total(s) > StripMean(s) + T.tail_rerise_nsigma * StripSigma(s))
        return kCutRerise;
      if (ev.Total(s) <= StripMean(s) + T.tail_return_nsigma * StripSigma(s))
        returned = kTRUE;
    }
  }
  // Persistence: the excess must hold for post_above_strips strips after the
  // reaction, each more than post_above_nsigma of its spread above the beam.
  if (T.post_above_nsigma > 0.0 && T.post_above_strips > 0)
    for (Int_t s = reac + 1; s <= TMath::Min(reac + T.post_above_strips, kLast);
         s++)
      if (StripSigma(s) > 0.0 &&
          !(ev.Total(s) > StripMean(s) + T.post_above_nsigma * StripSigma(s)))
        return kCutPostAbove;
  return kTagPass;
}

TagCut StripSumScatter::TailReason(const EnergyView &ev, Int_t reac) {
  return TailReason(ev, reac, NominalThresholds());
}

Bool_t StripSumScatter::PassesTail(const EnergyView &ev, Int_t reac) {
  return TailReason(ev, reac) == kTagPass;
}

TagCut StripSumScatter::RejectReason(const EnergyView &ev, Int_t reac,
                                     const TagThresholds &T) {
  const Double_t kReacJumpMin = T.jump_nsigma * StripSigma(reac);
  const Int_t kLast = Constants::cfg.IGNORE_STRIP_17 ? 16 : 17;

  if (!AllStripsFired(ev))
    return kCutAllStrips;
  // Otherwise a reaction at an earlier strip can pass this strip's jump gate
  // on a noise fluctuation and be counted here as well.
  if (!BeamUpstreamOf(ev, reac, T))
    return kCutUpstream;
  Double_t reac_jump = ev.Total(reac) - ev.Total(reac - 1);
  if (!(reac_jump > kReacJumpMin))
    return kCutJump;
  if (!(ev.Total(reac) > StripMean(reac) + kReacJumpMin))
    return kCutReacLevel;
  if (T.smooth_nsigma > 0.0)
    for (Int_t s = reac + 1; s <= kLast; s++)
      if (StripSigma(s) > 0.0 && TMath::Abs(ev.Total(s) - ev.Total(s - 1)) >
                                     T.smooth_nsigma * StripSigma(s))
        return kCutSmooth;
  const TagCut tail = TailReason(ev, reac, T);
  if (tail != kTagPass)
    return tail;
  Int_t end_strip = 17;
  if (Constants::cfg.IGNORE_STRIP_17 ||
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REQUIRE_STRIP_16_BELOW_BEAM)
    end_strip = 16;
  // The end strip below the beam by n sigma of its own spread: a residue
  // that has stopped (or is stopping) reads low there, the beam does not.
  if (!(ev.Total(end_strip) <
        StripMean(end_strip) - T.end_strip_nsigma * StripSigma(end_strip)))
    return kCutEndStrip;
  // A stop is gradual: the last step may take only part of the fall from the
  // peak. The elastic class keeps its excess to end_strip-1 and drops there.
  if (T.cliff_max > 0.0 && end_strip - 1 > reac) {
    Double_t peak = ev.Total(reac);
    for (Int_t s = reac + 1; s < end_strip; s++)
      peak = TMath::Max(peak, ev.Total(s));
    const Double_t fall = peak - ev.Total(end_strip);
    const Double_t last = ev.Total(end_strip - 1) - ev.Total(end_strip);
    if (fall > 0.0 && last / fall > T.cliff_max)
      return kCutCliff;
  }
  return kTagPass;
}

TagCut StripSumScatter::RejectReason(const EnergyView &ev, Int_t reac) {
  return RejectReason(ev, reac, NominalThresholds());
}

// Sigma-clipped mean and width of a sample, started from the median and the
// MAD. The start matters even on a beam-gated sample: what the gate lets
// through still carries the odd partial pile-up and reaction, and a clip that
// starts from the plain RMS of a sample with a second population never sheds
// it (before the gate was used here it settled at 0.41 a.u. against a beam
// width of 0.06 on 37Cl run 97, which turned every sigma-scaled cut into a
// no-op). The MAD sees such a population as a minority and starts inside the
// beam.
static Double_t ClippedWidth(const std::vector<Double_t> &values,
                             Double_t *mean_out = nullptr) {
  const Int_t kClipPasses = 3;
  const Double_t kClipNSigma = 3.0;
  if (mean_out)
    *mean_out = 0.0;
  if (values.size() < 2)
    return 0.0;
  std::vector<Double_t> sorted(values);
  std::sort(sorted.begin(), sorted.end());
  Double_t mean = sorted[sorted.size() / 2];
  for (size_t k = 0; k < sorted.size(); k++)
    sorted[k] = TMath::Abs(sorted[k] - mean);
  std::sort(sorted.begin(), sorted.end());
  Double_t width = 1.4826 * sorted[sorted.size() / 2];
  if (!(width > 0.0))
    return 0.0;
  for (Int_t pass = 0; pass < kClipPasses; pass++) {
    Double_t sum = 0.0, sum2 = 0.0;
    Long64_t kept = 0;
    for (Int_t k = 0; k < Int_t(values.size()); k++) {
      const Double_t v = values[k];
      if (TMath::Abs(v - mean) > kClipNSigma * width)
        continue;
      sum += v;
      sum2 += v * v;
      kept++;
    }
    if (kept < 2)
      break;
    mean = sum / Double_t(kept);
    width = TMath::Sqrt(TMath::Max(0.0, sum2 / Double_t(kept) - mean * mean));
  }
  if (mean_out)
    *mean_out = mean;
  return width;
}

Bool_t StripSumScatter::MeasureBeamNoise(TChain *chain, const BeamEllipses &be,
                                         Double_t *strip_mean,
                                         Double_t *strip_sigma) {
  // Enough gated events for a percent-level width on every strip; the gate
  // passes about half of what it sees.
  const Long64_t kMaxGated = 100000;
  const Long64_t kMinEvents = 1000;
  for (Int_t s = 0; s < 18; s++) {
    strip_mean[s] = 1.0;
    strip_sigma[s] = 0.0;
  }
  if (!chain || !be.ok)
    return kFALSE;

  EnergyView ev;
  ev.Attach(chain);
  EnableEventBranches(chain);
  const Long64_t n = chain->GetEntries();
  std::vector<std::vector<Double_t>> deposit(18);
  for (Long64_t j = 0; j < n && Long64_t(deposit[1].size()) < kMaxGated; j++) {
    chain->GetEntry(j);
    ev.Decode();
    // The beam sample: every strip fired and inside the entrance and exit
    // ellipses, the same selection the scatters call pure beam.
    if (!IsPureBeam(ev, be))
      continue;
    for (Int_t s = 0; s < 18; s++)
      deposit[s].push_back(ev.Total(s));
  }
  chain->ResetBranchAddresses();
  if (Long64_t(deposit[1].size()) < kMinEvents)
    return kFALSE;

  // Per strip, the mean and width of that sample, clipped for the residual
  // partial pile-up and reactions the end gates do not see.
  for (Int_t s = 0; s < 18; s++) {
    Double_t mean = 0.0;
    strip_sigma[s] = ClippedWidth(deposit[s], &mean);
    if (mean > 0.0)
      strip_mean[s] = mean;
  }
  for (Int_t s = 1; s <= 16; s++)
    if (!(strip_sigma[s] > 0.0))
      return kFALSE;
  return kTRUE;
}

void StripSumScatter::EnableEventBranches(TChain *chain) {
  chain->SetBranchStatus("*", 0);
  chain->SetBranchStatus("LeftdE", 1);
  chain->SetBranchStatus("RightdE", 1);
  chain->SetBranchStatus("Strip0dE", 1);
  chain->SetBranchStatus("Strip17dE", 1);
  chain->SetBranchStatus("Cathode", 1);
  // Absent in files built before the branch existed; harmless to enable.
  if (chain->GetBranch("SeedTs"))
    chain->SetBranchStatus("SeedTs", 1);
}

Bool_t StripSumScatter::AllStripsFired(const EnergyView &ev) {
  if (!Constants::cfg.IGNORE_STRIP_0 && !(ev.strip0 > 0.0))
    return kFALSE;
  if (!Constants::cfg.IGNORE_STRIP_17 && !(ev.strip17 > 0.0))
    return kFALSE;
  for (Int_t s = 1; s <= 16; s++)
    if (!(ev.Total(s) > 0.0))
      return kFALSE;
  return kTRUE;
}

Bool_t StripSumScatter::PassesReaction(const EnergyView &ev, Int_t reac) {
  return RejectReason(ev, reac) == kTagPass;
}

Bool_t StripSumScatter::IsPureBeam(const EnergyView &ev,
                                   const BeamEllipses &be) {
  if (!be.ok)
    return kFALSE;
  if (!AllStripsFired(ev))
    return kFALSE;
  // Must pass BOTH the entrance AND exit ellipses.
  if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PURE_BEAM_GATE ==
      StripSumScatterConfig::PURE_BEAM_GATE_S1_S2) {
    if (!PassesGate(be.s1_s2, ev, 1, 2))
      return kFALSE;
  } else {
    if (!PassesGate(be.s0_s1, ev, 0, 1))
      return kFALSE;
  }
  if (be.use_s15_s16) {
    if (!PassesGate(be.s15_s16, ev, 15, 16))
      return kFALSE;
  } else {
    if (!PassesGate(be.s16_s17, ev, 16, 17))
      return kFALSE;
  }
  return kTRUE;
}

Bool_t StripSumScatter::IsPileup(const EnergyView &ev) {
  const Double_t kNSigma =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PILEUP_NSIGMA;
  const Int_t kMinStrips =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PILEUP_MIN_STRIPS;
  if (!(kNSigma > 0.0))
    return kFALSE;
  Int_t n = 0;
  // A strip whose spread is not yet measured cannot be judged: the cut is
  // off there, which is what makes MeasureBeamNoise itself unaffected by it.
  for (Int_t s = 1; s <= 16; s++)
    if (StripSigma(s) > 0.0 &&
        ev.Total(s) >= StripMean(s) + kNSigma * StripSigma(s) &&
        ++n >= kMinStrips)
      return kTRUE;
  return kFALSE;
}

Bool_t StripSumScatter::IsNoise(const EnergyView &ev) {
  const Double_t kNSigma = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.NOISE_NSIGMA;
  const Int_t kMinStrips =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.NOISE_MIN_STRIPS;
  if (!(kNSigma > 0.0))
    return kFALSE;
  Int_t n = 0;
  for (Int_t s = 1; s <= 16; s++)
    if (StripSigma(s) > 0.0 &&
        ev.Total(s) <= StripMean(s) - kNSigma * StripSigma(s) &&
        ++n >= kMinStrips)
      return kTRUE;
  return kFALSE;
}

Double_t StripSumScatter::SumRange(const Double_t *total, Int_t lo, Int_t hi) {
  Double_t sum = 0.0;
  for (Int_t s = lo; s <= hi; s++)
    sum += total[s];
  return sum;
}

void StripSumScatter::PlaneXY(const Double_t *total, Int_t reac, Double_t &x,
                              Double_t &y) {
  x = SumRange(total, Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_LO,
               Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_HI);
  y = SumRange(total, YLoOf(reac), YHiOf(reac));
  // The published macros' `ratio`: the post window over the event's own
  // upstream mean, so per-event jitter common to every strip cancels.
  if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.Y_RATIO_TO_UPSTREAM && reac > 1) {
    const Double_t up = SumRange(total, 1, reac - 1) / Double_t(reac - 1);
    if (up > 0.0)
      y /= up;
  }
}

std::vector<GateSpec> StripSumScatter::ActiveGates() {
  std::vector<GateSpec> gates;
  GateSpec g;
  g.sx = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_STRIP_X;
  g.sy = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_STRIP_Y;
  gates.push_back(g);
  if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REQUIRE_GATE_S3_S4) {
    g.sx = 3;
    g.sy = 4;
    gates.push_back(g);
  }
  if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REQUIRE_GATE_S5_S6) {
    g.sx = 5;
    g.sy = 6;
    gates.push_back(g);
  }
  return gates;
}

TString StripSumScatter::CacheName() {
  TString name = "StripSumScatter_cache";
  if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REQUIRE_GATE_S3_S4)
    name += "_g34";
  if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REQUIRE_GATE_S5_S6)
    name += "_g56";
  name += ".root";
  return name;
}

Bool_t StripSumScatter::PassesGate(const BeamFit2D &gate, const EnergyView &ev,
                                   Int_t sx, Int_t sy) {
  Double_t g0 = ev.Total(sx);
  Double_t g1 = ev.Total(sy);
  if (!(g0 > 0.0 && g1 > 0.0))
    return kFALSE;
  const Double_t kGateNSigmaX =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_X;
  const Double_t kGateNSigmaY =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_Y;
  return BeamFitUtils::InEllipseXY(gate, g0, g1, kGateNSigmaX, kGateNSigmaY);
}

BeamFit2D
StripSumScatter::FindBeamGate(TChain *chain, Int_t sx, Int_t sy,
                              const std::vector<GateSpec> &prior_specs,
                              const std::vector<BeamFit2D> &prior_gates,
                              const TString &tag, const TString &subdir) {
  const Int_t kGateBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_BINS;
  const Double_t kGateMin = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_MIN;
  const Double_t kGateMax = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_MAX;
  const Int_t kSeedHalfBins =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.SEED_HALF_BINS;
  const Double_t kSeedFrac = 0.3;
  const Long64_t kSampleMaxPoints =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.SAMPLE_MAX_POINTS;

  BeamFit2D out;
  EnergyView ev;
  ev.Attach(chain);
  EnableEventBranches(chain);
  TH2F *h =
      new TH2F(Form("h2_beamgate_s%d_s%d_%s", sx, sy, tag.Data()),
               Form(";#DeltaE strip %d [a.u.];#DeltaE strip %d [a.u.]", sx, sy),
               kGateBins, kGateMin, kGateMax, kGateBins, kGateMin, kGateMax);
  h->SetDirectory(nullptr);
  Long64_t n = chain->GetEntries();
  Long64_t stride = FileSet::SampleStride(n, kSampleMaxPoints);
  for (Long64_t j = 0; j < n; j += stride) {
    chain->GetEntry(j);
    ev.Decode();
    // Series gating: only events passing every prior gate feed this fit.
    Bool_t prior_ok = kTRUE;
    for (Int_t gi = 0; gi < Int_t(prior_specs.size()); gi++)
      if (!PassesGate(prior_gates[gi], ev, prior_specs[gi].sx,
                      prior_specs[gi].sy)) {
        prior_ok = kFALSE;
        break;
      }
    if (!prior_ok)
      continue;
    Double_t x = ev.Total(sx);
    Double_t y = ev.Total(sy);
    if (x > 0.0 && y > 0.0)
      h->Fill(x, y);
  }
  if (h->GetEntries() < 100) {
    delete h;
    return out;
  }
  Double_t bw_x = h->GetXaxis()->GetBinWidth(1);
  Double_t bw_y = h->GetYaxis()->GetBinWidth(1);
  Int_t bx = 0, by = 0, bz = 0;
  h->GetMaximumBin(bx, by, bz);
  Double_t peak_val = h->GetBinContent(bx, by);
  Int_t lo_bx = std::max(1, bx - kSeedHalfBins);
  Int_t hi_bx = std::min(h->GetNbinsX(), bx + kSeedHalfBins);
  Int_t lo_by = std::max(1, by - kSeedHalfBins);
  Int_t hi_by = std::min(h->GetNbinsY(), by + kSeedHalfBins);
  Moments2D m = BeamFitUtils::ComputeMoments(h, lo_bx, hi_bx, lo_by, hi_by,
                                             kSeedFrac * peak_val, bw_x, bw_y);
  if (m.weight <= 0) {
    delete h;
    return out;
  }
  out.amp = peak_val;
  out.mu_x = m.mu_x;
  out.mu_y = m.mu_y;
  out.sigma_x = m.sigma_x;
  out.sigma_y = m.sigma_y;
  out.rho = m.rho;
  out.ok = kTRUE;

  const Double_t kGateNSigmaX =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_X;
  const Double_t kGateNSigmaY =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_Y;

  // The per-run beam-gate figures are the only thing written under
  // strip_sum_scatter/run<N>; skipping them skips those folders.
  if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.SKIP_RUN_PLOTS) {
    delete h;
    return out;
  }
  std::lock_guard<std::mutex> lock(g_plot_mutex);
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDraw2DHistogram(h, c);
  // Correlated 2D Gaussian ellipse matching the InEllipseXY gate.
  Double_t sxx = out.sigma_x * out.sigma_x;
  Double_t syy = out.sigma_y * out.sigma_y;
  Double_t sxy = out.rho * out.sigma_x * out.sigma_y;
  Double_t sum = sxx + syy;
  Double_t diff = sxx - syy;
  Double_t det = TMath::Sqrt(diff * diff + 4.0 * sxy * sxy);
  Double_t lambda1 = 0.5 * (sum + det);
  Double_t lambda2 = 0.5 * (sum - det);
  Double_t theta = 0.5 * TMath::ATan2(2.0 * sxy, diff) * 180.0 / TMath::Pi();
  TEllipse *e =
      new TEllipse(out.mu_x, out.mu_y, kGateNSigmaX * TMath::Sqrt(lambda1),
                   kGateNSigmaX * TMath::Sqrt(lambda2), 0, 360, theta);
  e->SetFillStyle(0);
  e->SetLineColor(kRed + 1);
  e->SetLineWidth(2);
  e->Draw();
  PlottingUtils::SaveFigure(c, Form("beam_gate_s%d_s%d", sx, sy), subdir,
                            PlotSaveOptions::kLINEAR);
  delete c;

  delete h;
  return out;
}

void StripSumScatter::DrawTraceSet(const std::vector<TGraph *> &traces,
                                   Int_t color) {
  for (Int_t i = 0; i < Int_t(traces.size()); i++) {
    traces[i]->SetLineColor(color);
    traces[i]->SetLineWidth(1);
    traces[i]->Draw("L SAME");
  }
}

TGraph *StripSumScatter::TraceFromTotal(const Double_t *total) {
  return EventsSummary::BuildTraceFromTotals(total);
}

void StripSumScatter::DrawRegionTraces(
    const TString &save_name, const TString &subdir,
    const std::vector<TGraph *> &beam, const std::vector<TGraph *> &aa,
    const std::vector<TGraph *> &an, Double_t y_min, Double_t y_max,
    const char *y_title, Bool_t beam_sigma_measured, const char *an_label) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);
  Int_t s_lo = Constants::cfg.IGNORE_STRIP_0 ? 1 : 0;
  Int_t s_hi = Constants::cfg.IGNORE_STRIP_17 ? 16 : 17;
  TH2F *frame =
      new TH2F("h_region_trace_frame", Form(";Strip;%s", y_title),
               s_hi - s_lo + 1, s_lo - 0.5, s_hi + 0.5, 100, y_min, y_max);
  frame->SetStats(0);
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  frame->Draw();

  // Beam: mean of the sampled beam traces with a +-1 sigma band, so the
  // region traces read against the noise the sigma cuts resolve through.
  TGraphErrors *beam_band = nullptr;
  if (!beam.empty()) {
    const Int_t npts = beam[0]->GetN();
    std::vector<Double_t> mean(npts, 0.0), m2(npts, 0.0);
    for (Int_t t = 0; t < Int_t(beam.size()); t++) {
      const Double_t *yv = beam[t]->GetY();
      for (Int_t p = 0; p < npts; p++) {
        mean[p] += yv[p];
        m2[p] += yv[p] * yv[p];
      }
    }
    const Double_t nt = Double_t(beam.size());
    const Double_t *xv = beam[0]->GetX();
    beam_band = new TGraphErrors(npts);
    for (Int_t p = 0; p < npts; p++) {
      mean[p] /= nt;
      const Double_t var = m2[p] / nt - mean[p] * mean[p];
      const Double_t rms = var > 0.0 ? TMath::Sqrt(var) : 0.0;
      const Int_t strip = Int_t(TMath::Nint(xv[p]));
      const Double_t measured = StripSigma(strip) * mean[p];
      const Double_t sigma =
          beam_sigma_measured && measured > 0.0 ? measured : rms;
      beam_band->SetPoint(p, xv[p], mean[p]);
      beam_band->SetPointError(p, 0.0, sigma);
    }
    // Dark mean line over a light band.
    beam_band->SetLineColor(kGray + 3);
    beam_band->SetLineWidth(3);
    beam_band->SetFillColorAlpha(kGray + 1, 0.35);
    beam_band->Draw("3 SAME");
  }
  DrawTraceSet(aa, kAzure + 2);
  DrawTraceSet(an, kRed + 1);
  if (beam_band)
    beam_band->Draw("LX SAME"); // mean line above the traces

  TGraph *p_aa = new TGraph(1);
  TGraph *p_an = new TGraph(1);
  TGraph *proxies[2] = {p_aa, p_an};
  Int_t pcol[2] = {kAzure + 2, kRed + 1};
  for (Int_t i = 0; i < 2; i++) {
    proxies[i]->SetPoint(0, -1e9, -1e9);
    proxies[i]->SetLineColor(pcol[i]);
    proxies[i]->SetLineWidth(3);
  }
  TLegend *leg = PlottingUtils::AddLegend(0.725, 0.875, 0.70, 0.86);
  if (beam_band)
    leg->AddEntry(beam_band, "Beam #pm 1#sigma", "lf");
  // An empty class is left out of the legend: the all-tagged overlay has
  // one class of traces against the beam.
  if (!aa.empty())
    leg->AddEntry(p_aa, "(#alpha,#alpha')", "l");
  if (!an.empty())
    leg->AddEntry(p_an, an_label ? an_label : "(#alpha,n)", "l");
  leg->Draw();

  PlottingUtils::SaveFigure(c, save_name, subdir, PlotSaveOptions::kLINEAR);
  delete c;
}

void StripSumScatter::DrawRegionMeanTraces(const TString &save_name,
                                           const TString &subdir,
                                           const std::vector<TGraph *> &beam,
                                           const std::vector<TGraph *> &aa,
                                           const std::vector<TGraph *> &an,
                                           Double_t y_min, Double_t y_max,
                                           const char *y_title) {
  std::lock_guard<std::mutex> lock(g_plot_mutex);
  TH2F *frame = new TH2F("h_region_mean_frame", Form(";Strip;%s", y_title), 18,
                         -0.5, 17.5, 100, y_min, y_max);
  frame->SetStats(0);
  frame->SetDirectory(nullptr);
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  frame->Draw();

  const std::vector<TGraph *> *regions[3] = {&beam, &aa, &an};
  // Dark mean lines over light bands of the same hue.
  Int_t colors[3] = {kGray + 3, kAzure + 3, kRed + 2};
  Int_t bands[3] = {kGray + 1, kAzure + 1, kRed - 7};
  const char *labels[3] = {"Beam", "(#alpha,#alpha')", "(#alpha,n)"};
  std::vector<TGraphErrors *> means;
  TLegend *leg = PlottingUtils::AddLegend(0.725, 0.875, 0.70, 0.86);

  for (Int_t r = 0; r < 3; r++) {
    const std::vector<TGraph *> &tr = *regions[r];
    if (tr.empty())
      continue;
    Int_t npts = tr[0]->GetN();
    std::vector<Double_t> mean(npts, 0.0), m2(npts, 0.0);
    for (Int_t t = 0; t < Int_t(tr.size()); t++) {
      Double_t *yv = tr[t]->GetY();
      for (Int_t p = 0; p < npts; p++) {
        mean[p] += yv[p];
        m2[p] += yv[p] * yv[p];
      }
    }
    Double_t nt = Double_t(tr.size());
    TGraphErrors *ge = new TGraphErrors(npts);
    Double_t *xv = tr[0]->GetX();
    for (Int_t p = 0; p < npts; p++) {
      mean[p] /= nt;
      Double_t var = m2[p] / nt - mean[p] * mean[p];
      ge->SetPoint(p, xv[p], mean[p]);
      ge->SetPointError(p, 0.0, var > 0.0 ? TMath::Sqrt(var) : 0.0);
    }
    ge->SetLineColor(colors[r]);
    ge->SetLineWidth(3);
    ge->SetFillColorAlpha(bands[r], 0.35);
    ge->Draw("3 SAME"); // +-1 RMS band (all bands first, behind the lines)
    means.push_back(ge);
    leg->AddEntry(ge, labels[r], "lf");
  }
  // Mean lines on top of every band.
  for (Int_t i = 0; i < Int_t(means.size()); i++)
    means[i]->Draw("LX SAME"); // mean line, no end caps
  leg->Draw();

  PlottingUtils::SaveFigure(c, save_name, subdir, PlotSaveOptions::kLINEAR);
  for (Int_t i = 0; i < Int_t(means.size()); i++)
    delete means[i];
  delete leg;
  delete c;
  delete frame;
}

void StripSumScatter::TraceYRange(const std::vector<TGraph *> &beam,
                                  const std::vector<TGraph *> &aa,
                                  const std::vector<TGraph *> &an,
                                  Double_t &y_min, Double_t &y_max) {
  y_min = std::numeric_limits<Double_t>::max();
  y_max = -std::numeric_limits<Double_t>::max();
  const std::vector<TGraph *> *sets[3] = {&beam, &aa, &an};
  for (Int_t si = 0; si < 3; si++) {
    const std::vector<TGraph *> &v = *sets[si];
    for (Int_t i = 0; i < Int_t(v.size()); i++) {
      Double_t x = 0.0, y = 0.0;
      for (Int_t k = 0; k < v[i]->GetN(); k++) {
        v[i]->GetPoint(k, x, y);
        if (x < 0.5 || x > 16.5)
          continue;
        if (y < y_min)
          y_min = y;
        if (y > y_max)
          y_max = y;
      }
    }
  }
  if (y_min > y_max) { // no in-range points sampled
    y_min = 0.0;
    y_max = 1.0;
  }
  Double_t pad = 0.05 * (y_max - y_min);
  if (pad <= 0.0)
    pad = 1.0;
  y_min -= pad;
  y_max += pad;
}

TCutG *StripSumScatter::PromptCut(TCanvas *c, const char *name,
                                  const char *label) {
  std::cout << "  >>> draw the " << label
            << " region: left-click vertices, double-click to close"
            << std::endl;
  c->cd();
  TCutG *cut = static_cast<TCutG *>(c->WaitPrimitive("CUTG", "CutG"));
  if (!cut) {
    std::cerr << "  no " << label << " cut drawn" << std::endl;
    return nullptr;
  }
  cut->SetName(name);
  cut->SetLineColor(kBlack);
  cut->SetLineWidth(2);
  return cut;
}

/// Saved region cuts.
///
/// Drawing the regions by hand is the only way to place them when the reaction
/// population's location is not yet known, but a hand-drawn cut that is not
/// stored makes the run unreproducible: the next pass gets a different polygon
/// and no two results are comparable. Persisting them separates "decide where
/// the region is", which needs a person once, from "apply it", which should be
/// automatic from then on -- and without a DISPLAY.
void StripSumScatter::SaveRegionCuts(Int_t reac, TCutG *cut_an, TCutG *cut_aa) {
  RegionCutStore::Save(reac, cut_an, cut_aa);
}

TCutG *StripSumScatter::LoadRegionCut(const char *name, Int_t reac) {
  return RegionCutStore::Load(name, reac);
}

void StripSumScatter::SmoothTrace(const Double_t *in, Double_t *out,
                                  Int_t width) {
  Int_t half = width / 2;
  for (Int_t s = 0; s < 18; s++) {
    Int_t lo = TMath::Max(0, s - half);
    Int_t hi = TMath::Min(17, s + half);
    Double_t sum = 0.0;
    for (Int_t t = lo; t <= hi; t++)
      sum += in[t];
    out[s] = sum / Double_t(hi - lo + 1);
  }
}

/// Savitzky-Golay smoothing: 3rd-degree polynomial, half-window of 2
/// (5-point convolution). Uses standard SG coefficients that sum to 1.
/// For a 5-point window with 3rd-degree polynomial, the smoothed value at
/// the center uses coefficients: [-3, 12, 17, 12, -3] / 35.
/// At edges, the window shrinks and coefficients are renormalised.
void StripSumScatter::SavitzkyGolay(const Double_t *in, Double_t *out) {
  static const Int_t K = 2; // half-width (5-point window)

  // Standard SG coefficients for 5-point, 3rd-degree polynomial (smoothed
  // value): These are translation-invariant - same for all center positions.
  static const Double_t sg_coeff[2 * K + 1] = {
      -3.0 / 35.0, // coefficient for t = s - 2
      12.0 / 35.0, // coefficient for t = s - 1
      17.0 / 35.0, // coefficient for t = s (center)
      12.0 / 35.0, // coefficient for t = s + 1
      -3.0 / 35.0  // coefficient for t = s + 2
  };

  for (Int_t s = 0; s < 18; s++) {
    Int_t lo = TMath::Max(0, s - K);
    Int_t hi = TMath::Min(17, s + K);
    Double_t val = 0.0;
    Double_t wsum = 0.0;

    // Apply SG coefficients for the positions within the clipped window
    for (Int_t t = lo; t <= hi; t++) {
      Int_t offset = t - s + K; // 0..4, position within 5-point window
      val += sg_coeff[offset] * in[t];
      wsum += sg_coeff[offset];
    }

    // Renormalise at edges where window shrinks
    out[s] = (wsum != 0.0) ? val / wsum : in[s];
  }
}

// CFD-style trigger finder: scan left-to-right for the first strip whose
// beam-subtracted signal beats a peak fraction AND beam sigma; else -1.

Int_t StripSumScatter::FindTrigger(const Double_t *td, const Double_t *base,
                                   Double_t beam_sigma) {
  const Int_t s_lo = 2;
  const Int_t s_hi = 16;

  const Double_t frac =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRIGGER_CFD_FRAC;

  const Double_t nsigma =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRIGGER_NSIGMA * beam_sigma;

  Double_t peak_signal = -1.0e30;

  for (Int_t s = s_lo; s <= s_hi; s++) {
    Double_t signal = td[s] - base[s];
    if (signal > peak_signal)
      peak_signal = signal;
  }

  if (peak_signal < nsigma)
    return -1;

  Double_t thresh = frac * peak_signal;

  for (Int_t s = s_lo; s <= s_hi; s++) {
    Double_t signal = td[s] - base[s];

    if (signal >= thresh && signal >= nsigma)
      return s;
  }

  return -1;
}

// Build a TGraph trace from a Savitzky-Golay-smoothed set of per-strip totals.
TGraph *StripSumScatter::SmoothedTraceFromTotal(const Double_t *total) {
  Double_t sgd[18];
  SavitzkyGolay(total, sgd);
  return EventsSummary::BuildTraceFromTotals(sgd);
}

void StripSumScatter::ClusterVarHists(Int_t reac, TCutG *cut_aa, TCutG *cut_an,
                                      const TString &subdir) {
  const Int_t NV = 9;
  const Int_t NC = 3;
  const char *vkey[NV] = {"energy",      "peak3",      "plateau",
                          "tail",        "reacstrip",  "mult",
                          "trigtaildev", "reacslope3", "beamdev"};
  const char *vtitle[NV] = {
      "#Sigma_{all strips}(#DeltaE#minus1) [a.u.]",
      "#Sigma_{trig#pm1}#DeltaE (0 if no trigger) [a.u.]",
      "Plateau Excess  #Sigma_{trig+1..trig+POST}(#DeltaE#minus1) [a.u.]",
      "#DeltaE(s17) [a.u.]",
      "Trigger Strip",
      "Both-side Multiplicity (strips 1-16)",
      "|#DeltaE#minusbeam| at trigger + at s17 [a.u.]",
      "#DeltaE(reac+3) #minus #DeltaE(reac#minus3) [a.u.]",
      "RMS_{8-17}(#DeltaE#minusbeam) [a.u.]"};
  const char *clabel[NC] = {"beam", "(a,a')", "(a,n)"};

  const Int_t kXLo = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_HI;
  const Int_t kClusterSmoothWindow =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.CLUSTER_SMOOTH_WINDOW;

  // Per (class, variable) value lists for raw traces.
  std::vector<Double_t> vals_raw[NC][NV];
  // Same structure, but cluster variables computed on SG-smoothed traces
  // (SG removes the sawtooth: peak/onset hit real peaks, not odd strips).
  std::vector<Double_t> vals_sg[NC][NV];
  UInt_t bit = (1u << ReacIndex(reac));

  // Per-strip beam baseline = mean of the beam-flat reservoir events; it
  // carries the L_odd/R_even sawtooth, so subtraction removes it exactly.
  Double_t base[18];
  for (Int_t s = 0; s < 18; s++)
    base[s] = 0.0;
  Long64_t nbeam = 0;
  for (Int_t k = 0; k < Int_t(m_reservoir.size()); k++)
    if (m_reservoir[k].beam_flat) {
      for (Int_t s = 0; s < 18; s++)
        base[s] += m_reservoir[k].Total(s);
      nbeam++;
    }
  for (Int_t s = 0; s < 18; s++)
    base[s] = (nbeam > 0) ? base[s] / Double_t(nbeam) : 1.0;

  // Pooled beam-noise RMS = sqrt(mean (total - base)^2 over beam-flat
  // events + strips); onset threshold = this many sigma; flat beam: no trigger.
  Double_t beam_sumsq = 0.0;
  Long64_t beam_npt = 0;
  for (Int_t k = 0; k < Int_t(m_reservoir.size()); k++)
    if (m_reservoir[k].beam_flat)
      for (Int_t s = 0; s < 18; s++) {
        Double_t d = m_reservoir[k].Total(s) - base[s];
        beam_sumsq += d * d;
        beam_npt++;
      }
  Double_t beam_sigma =
      (beam_npt > 0) ? TMath::Sqrt(beam_sumsq / beam_npt) : 0.0;

  // Count triggers over ALL reservoir events (not just classified subset) so
  // the numbers are directly comparable to the Python pipeline.
  Long64_t triggered = 0, no_trigger = 0;
  Double_t td_all[18];
  for (Int_t k = 0; k < Int_t(m_reservoir.size()); k++) {
    m_reservoir[k].Totals(td_all);
    if (FindTrigger(td_all, base, beam_sigma) >= 0)
      triggered++;
    else
      no_trigger++;
  }
  const Double_t reac_onset_gate =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRIGGER_NSIGMA * beam_sigma;
  const Double_t cf_frac =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRIGGER_CFD_FRAC;
  std::cout << "  beam reference: mean+RMS of " << nbeam
            << " pure-beam events (fitted s0,s1 & s16,s17 ellipses); "
            << Form("noise sigma=%.4f", beam_sigma) << std::endl;
  std::cout << Form("  reaction onset: gate %g-sigma = %.4f, CF fraction %g; "
                    "triggered %lld of %lld (no trigger: %lld)",
                    Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRIGGER_NSIGMA,
                    reac_onset_gate, cf_frac, Long64_t(triggered),
                    Long64_t(m_reservoir.size()), Long64_t(no_trigger))
            << std::endl;

  for (Int_t k = 0; k < Int_t(m_reservoir.size()); k++) {
    const TraceEvt &e = m_reservoir[k];
    Double_t td[18];
    e.Totals(td);
    Int_t cls = -1;
    if (e.beam_flat)
      cls = 0;
    else if (e.reac_mask & bit) {
      Double_t x = 0.0, y = 0.0;
      PlaneXY(td, reac, x, y);
      if (cut_aa && cut_aa->IsInside(x, y))
        cls = 1;
      else if (cut_an && cut_an->IsInside(x, y))
        cls = 2;
    }
    if (cls < 0)
      continue;
    // The same five variables the blind clustering uses (normed total, beam
    // at 1 per strip; strips 0/17 included).
    Double_t energy = 0.0;
    for (Int_t s = 0; s < 18; s++) {
      energy += td[s] - 1.0;
    }

    Int_t trigger_strip = FindTrigger(td, base, beam_sigma);
    Bool_t has_trig = (trigger_strip >= 0);
    // Plateau excess: sum over the sliding window tracking the trigger
    // (reacstrip+1 .. +PLATEAU_POST); 0 if no trigger; out-of-range dropped.
    const Int_t kPlateauPost =
        Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PLATEAU_POST;
    Double_t plateau = 0.0;
    if (has_trig)
      for (Int_t d = 1; d <= kPlateauPost; d++) {
        Int_t s = trigger_strip + d;
        if (s >= 0 && s < 18)
          plateau += td[s] - 1.0;
      }
    // peak-3 sum centered on the TRIGGER strip (not the argmax). 0 when there
    // is no trigger, which is the discriminator: flat beam scores 0.
    Double_t peak3 = 0.0;
    if (has_trig)
      for (Int_t s = trigger_strip - 1; s <= trigger_strip + 1; s++)
        if (s >= 0 && s < 18)
          peak3 += td[s];
    // |deviation from beam| at the TRIGGER strip plus at the end strip s17 --
    // the (a,n) signature is a rise at the trigger AND a collapse at s17.
    Double_t trigtaildev =
        has_trig ? TMath::Abs(td[trigger_strip] - base[trigger_strip]) +
                       TMath::Abs(td[17] - base[17])
                 : 0.0;
    // Slopes across the trigger: dE(reac+n) - dE(reac-n); +-n share parity
    // so the sawtooth cancels; filled only with a trigger, endpoints in range.
    Bool_t ok3 =
        has_trig && (trigger_strip - 3 >= 0) && (trigger_strip + 3 <= 17);
    Double_t reacslope3 =
        ok3 ? td[trigger_strip + 3] - td[trigger_strip - 3] : 0.0;

    // Beam-likeness of the back half (strips 8-17): RMS of trace minus beam
    // baseline; low = beam-like, high = plateau/collapse/pileup; always filled.
    Double_t beamdev = 0.0;
    Int_t n_bl = 0;
    for (Int_t s = 8; s <= 17; s++) {
      Double_t d = td[s] - base[s];
      beamdev += d * d;
      n_bl++;
    }
    beamdev = TMath::Sqrt(beamdev / Double_t(n_bl));
    Double_t v[NV] = {energy,
                      peak3,
                      plateau,
                      td[17], // raw end strip (was td[17] - 1.0)
                      Double_t(trigger_strip),
                      Double_t(e.both_mult),
                      trigtaildev,
                      reacslope3,
                      beamdev};
    // peak3 is filled even without a trigger (it scores 0, the
    // discriminator); the others skip without a trigger or an in-range window.
    Bool_t vok[NV] = {kTRUE, kTRUE,    kTRUE, kTRUE, kTRUE,
                      kTRUE, has_trig, ok3,   kTRUE};
    for (Int_t iv = 0; iv < NV; iv++)
      if (vok[iv])
        vals_raw[cls][iv].push_back(v[iv]);

    // Savitzky-Golay-smoothed trace: recompute trigger and cluster variables
    // on the SG copy; removing the sawtooth puts features on real structure.
    Double_t sgd[18];
    SavitzkyGolay(td, sgd);
    Double_t energy_sg = 0.0;
    for (Int_t s = 0; s < 18; s++)
      energy_sg += sgd[s] - 1.0;
    Double_t ex_sg[18], sm_ex_sg[18];
    for (Int_t s = 0; s < 18; s++)
      ex_sg[s] = sgd[s] - base[s];
    SmoothTrace(ex_sg, sm_ex_sg, kClusterSmoothWindow);
    Int_t reacstrip_sg = FindTrigger(sgd, base, beam_sigma);
    Bool_t has_trig_sg = (reacstrip_sg >= 0);
    // Plateau excess on the SG trace: sliding window over reacstrip_sg+1
    // .. reacstrip_sg+PLATEAU_POST; 0 when no trigger; out-of-range dropped.
    Double_t plateau_sg = 0.0;
    if (has_trig_sg)
      for (Int_t d = 1; d <= kPlateauPost; d++) {
        Int_t s = reacstrip_sg + d;
        if (s >= 0 && s < 18)
          plateau_sg += sgd[s] - 1.0;
      }
    Double_t peak3_sg = 0.0;
    if (has_trig_sg)
      for (Int_t s = reacstrip_sg - 1; s <= reacstrip_sg + 1; s++)
        if (s >= 0 && s < 18)
          peak3_sg += sgd[s];
    Double_t trigtaildev_sg =
        has_trig_sg ? TMath::Abs(sgd[reacstrip_sg] - base[reacstrip_sg]) +
                          TMath::Abs(sgd[17] - base[17])
                    : 0.0;
    Bool_t ok3_sg =
        has_trig_sg && (reacstrip_sg - 3 >= 0) && (reacstrip_sg + 3 <= 17);
    Double_t reacslope3_sg =
        ok3_sg ? sgd[reacstrip_sg + 3] - sgd[reacstrip_sg - 3] : 0.0;
    Double_t beamdev_sg = 0.0;
    Int_t n_bl_sg = 0;
    for (Int_t s = 8; s <= 17; s++) {
      Double_t d = sgd[s] - base[s];
      beamdev_sg += d * d;
      n_bl_sg++;
    }
    beamdev_sg = TMath::Sqrt(beamdev_sg / Double_t(n_bl_sg));
    Double_t v_sg[NV] = {energy_sg,
                         peak3_sg,
                         plateau_sg,
                         sgd[17],
                         Double_t(reacstrip_sg),
                         Double_t(e.both_mult),
                         trigtaildev_sg,
                         reacslope3_sg,
                         beamdev_sg};
    Bool_t vok_sg[NV] = {kTRUE, kTRUE,       kTRUE,  kTRUE, kTRUE,
                         kTRUE, has_trig_sg, ok3_sg, kTRUE};
    for (Int_t iv = 0; iv < NV; iv++)
      if (vok_sg[iv])
        vals_sg[cls][iv].push_back(v_sg[iv]);
  }

  std::cout << "cluster-var hists (reac " << reac
            << "): beam=" << vals_raw[0][0].size()
            << " (a,a')=" << vals_raw[1][0].size()
            << " (a,n)=" << vals_raw[2][0].size() << std::endl;

  std::vector<Int_t> colors = PlottingUtils::GetDefaultColors();
  // Two smoothing passes: raw traces, then Savitzky-Golay smoothed.
  const Int_t kNP = 2;
  const char *pass_label[kNP] = {"raw", "sg"};

  for (Int_t ip = 0; ip < kNP; ip++) {
    if (ip == 1 && Constants::cfg.STRIP_SUM_SCATTER_CONFIG.SKIP_SAVGOL_PLOTS)
      continue;
    std::vector<Double_t>(*vals)[NC][NV] = (ip == 0) ? &vals_raw : &vals_sg;
    for (Int_t iv = 0; iv < NV; iv++) {
      Double_t lo = 1.0e30, hi = -1.0e30;
      Double_t mean[NC] = {0.0, 0.0, 0.0};
      for (Int_t ic = 0; ic < NC; ic++) {
        for (Int_t m = 0; m < Int_t((*vals)[ic][iv].size()); m++) {
          lo = TMath::Min(lo, (*vals)[ic][iv][m]);
          hi = TMath::Max(hi, (*vals)[ic][iv][m]);
          mean[ic] += (*vals)[ic][iv][m];
        }
        if (!(*vals)[ic][iv].empty())
          mean[ic] /= Double_t((*vals)[ic][iv].size());
      }
      if (ip == 0) {
        std::cout << "  [" << pass_label[ip] << "] " << vkey[iv]
                  << ": mean beam=" << mean[0] << " (a,a')=" << mean[1]
                  << " (a,n)=" << mean[2] << " [range " << lo << ".." << hi
                  << "]" << std::endl;
      }
      Int_t nbins = 80;
      if (iv == 4) { // reaction strip: integer bins 0..17
        lo = -1.5;
        hi = 17.5;
        nbins = 19;
      } else if (iv == 5) { // both-channel multiplicity: integer bins 0..16
        lo = -0.5;
        hi = 16.5;
        nbins = 17;
      } else {
        if (!(hi > lo)) { // constant -> give it a drawable range
          lo -= 0.5;
          hi += 0.5;
        }
        Double_t pad = 0.05 * (hi - lo);
        lo -= pad;
        hi += pad;
      }

      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kTRUE);
      TString axis = Form(";%s;Counts", vtitle[iv]);
      std::vector<TH1F *> hs;
      Double_t ymax = 0.0;
      for (Int_t ic = 0; ic < NC; ic++) {
        TH1F *h = new TH1F(
            Form("h_cv_%s_%s_c%d_r%d", vkey[iv], pass_label[ip], ic, reac),
            axis, nbins, lo, hi);
        h->SetDirectory(nullptr);
        for (Int_t m = 0; m < Int_t((*vals)[ic][iv].size()); m++)
          h->Fill((*vals)[ic][iv][m]);

        PlottingUtils::ConfigureHistogram(h, colors[ic % Int_t(colors.size())],
                                          axis);
        h->SetStats(0);
        ymax = TMath::Max(ymax, h->GetMaximum());
        hs.push_back(h);
      }
      if (ymax <= 0.0)
        ymax = 1.0;
      TLegend *leg = PlottingUtils::AddLegend(0.775, 0.875, 0.70, 0.86);

      for (Int_t ic = 0; ic < NC; ic++) {
        if (ic == 0) {
          hs[0]->SetMaximum(3.0 * ymax);
          hs[0]->SetMinimum(1.0e-1);
          hs[0]->Draw("HIST");
        } else {
          hs[ic]->Draw("HIST SAME");
        }
        leg->AddEntry(hs[ic], clabel[ic], "l");
      }
      leg->Draw();

      TString sub_subdir = subdir + "/clusters_" + pass_label[ip];

      PlottingUtils::SaveFigure(
          c, Form("cluster_var_%s_%s_reac%d", vkey[iv], pass_label[ip], reac),
          sub_subdir, PlotSaveOptions::kLOG);
      for (Int_t m = 0; m < Int_t(hs.size()); m++)
        delete hs[m];
      delete c;
    }
  }
}

TString StripSumScatter::BuildFingerprint(const FileSet::GateGroups &groups) {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Double_t kReacJumpNSigma =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REAC_JUMP_NSIGMA;
  const Double_t kSmoothNSigma =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.SMOOTHNESS_NSIGMA;
  const Double_t kEndStripMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.END_STRIP_NSIGMA;

  const Int_t kGateStripX =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_STRIP_X;
  const Int_t kGateStripY =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_STRIP_Y;
  const Double_t kGateNSigmaX =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_X;
  const Double_t kGateNSigmaY =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_Y;
  const Int_t kGateBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_BINS;
  const Double_t kGateMin = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_MIN;
  const Double_t kGateMax = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_MAX;

  const Int_t kXBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.XBINS;
  const Int_t kYBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.YBINS;

  // Two parts split by the bar. Before: what decides tagging and keeping; a
  // change there refills. After: plane-only, re-projected from the reservoir.
  TString s = Form(
      "v33 reac[%d,%d] bmult[%d,%d] pileup=%.2fsig,%d noise=%.2fsig,%d "
      "jump=%.2fsig smooth=%.2fsig "
      "end=%.2fsig gate[s%d,s%d,%.2f,%.2f,%d,%.3f,%.3f]",
      kReacMin, kReacMax, Constants::cfg.STRIP_SUM_SCATTER_CONFIG.BOTH_MULT_MAX,
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.BOTH_MULT_COUNT_TO,
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PILEUP_NSIGMA,
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PILEUP_MIN_STRIPS,
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.NOISE_NSIGMA,
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.NOISE_MIN_STRIPS, kReacJumpNSigma,
      kSmoothNSigma, kEndStripMax, kGateStripX, kGateStripY, kGateNSigmaX,
      kGateNSigmaY, kGateBins, kGateMin, kGateMax);
  // The tail-shape conditions (PassesTail), in sigma, and the measured noise
  // every sigma-unit cut resolves through: a drift in the noise refills.
  {
    const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;
    s += Form(" tail[fall=%d,%.2fsig above=%d,%.2fsig rerise=%.2fsig,%.2fsig "
              "cliff=%.2f]",
              C.TAIL_FALL_FROM_STRIP, C.TAIL_RISE_NSIGMA, C.POST_ABOVE_STRIPS,
              C.POST_ABOVE_NSIGMA, C.TAIL_RETURN_NSIGMA, C.TAIL_RERISE_NSIGMA,
              C.TAIL_CLIFF_MAX_FRACTION);
    // The variation steps: the variant counts in the cache depend on them.
    if (C.CUT_VARIATION)
      s += Form(" var[%.2fsig,%.2f]", C.CUT_VARIATION_NSIGMA_STEP,
                C.CUT_VARIATION_CLIFF_STEP);
    s += " beam[";
    for (Int_t strip = 1; strip <= 17; strip++)
      s += Form("%s%.4f/%.4f", strip == 1 ? "" : ",", StripMean(strip),
                StripSigma(strip));
    s += "]";
  }
  // The gate resolves through the measured noise, so the thresholds actually
  // applied are stamped too: a drift in the noise refills.
  s += " jmin[";
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
    s += Form("%s%.4f", reac == kReacMin ? "" : ",", JumpMin(reac));
  s += "]";
  // The upstream-beam precondition changes which events are tagged; like the
  // jump gate, it resolves through measured noise, so tolerances stamped.
  s += Form(" up=%d,%.2fsig[",
            Int_t(Constants::cfg.STRIP_SUM_SCATTER_CONFIG
                      .REQUIRE_BEAM_UPSTREAM_OF_REAC),
            Constants::cfg.STRIP_SUM_SCATTER_CONFIG.BEAM_UPSTREAM_NSIGMA);
  for (Int_t strip = 1; strip < kReacMax; strip++)
    s += Form("%s%.4f", strip == 1 ? "" : ",",
              Constants::cfg.STRIP_SUM_SCATTER_CONFIG.BEAM_UPSTREAM_NSIGMA *
                  StripSigma(strip));
  s += "]";
  // Active beam gates (also keyed by cache filename, but folded in here too so
  // a mismatch never silently reuses a stale same-named cache).
  std::vector<GateSpec> gates = ActiveGates();
  for (Int_t i = 0; i < Int_t(gates.size()); i++)
    s += Form(" g[s%d,s%d]", gates[i].sx, gates[i].sy);

  // Display windows are deliberately absent: they no longer change what is
  // built, so retuning them must not invalidate the cache.
  for (Int_t i = 0; i < Int_t(groups.order.size()); i++) {
    const Int_t key = groups.order[i];
    s += Form(" %s:%lld", groups.label.find(key)->second.Data(),
              groups.chain.find(key)->second->GetEntries());
  }

  // The built quantity: build range, binning, the x range and each strip's
  // y window, so a change to the window rule re-projects.
  TString plane =
      Form("buildx[%.3f,%.3f] buildy[%.3f,%.3f] bins[%d,%d] x[%d,%d]",
           ScatterBuildRange::kXMin, ScatterBuildRange::kXMax,
           ScatterBuildRange::kYMin, ScatterBuildRange::kYMax, kXBins, kYBins,
           Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_LO,
           Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_HI);
  plane += " y";
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
    plane +=
        Form("%s%d-%d", reac == kReacMin ? "[" : ",", YLoOf(reac), YHiOf(reac));
  plane += "]";
  plane +=
      Form(" yratio=%d",
           Int_t(Constants::cfg.STRIP_SUM_SCATTER_CONFIG.Y_RATIO_TO_UPSTREAM));
  return s + " | " + plane;
}

// The tagging half of a fingerprint; a cache from before the split has no
// bar and never matches.
static TString TagPart(const TString &fingerprint) {
  const Ssiz_t bar = fingerprint.Index(" | ");
  return bar < 0 ? TString("") : TString(fingerprint(0, bar));
}

/// Per-reaction-strip y-axis bounds straight from
/// Constants::cfg.STRIP_SUM_SCATTER_CONFIG.Y_DISPLAY_RANGE (tunable per
/// dataset, per strip); strips absent from the map fall back to
/// Y_DISPLAY_MIN/Y_DISPLAY_MAX. x stays fixed (strip-independent).
void StripSumScatter::YBounds(Double_t *y_lo, Double_t *y_hi) {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    Int_t ri = reac - kReacMin;
    std::map<Int_t, std::pair<Double_t, Double_t>>::const_iterator it =
        Constants::cfg.STRIP_SUM_SCATTER_CONFIG.Y_DISPLAY_RANGE.find(reac);
    if (it != Constants::cfg.STRIP_SUM_SCATTER_CONFIG.Y_DISPLAY_RANGE.end()) {
      y_lo[ri] = it->second.first;
      y_hi[ri] = it->second.second;
    } else {
      y_lo[ri] = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.Y_DISPLAY_MIN;
      y_hi[ri] = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.Y_DISPLAY_MAX;
    }
  }
}

TString StripSumScatter::PrettyLabel(const TString &tag) {
  TString base = RemixSim::TagWithoutStrip(tag);
  base.ReplaceAll("_eres", "");
  if (base == "aa")
    return "(#alpha,#alpha')";
  if (base == "an")
    return "(#alpha,n)";
  if (base == "beam")
    return "Beam";
  return base;
}

/// Per-strip sim normalization gains: read the sim beam file and, for each
/// strip, average the (unit-gain) per-strip beam total, then set gain[s] = 1 /
/// mean[s] so every strip's sim beam lands on 1 a.u. This flattens the sim beam
/// the SAME way the per-channel data normalization flattens the experimental
/// beam
/// -- a single global factor would not, since it preserves the sim's per-strip
/// structure. Strips with no beam signal keep gain 0 (drop out like an
/// uncalibrated channel).
Bool_t StripSumScatter::SimBeamGains(Double_t *gain) {
  const Long64_t kSampleMaxPoints =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.SAMPLE_MAX_POINTS;

  for (Int_t s = 0; s < 18; s++)
    gain[s] = 0.0;
  // Reference the ERES beam file (same type as the SimOverlay populations),
  // so the eres beam lands exactly on 1; else the non-eres SIM_BEAM_FILE.
  TString file;
  std::vector<RemixSim::SimFileSpec> specs = RemixSim::BuildFileSpecs();
  for (Int_t i = 0; i < Int_t(specs.size()); i++) {
    if (!RemixSim::IsEresTag(specs[i].tag))
      std::cout << "No eres sim file for " << specs[i].tag << ", using standard"
                << std::endl;

    TString base = RemixSim::TagWithoutStrip(specs[i].tag);
    base.ReplaceAll("_eres", "");
    if (base == "beam") {
      file = RemixSim::SimRootPath(specs[i]);
      break;
    }
  }
  if (file.Length() == 0)
    file =
        Paths::DatasetDir() + "/sim_root_files/" + Constants::cfg.SIM_BEAM_FILE;
  TFile *f = IO::OpenForReading(file);
  if (!f || f->IsZombie()) {
    std::cerr << "strip-sum-scatter: cannot open sim beam file " << file
              << "; sim overlay stays in raw sim units." << std::endl;
    if (f)
      delete f;
    return kFALSE;
  }
  TTree *t = static_cast<TTree *>(f->Get("events_MeV"));
  if (!t) {
    std::cerr << "strip-sum-scatter: no events_MeV tree in sim beam file "
              << file << "; sim overlay stays in raw sim units." << std::endl;
    f->Close();
    delete f;
    return kFALSE;
  }
  RemixSim::Event e;
  if (!e.Attach(t)) {
    f->Close();
    delete f;
    return kFALSE;
  }
  Long64_t n = t->GetEntries();
  Long64_t stride = FileSet::SampleStride(n, kSampleMaxPoints);
  Double_t sum[18] = {0};
  Long64_t cnt[18] = {0};
  // Unit gains so SimTotal yields the raw per-strip beam total (IGNORE_SHORT
  // aware), which is exactly the quantity these gains will later normalize.
  Double_t unit[18];
  for (Int_t s = 0; s < 18; s++)
    unit[s] = 1.0;
  for (Long64_t j = 0; j < n; j += stride) {
    t->GetEntry(j);
    Double_t total[18];
    SimTotal(e, unit, total);
    for (Int_t s = 0; s < 18; s++)
      if (total[s] > 0.0) {
        sum[s] += total[s];
        cnt[s]++;
      }
  }
  f->Close();
  delete f;
  Int_t n_set = 0;
  for (Int_t s = 0; s < 18; s++) {
    if (cnt[s] > 0 && sum[s] > 0.0) {
      gain[s] = 1.0 / (sum[s] / Double_t(cnt[s]));
      n_set++;
    }
  }
  if (n_set == 0)
    return kFALSE;
  std::cout << "strip-sum-scatter: sim per-strip beam normalization to 1 a.u. ("
            << n_set << " strips)." << std::endl;
  return kTRUE;
}

void StripSumScatter::SimTotal(const RemixSim::Event &e, const Double_t *gain,
                               Double_t *total) {
  for (Int_t s = 0; s < 18; s++)
    total[s] = gain[s] * e.Total(s);
  if (Constants::cfg.IGNORE_SHORT_STRIPS)
    for (Int_t s = 1; s <= 16; s++)
      total[s] = gain[s] * ((s % 2) != 0 ? e.Left(s) : e.Right(s));
}

TGraph *StripSumScatter::SimPopScatter(const TString &file, Int_t reac,
                                       const Double_t *gain,
                                       Long64_t max_points) {
  const Int_t kXLo = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_HI;

  TFile *f = IO::OpenForReading(file);
  if (!f || f->IsZombie()) {
    if (f)
      delete f;
    return nullptr;
  }
  TTree *t = static_cast<TTree *>(f->Get("events_MeV"));
  if (!t) {
    std::cerr << "  no events_MeV tree in " << file << std::endl;
    f->Close();
    delete f;
    return nullptr;
  }
  Int_t y_lo = reac + 1;
  Int_t y_hi = TMath::Min(reac + 6, 17);
  RemixSim::Event e;
  if (!e.Attach(t)) {
    f->Close();
    delete f;
    return nullptr;
  }
  Long64_t n = t->GetEntries();
  Long64_t stride = FileSet::SampleStride(n, max_points);
  TGraph *g = new TGraph();
  Long64_t k = 0;
  for (Long64_t j = 0; j < n; j += stride) {
    t->GetEntry(j);
    Double_t total[18];
    SimTotal(e, gain, total);
    Double_t x = SumRange(total, kXLo, kXHi);
    Double_t y = SumRange(total, y_lo, y_hi);
    if (x > 0.0)
      g->SetPoint(k++, x, y);
  }
  g->Set(k);
  f->Close();
  delete f;
  return g;
}

/// Up to max_traces per-strip dE-profile traces, stride-sampled across one sim
/// file. Sim energies are arbitrary-unit floats; total[s] is the per-strip
/// normalized gain[s]*(left[s] + right[s]) so the traces share the data's axis
/// (and the sim beam is flat at NORM, like the data).
std::vector<TGraph *> StripSumScatter::SimPopTraces(const TString &file,
                                                    const Double_t *gain,
                                                    Long64_t max_traces) {
  std::vector<TGraph *> traces;
  TFile *f = IO::OpenForReading(file);
  if (!f || f->IsZombie()) {
    if (f)
      delete f;
    return traces;
  }
  TTree *t = static_cast<TTree *>(f->Get("events_MeV"));
  if (!t) {
    f->Close();
    delete f;
    return traces;
  }
  RemixSim::Event e;
  if (!e.Attach(t)) {
    f->Close();
    delete f;
    return traces;
  }
  Long64_t n = t->GetEntries();
  Long64_t stride = FileSet::SampleStride(n, max_traces);
  for (Long64_t j = 0; j < n && Int_t(traces.size()) < max_traces;
       j += stride) {
    t->GetEntry(j);
    Double_t total[18];
    SimTotal(e, gain, total);
    traces.push_back(EventsSummary::BuildTraceFromTotals(total));
  }
  f->Close();
  delete f;
  return traces;
}

/// Per reaction strip, overlay TRACES_PER_CLASS sampled per-strip traces of
/// each sim population in the experimental DrawRegionTraces style (beam grey,
/// (a,a') azure, (a,n) red). The beam reference is the same for every strip.
/// Sampled fresh each run (40 traces/file is trivial).
void StripSumScatter::SimTraceOverlay() {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kTracesPerRegion =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRACES_PER_CLASS;

  std::vector<RemixSim::SimFileSpec> specs = RemixSim::BuildFileSpecs();
  if (specs.empty())
    return;
  // One file per class and strip, the _eres twin preferred where it exists
  // (alphabetical order would otherwise leave the plain file in the map).
  std::map<Int_t, TString> aa_file, an_file; // reaction strip -> sim file
  std::map<Int_t, Bool_t> aa_eres, an_eres;
  std::vector<TString> beam_files;
  Bool_t beam_eres = kFALSE;
  for (Int_t i = 0; i < Int_t(specs.size()); i++) {
    TString base = RemixSim::TagWithoutStrip(specs[i].tag);
    base.ReplaceAll("_eres", "");
    TString file = RemixSim::SimRootPath(specs[i]);
    Int_t strip = RemixSim::ReactionStripOf(specs[i].tag);
    const Bool_t eres = RemixSim::IsEresTag(specs[i].tag);
    if (base == "beam") {
      if (beam_eres && !eres)
        continue;
      if (eres && !beam_eres)
        beam_files.clear();
      beam_files.push_back(file);
      beam_eres = beam_eres || eres;
    } else if (strip >= kReacMin && strip <= kReacMax) {
      std::map<Int_t, TString> &files = base == "aa" ? aa_file : an_file;
      std::map<Int_t, Bool_t> &have_eres = base == "aa" ? aa_eres : an_eres;
      if (base != "aa" && base != "an")
        continue;
      if (files.count(strip) && (have_eres[strip] || !eres))
        continue;
      files[strip] = file;
      have_eres[strip] = eres;
    }
  }

  Double_t gain[18];
  if (!SimBeamGains(gain))
    for (Int_t s = 0; s < 18; s++)
      gain[s] = 1.0;

  std::vector<TGraph *> beam_traces;
  for (Int_t i = 0; i < Int_t(beam_files.size()) &&
                    Int_t(beam_traces.size()) < kTracesPerRegion;
       i++) {
    std::vector<TGraph *> t = SimPopTraces(
        beam_files[i], gain, kTracesPerRegion - Int_t(beam_traces.size()));
    for (Int_t k = 0; k < Int_t(t.size()); k++)
      beam_traces.push_back(t[k]);
  }

  for (Int_t r = kReacMin; r <= kReacMax; r++) {
    std::vector<TGraph *> aa_traces, an_traces;
    if (aa_file.find(r) != aa_file.end())
      aa_traces = SimPopTraces(aa_file[r], gain, kTracesPerRegion);
    if (an_file.find(r) != an_file.end())
      an_traces = SimPopTraces(an_file[r], gain, kTracesPerRegion);
    if (aa_traces.empty() && an_traces.empty())
      continue;
    DrawRegionTraces(Form("sim_region_traces_reac%d", r), "sim_scatter",
                     beam_traces, aa_traces, an_traces, 0.6, 1.6,
                     "#DeltaE [a.u.]");
    for (Int_t i = 0; i < Int_t(aa_traces.size()); i++)
      delete aa_traces[i];
    for (Int_t i = 0; i < Int_t(an_traces.size()); i++)
      delete an_traces[i];
  }
  for (Int_t i = 0; i < Int_t(beam_traces.size()); i++)
    delete beam_traces[i];
}

/// Fingerprint of the sim inputs + window geometry: each eres file's size+mtime
/// (cheap, no open) plus the reaction-strip range and x window. Regenerating
/// the sim (new mtimes) or changing the windows invalidates the cached overlay.
TString StripSumScatter::SimFingerprint(
    const std::vector<RemixSim::SimFileSpec> &specs) {
  // v2: per-strip norm via SimBeamGains(); the eres gain file is stamped by
  // the per-spec loop; v3: norm hardcoded to 1 a.u. (NORM_MUSIC_MEV gone).
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kXLo = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_HI;

  // v4: one population per class and strip (the eres twin preferred), so
  // caches that held every population twice are rebuilt.
  TString s = Form("v4 reac[%d,%d] x[%d,%d]", kReacMin, kReacMax, kXLo, kXHi);
  for (Int_t i = 0; i < Int_t(specs.size()); i++) {
    TString f = RemixSim::SimRootPath(specs[i]);
    Long_t id = 0, flags = 0, mtime = 0;
    Long64_t size = -1;
    if (gSystem->GetPathInfo(f, &id, &size, &flags, &mtime) != 0) {
      size = -1;
      mtime = 0;
    }
    s += Form(" %s:%lld:%ld", specs[i].tag.Data(), size, mtime);
  }
  return s;
}

/// Reload cached sim scatter graphs (grouped by reaction strip; each graph's
/// title holds its population label) if the fingerprint matches. Caller owns
/// the returned graphs.
Bool_t StripSumScatter::LoadSimCache(
    const TString &fp, std::map<Int_t, std::vector<TGraph *>> &by_strip) {
  TString full = IO::GetRootFilesBaseDir() + TString("/") +
                 "StripSumScatter_simcache.root";
  if (gSystem->AccessPathName(full))
    return kFALSE;
  TFile *f = IO::OpenForReading("StripSumScatter_simcache.root");
  if (!f || f->IsZombie()) {
    if (f)
      delete f;
    return kFALSE;
  }
  TNamed *cfp = static_cast<TNamed *>(f->Get("sim_fingerprint"));
  if (!cfp || fp != cfp->GetTitle()) {
    f->Close();
    delete f;
    return kFALSE;
  }
  TIter next(f->GetListOfKeys());
  TKey *key;
  while ((key = static_cast<TKey *>(next()))) {
    TString name = key->GetName();
    if (!name.BeginsWith("simg_r"))
      continue;
    TString rest = name(6, name.Length() - 6); // after "simg_r": <strip>_p<idx>
    Int_t us = rest.Index("_p");
    if (us < 0)
      continue;
    Int_t r = TString(rest(0, us)).Atoi();
    TGraph *g = static_cast<TGraph *>(f->Get(name));
    if (!g)
      continue;
    by_strip[r].push_back(static_cast<TGraph *>(g->Clone()));
  }
  f->Close();
  delete f;
  return kTRUE;
}

void StripSumScatter::WriteSimCache(
    const TString &fp, const std::map<Int_t, std::vector<TGraph *>> &by_strip) {
  TFile *out = IO::OpenForWriting("StripSumScatter_simcache.root", "RECREATE");
  if (!out || out->IsZombie()) {
    if (out)
      delete out;
    return;
  }
  out->cd();
  TNamed cfp("sim_fingerprint", fp.Data());
  cfp.Write();
  std::map<Int_t, std::vector<TGraph *>>::const_iterator it;
  for (it = by_strip.begin(); it != by_strip.end(); ++it)
    for (Int_t i = 0; i < Int_t(it->second.size()); i++)
      it->second[i]->Write(Form("simg_r%d_p%d", it->first, Int_t(i)));
  out->Close();
  delete out;
}

/// Sim-only comparison plots: one per reaction strip, each sim population a
/// coloured+labelled point cloud on the same axes as that strip's data scatter,
/// for side-by-side comparison with the data PID scatters. Beam (no reaction
/// strip) overlays on every strip. The scatter graphs are fingerprint-cached
/// (sim file sizes/mtimes + window geometry), so re-runs reload them instead of
/// rescanning the sim files.
void StripSumScatter::SimOverlay() {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kXLo = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_HI;

  std::vector<RemixSim::SimFileSpec> specs = RemixSim::BuildFileSpecs();
  if (specs.empty()) {
    std::cerr
        << "strip-sum-scatter: no sim control files; skipping sim overlay."
        << std::endl;
    return;
  }
  TString fp = SimFingerprint(specs);

  std::map<Int_t, std::vector<TGraph *>>
      by_strip; // strip -> graphs (title=label)
  Bool_t loaded = LoadSimCache(fp, by_strip);
  if (!loaded) {
    // One population per class+strip; each has two control files, plain and
    // its _eres twin (same label); the eres file wins where it exists.
    std::map<std::pair<TString, Int_t>, std::pair<SimPop, Bool_t>> chosen;
    for (Int_t i = 0; i < Int_t(specs.size()); i++) {
      TString base = RemixSim::TagWithoutStrip(specs[i].tag);
      base.ReplaceAll("_eres", "");
      const Int_t strip = RemixSim::ReactionStripOf(specs[i].tag);
      const Bool_t eres = RemixSim::IsEresTag(specs[i].tag);
      std::pair<TString, Int_t> key(base, strip);
      if (chosen.count(key) && (chosen[key].second || !eres))
        continue;
      SimPop p;
      p.file = RemixSim::SimRootPath(specs[i]);
      p.label = PrettyLabel(specs[i].tag);
      chosen[key] = std::make_pair(p, eres);
    }
    std::map<Int_t, std::vector<SimPop>> reacted;
    std::vector<SimPop> refs;
    for (std::map<std::pair<TString, Int_t>,
                  std::pair<SimPop, Bool_t>>::const_iterator it =
             chosen.begin();
         it != chosen.end(); ++it) {
      if (!it->second.second)
        std::cout << "No eres sim file for " << it->first.first
                  << (it->first.second >= 0 ? Form("_s%d", it->first.second)
                                            : "")
                  << ", using standard" << std::endl;
      if (it->first.second < 0)
        refs.push_back(it->second.first);
      else
        reacted[it->first.second].push_back(it->second.first);
    }
    Double_t gain[18];
    if (!SimBeamGains(gain))
      for (Int_t s = 0; s < 18; s++)
        gain[s] = 1.0;
    const Long64_t kSimMaxPoints = 25000;
    for (Int_t r = kReacMin; r <= kReacMax; r++) {
      std::vector<SimPop> group = reacted[r];
      for (Int_t i = 0; i < Int_t(refs.size()); i++)
        group.push_back(refs[i]);
      for (Int_t i = 0; i < Int_t(group.size()); i++) {
        TGraph *g = SimPopScatter(group[i].file, r, gain, kSimMaxPoints);
        if (!g || g->GetN() == 0) {
          if (g)
            delete g;
          continue;
        }
        g->SetTitle(group[i].label);
        by_strip[r].push_back(g);
      }
    }
    Int_t n_graphs = 0;
    std::map<Int_t, std::vector<TGraph *>>::const_iterator cit;
    for (cit = by_strip.begin(); cit != by_strip.end(); ++cit)
      n_graphs += Int_t(cit->second.size());
    if (n_graphs == 0) {
      std::cerr << "strip-sum-scatter: no sim data found (regenerate "
                   "sim_root_files); skipping sim overlay."
                << std::endl;
      return;
    }
    WriteSimCache(fp, by_strip);
    std::cout << "strip-sum-scatter: built + cached sim overlay (" << n_graphs
              << " population graphs)." << std::endl;
  } else {
    std::cout
        << "strip-sum-scatter: loaded cached sim overlay (fingerprint match)."
        << std::endl;
  }

  std::map<Int_t, std::vector<TGraph *>>::iterator it;
  for (it = by_strip.begin(); it != by_strip.end(); ++it) {
    Int_t r = it->first;
    std::map<Int_t, TH2F *>::const_iterator sit = m_scatter.find(r);
    if (sit == m_scatter.end())
      continue;
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    TH2F *ref = sit->second;
    TH2F *frame =
        new TH2F(Form("sim_frame_r%d", r), "", 10, ref->GetXaxis()->GetXmin(),
                 ref->GetXaxis()->GetXmax(), 10, ref->GetYaxis()->GetXmin(),
                 ref->GetYaxis()->GetXmax());
    frame->SetStats(0);
    frame->GetXaxis()->SetTitle(ref->GetXaxis()->GetTitle());
    frame->GetYaxis()->SetTitle(ref->GetYaxis()->GetTitle());
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    c->SetLeftMargin(0.18);
    frame->Draw();
    // Match the experimental region-traces legend placement (top-right).
    TLegend *leg = PlottingUtils::AddLegend(0.725, 0.875, 0.70, 0.86);
    for (Int_t i = 0; i < Int_t(it->second.size()); i++) {
      TGraph *g = it->second[i];
      // Match the experimental region-trace colours (DrawRegionTraces): beam
      // grey, (a,a') azure, (a,n) red -- keyed off the population label.
      TString lab = g->GetTitle();
      Int_t color = kBlack;
      if (lab == "Beam")
        color = kGray + 2;
      else if (lab == "(#alpha,#alpha')")
        color = kAzure + 2;
      else if (lab == "(#alpha,n)")
        color = kRed + 1;
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.3);
      g->SetMarkerColorAlpha(color, 0.35);
      g->SetLineColor(color);
      g->Draw("P SAME");
      leg->AddEntry(g, g->GetTitle(), "p");
    }
    leg->Draw();
    PlottingUtils::SaveFigure(c,
                              Form("sim_normsumE_reac%d_s%d_%d_vs_s%d_%d", r,
                                   YLoOf(r), YHiOf(r), kXLo, kXHi),
                              "sim_scatter", PlotSaveOptions::kLINEAR);
    delete leg;
    delete c;
    delete frame;
  }

  std::map<Int_t, std::vector<TGraph *>>::iterator dit;
  for (dit = by_strip.begin(); dit != by_strip.end(); ++dit)
    for (Int_t i = 0; i < Int_t(dit->second.size()); i++)
      delete dit->second[i];
}

Bool_t StripSumScatter::TryLoadCache(const TString &cacheName,
                                     const TString &fingerprint) {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;

  TString cache_full = IO::GetRootFilesBaseDir() + TString("/") + cacheName;
  if (gSystem->AccessPathName(cache_full)) {
    std::cout << "strip-sum-scatter: no cache file found; will rebuild."
              << std::endl;
    return kFALSE;
  }

  TFile *cf = IO::OpenForReading(cacheName);
  if (!cf || cf->IsZombie()) {
    if (cf)
      delete cf;
    std::cout << "strip-sum-scatter: cache file unreadable; rebuilding."
              << std::endl;
    return kFALSE;
  }

  TNamed *fp = static_cast<TNamed *>(cf->Get("fingerprint"));
  const Bool_t exact = fp && fingerprint == fp->GetTitle();
  // Same tagging, different plane: the reservoir has every tagged event, so
  // the scatters are rebuilt from it below and the cache rewritten.
  const Bool_t reproject = !exact && fp && !TagPart(fingerprint).IsNull() &&
                           TagPart(fingerprint) == TagPart(fp->GetTitle());
  if (!exact && !reproject) {
    std::cout << "strip-sum-scatter: cache present but stale; rebuilding."
              << std::endl;
    std::cout << "  cached: " << (fp ? fp->GetTitle() : "(none)") << std::endl;
    std::cout << "  wanted: " << fingerprint << std::endl;
    cf->Close();
    delete cf;
    return kFALSE;
  }

  Bool_t ok = kTRUE;
  for (Int_t reac = kReacMin; reac <= kReacMax && ok && exact; reac++) {
    TH2F *h = static_cast<TH2F *>(cf->Get(Form("scatter_r%d", reac)));
    if (!h) {
      ok = kFALSE;
      break;
    }
    TH2F *hc = static_cast<TH2F *>(h->Clone());
    hc->SetDirectory(nullptr);
    m_scatter[reac] = hc;
  }

  // The normalization counts ride along, so a re-projected cache keeps them.
  if (TParameter<Long64_t> *p =
          dynamic_cast<TParameter<Long64_t> *>(cf->Get("n_seen")))
    m_nSeen = p->GetVal();
  if (TParameter<Long64_t> *p =
          dynamic_cast<TParameter<Long64_t> *>(cf->Get("n_normed")))
    m_nNormed = p->GetVal();
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    if (TParameter<Long64_t> *p = dynamic_cast<TParameter<Long64_t> *>(
            cf->Get(Form("n_normed_r%d", reac))))
      m_normedAt[ReacIndex(reac)] = p->GetVal();
    if (TParameter<Long64_t> *p = dynamic_cast<TParameter<Long64_t> *>(
            cf->Get(Form("n_tagged_r%d", reac))))
      m_tagged[ReacIndex(reac)] = p->GetVal();
  }
  // The cut-variation counts, by the variant names the cache lists; a cache
  // written without them (or with the variation off) carries none.
  m_variantNames.clear();
  m_taggedVar.clear();
  if (TNamed *vn = dynamic_cast<TNamed *>(cf->Get("cut_variants"))) {
    TString names = vn->GetTitle();
    Int_t from = 0;
    TString tok;
    while (names.Tokenize(tok, from, ",")) {
      if (tok.IsNull())
        continue;
      m_variantNames.push_back(tok);
      std::vector<Long64_t> counts(m_tagged.size(), 0);
      for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
        if (TParameter<Long64_t> *p = dynamic_cast<TParameter<Long64_t> *>(
                cf->Get(Form("n_tagged_r%d_%s", reac, tok.Data()))))
          counts[ReacIndex(reac)] = p->GetVal();
      m_taggedVar.push_back(counts);
    }
  }

  TTree *tt = static_cast<TTree *>(cf->Get("traces"));
  if (reproject && !tt)
    ok = kFALSE;
  if (ok && tt) {
    TraceEvt e;
    tt->SetBranchAddress("leftdE", e.leftdE);
    tt->SetBranchAddress("rightdE", e.rightdE);
    tt->SetBranchAddress("strip0dE", &e.strip0dE);
    tt->SetBranchAddress("strip17dE", &e.strip17dE);
    tt->SetBranchAddress("leftdE_adc", e.leftdE_adc);
    tt->SetBranchAddress("rightdE_adc", e.rightdE_adc);
    tt->SetBranchAddress("strip0_adc", &e.strip0_adc);
    tt->SetBranchAddress("strip17_adc", &e.strip17_adc);
    tt->SetBranchAddress("reac_mask", &e.reac_mask);
    tt->SetBranchAddress("beam_flat", &e.beam_flat);
    tt->SetBranchAddress("both_mult", &e.both_mult);
    e.seed_ts = 0;
    if (tt->GetBranch("seed_ts"))
      tt->SetBranchAddress("seed_ts", &e.seed_ts);
    Long64_t nt = tt->GetEntries();
    m_reservoir.reserve(nt);
    for (Long64_t j = 0; j < nt; j++) {
      tt->GetEntry(j);
      m_reservoir.push_back(e);
    }
  }

  cf->Close();
  delete cf;

  if (!ok) {
    std::cout << "strip-sum-scatter: cache partially corrupt; rebuilding."
              << std::endl;
    m_reservoir.clear();
    return kFALSE;
  }
  if (reproject) {
    std::cout << "strip-sum-scatter: cache tagging matches but the plane "
                 "changed; re-projecting "
              << m_reservoir.size() << " reservoir events." << std::endl;
    AllocateScatters();
    ReprojectFromReservoir();
    WriteCache(cacheName, fingerprint);
  } else {
    std::cout << "strip-sum-scatter: loaded cached scatters + "
              << m_reservoir.size() << " reservoir events (fingerprint match)."
              << std::endl;
  }
  return kTRUE;
}

void StripSumScatter::WriteCache(const TString &cacheName,
                                 const TString &fingerprint) {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;

  TFile *out = IO::OpenForWriting(cacheName, "RECREATE");
  if (!out || out->IsZombie()) {
    if (out)
      delete out;
    return;
  }
  out->cd();
  TNamed fp("fingerprint", fingerprint.Data());
  fp.Write();
  // Normalization counts travel with the scatters they describe, so a cross
  // section never has to re-derive them from a different pass over the data.
  TParameter<Long64_t>("n_seen", m_nSeen).Write();
  TParameter<Long64_t>("n_normed", m_nNormed).Write();
  // The beam level and noise every sigma-unit cut was resolved against.
  for (Int_t s = 0; s < 18; s++) {
    TParameter<Double_t>(Form("strip_mean_s%d", s), s_stripMean[s]).Write();
    TParameter<Double_t>(Form("strip_sigma_s%d", s), s_stripSigma[s]).Write();
  }
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
    TParameter<Long64_t>(Form("n_normed_r%d", reac),
                         m_normedAt[ReacIndex(reac)])
        .Write();
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
    TParameter<Long64_t>(Form("n_tagged_r%d", reac), m_tagged[ReacIndex(reac)])
        .Write();
  // The cut-variation counts: the variant names as one list, then one count
  // per variant and strip.
  {
    TString names;
    for (size_t v = 0; v < m_variantNames.size(); v++)
      names += (v ? "," : "") + m_variantNames[v];
    TNamed("cut_variants", names.Data()).Write();
    for (size_t v = 0; v < m_variantNames.size(); v++)
      for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
        TParameter<Long64_t>(
            Form("n_tagged_r%d_%s", reac, m_variantNames[v].Data()),
            m_taggedVar[v][ReacIndex(reac)])
            .Write();
  }
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
    m_scatter[reac]->Write(Form("scatter_r%d", reac));

  TTree *tt = new TTree("traces", "strip-sum trace reservoir");
  TraceEvt e;
  tt->Branch("leftdE", e.leftdE, "leftdE[16]/F");
  tt->Branch("rightdE", e.rightdE, "rightdE[16]/F");
  tt->Branch("strip0dE", &e.strip0dE, "strip0dE/F");
  tt->Branch("strip17dE", &e.strip17dE, "strip17dE/F");
  tt->Branch("leftdE_adc", e.leftdE_adc, "leftdE_adc[16]/F");
  tt->Branch("rightdE_adc", e.rightdE_adc, "rightdE_adc[16]/F");
  tt->Branch("strip0_adc", &e.strip0_adc, "strip0_adc/F");
  tt->Branch("strip17_adc", &e.strip17_adc, "strip17_adc/F");
  tt->Branch("reac_mask", &e.reac_mask, "reac_mask/i");
  tt->Branch("beam_flat", &e.beam_flat, "beam_flat/O");
  tt->Branch("both_mult", &e.both_mult, "both_mult/I");
  tt->Branch("seed_ts", &e.seed_ts, "seed_ts/l");
  for (Int_t k = 0; k < Int_t(m_reservoir.size()); k++) {
    e = m_reservoir[k];
    tt->Fill();
  }
  tt->Write();
  out->Close();
  delete out;
  std::cout << "strip-sum-scatter: wrote cache " << cacheName << std::endl;
}

/// One run's beam ellipses and series gates. Split out so runs can be fitted in
/// parallel: each call touches only its own chain and returns its own result,
/// with no shared state to guard.
SingleRunFitResult
StripSumScatter::FitRunGates(Int_t key, const TString &label, TChain *chain,
                             const std::vector<GateSpec> &activeGates) {
  (void)key;
  const TString run = label; // the group's name in every message
  SingleRunFitResult res;
  if (!chain || chain->GetEntries() == 0)
    return res;

  // --- Beam classification ellipses (no gating) ---
  BeamEllipses be;
  be.ok = kFALSE;
  {
    std::vector<GateSpec> emptyPrior;
    std::vector<BeamFit2D> emptyGates;
    const TString tag = label;
    const TString subdir = "strip_sum_scatter/" + label;
    Int_t ent_sx = 0, ent_sy = 1;
    const Char_t *ent_tag = "s0/s1";
    if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PURE_BEAM_GATE ==
        StripSumScatterConfig::PURE_BEAM_GATE_S1_S2) {
      ent_sx = 1;
      ent_sy = 2;
      ent_tag = "s1/s2";
    }
    BeamFit2D ent_ell = FindBeamGate(chain, ent_sx, ent_sy, emptyPrior,
                                     emptyGates, tag, subdir);
    if (ent_sx == 0)
      be.s0_s1 = ent_ell;
    else
      be.s1_s2 = ent_ell;
    if (ent_ell.ok) {
      std::lock_guard<std::mutex> lk(g_log_mutex);
      std::cout << "  " << run << " beam ellipse " << ent_tag << ": mu=("
                << ent_ell.mu_x << "," << ent_ell.mu_y << ")" << std::endl;
    } else {
      std::lock_guard<std::mutex> lk(g_log_mutex);
      std::cerr << "  " << run << " beam ellipse " << ent_tag
                << " failed; skipping run" << std::endl;
      return res;
    }
    if (Constants::cfg.IGNORE_STRIP_17) {
      be.use_s15_s16 = kTRUE;
      be.s15_s16 =
          FindBeamGate(chain, 15, 16, emptyPrior, emptyGates, tag, subdir);
      if (be.s15_s16.ok) {
        std::lock_guard<std::mutex> lk(g_log_mutex);
        std::cout << "  " << run << " beam ellipse s15/s16: mu=("
                  << be.s15_s16.mu_x << "," << be.s15_s16.mu_y << ")"
                  << std::endl;
      } else {
        std::lock_guard<std::mutex> lk(g_log_mutex);
        std::cerr << "  " << run << " beam ellipse s15/s16 failed; skipping run"
                  << std::endl;
        return res;
      }
    } else {
      be.use_s15_s16 = kFALSE;
      be.s16_s17 =
          FindBeamGate(chain, 16, 17, emptyPrior, emptyGates, tag, subdir);
      if (be.s16_s17.ok) {
        std::lock_guard<std::mutex> lk(g_log_mutex);
        std::cout << "  " << run << " beam ellipse s16/s17: mu=("
                  << be.s16_s17.mu_x << "," << be.s16_s17.mu_y << ")"
                  << std::endl;
      } else {
        std::lock_guard<std::mutex> lk(g_log_mutex);
        std::cerr << "  " << run << " beam ellipse s16/s17 failed; skipping run"
                  << std::endl;
        return res;
      }
    }
    be.ok = kTRUE;
  }
  res.pure_beam = be;

  // --- Scatter filter gates (series gating) ---
  std::vector<BeamFit2D> runGates;
  std::vector<GateSpec> priorSpecs;
  Bool_t allOk = kTRUE;
  for (Int_t gi = 0; gi < Int_t(activeGates.size()); gi++) {
    BeamFit2D g =
        FindBeamGate(chain, activeGates[gi].sx, activeGates[gi].sy, priorSpecs,
                     runGates, label, "strip_sum_scatter/" + label);
    if (g.ok) {
      std::lock_guard<std::mutex> lk(g_log_mutex);
      std::cout << "  " << run << " beam gate s" << activeGates[gi].sx << "/s"
                << activeGates[gi].sy << ": mu=(" << g.mu_x << "," << g.mu_y
                << ")" << std::endl;
    } else {
      std::lock_guard<std::mutex> lk(g_log_mutex);
      std::cerr << "  " << run << " beam gate s" << activeGates[gi].sx << "/s"
                << activeGates[gi].sy << " failed; skipping run" << std::endl;
      allOk = kFALSE;
    }
    runGates.push_back(g);
    priorSpecs.push_back(activeGates[gi]);
  }
  res.series_gates = runGates;
  res.ok = allOk;
  return res;
}

/// One run's scatter fill. Each call builds PRIVATE scatter histograms and its
/// own reservoir slice, so the runs never touch shared state; the caller merges
/// them in run order, which makes the threaded result identical to sequential.
SingleRunFillResult
StripSumScatter::FillRunScatters(Int_t key, const TString &label, TChain *chain,
                                 const std::vector<GateSpec> &activeGates,
                                 const std::vector<BeamFit2D> &runGates,
                                 const BeamEllipses &runBeam) {
  const Int_t run = key; // unique per task: names the private histograms
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kXLo = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_HI;
  const Int_t kXBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.XBINS;
  const Int_t kYBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.YBINS;
  const Int_t kBeamReservoirCap =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRACES_PER_CLASS * 10;
  const Int_t nReacStrips = kReacMax - kReacMin + 1;

  SingleRunFillResult res;
  res.scatters.assign(nReacStrips, nullptr);
  if (!chain)
    return res;
  // Private, directory-less clones over the same fixed build range as the
  // merged ones.
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    TH2F *h =
        new TH2F(Form("scatter_r%d_run%d", reac, run), "", kXBins,
                 ScatterBuildRange::kXMin, ScatterBuildRange::kXMax, kYBins,
                 ScatterBuildRange::kYMin, ScatterBuildRange::kYMax);
    h->SetDirectory(nullptr);
    res.scatters[ReacIndex(reac)] = h;
  }
  Long64_t totalGated = 0, totalSeen = 0, totalNormed = 0;
  res.tagged.assign(nReacStrips, 0);
  res.normed_at.assign(nReacStrips, 0);
  res.cut_counts.assign(nReacStrips * kNTagCuts, 0);
  res.pre_counts.assign(kNPreCuts, 0);
  // The cut-variation set, fixed for the run: each event past the event-level
  // cuts is tagged once per variant as well.
  const std::vector<std::pair<TString, TagThresholds>> variants =
      ThresholdVariants();
  res.tagged_var.assign(variants.size(), std::vector<Long64_t>(nReacStrips, 0));
  Int_t nBeamKept = 0;
  EnergyView ev;
  ev.Attach(chain);
  EnableEventBranches(chain);
  ULong64_t seed_ts_in = 0;
  if (chain->GetBranch("SeedTs"))
    chain->SetBranchAddress("SeedTs", &seed_ts_in);
  Long64_t n = chain->GetEntries();
  Int_t nReac = kReacMax - kReacMin + 1;
  {
    std::lock_guard<std::mutex> lk(g_log_mutex);
    std::cout << "  " << label << ": filling " << nReac
              << " reaction-strip scatters over " << n << " events..."
              << std::endl;
  }

  for (Long64_t j = 0; j < n; j++) {
    chain->GetEntry(j);
    ev.Decode();
    totalSeen++;

    Bool_t passesAll = kTRUE;
    for (Int_t gi = 0; gi < Int_t(activeGates.size()); gi++)
      if (!PassesGate(runGates[gi], ev, activeGates[gi].sx,
                      activeGates[gi].sy)) {
        passesAll = kFALSE;
        break;
      }
    if (!passesAll) {
      res.pre_counts[kPreGate]++;
      continue;
    }
    if (IsPileup(ev)) { // reject overlapping-beam pileup
      res.pre_counts[kPrePileup]++;
      continue;
    }
    if (IsNoise(ev)) {
      res.pre_counts[kPreNoise]++;
      continue;
    }
    // Both-ends multiplicity: counted on raw ADC so it sees the short end
    // even when IGNORE_SHORT_STRIPS zeroes it in the decode.
    if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.BOTH_MULT_MAX >= 0) {
      const Int_t hi = TMath::Min(
          16, Constants::cfg.STRIP_SUM_SCATTER_CONFIG.BOTH_MULT_COUNT_TO);
      Int_t nboth = 0;
      for (Int_t s = 1; s <= hi; s++)
        if (ev.leftdE_adc[s - 1] > 0 && ev.rightdE_adc[s - 1] > 0)
          nboth++;
      if (nboth > Constants::cfg.STRIP_SUM_SCATTER_CONFIG.BOTH_MULT_MAX) {
        res.pre_counts[kPreBothMult]++;
        continue;
      }
    }
    res.pre_counts[kPrePass]++;
    // The last pre-reaction count: here, not at `seen`, is what makes the
    // tag-count ratio a cross section: same gate + quality efficiencies.
    totalNormed++;
    // Per-strip denominator: beam counts toward `reac` only if it met the
    // conditions a reaction there must meet: efficiencies cancel in the ratio.
    if (AllStripsFired(ev)) {
      for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
        if (BeamUpstreamOf(ev, reac))
          res.normed_at[ReacIndex(reac)]++;
    }

    UInt_t mask = 0;
    Double_t totals[18];
    ev.Totals(totals);
    for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
      const TagCut why = RejectReason(ev, reac);
      res.cut_counts[ReacIndex(reac) * kNTagCuts + why]++;
      // The variants: a strip that did not fire fails every set alike, so
      // only events past that are worth re-tagging.
      if (why != kCutAllStrips)
        for (size_t v = 0; v < variants.size(); v++)
          if (RejectReason(ev, reac, variants[v].second) == kTagPass)
            res.tagged_var[v][ReacIndex(reac)]++;
      if (why != kTagPass)
        continue;
      mask |= (1u << ReacIndex(reac));
      res.tagged[ReacIndex(reac)]++;
      Double_t x = 0.0, y = 0.0;
      PlaneXY(totals, reac, x, y);
      res.scatters[ReacIndex(reac)]->Fill(x, y);
    }

    // Keep reaction-passing events for traces; cap pure-beam (mutually
    // exclusive, no reaction jump); a per-task bound: the merge re-caps it.
    Bool_t beam = (mask == 0) && IsPureBeam(ev, runBeam);
    if (mask == 0 && !(beam && nBeamKept < kBeamReservoirCap))
      continue;
    if (beam)
      nBeamKept++;

    TraceEvt e;
    for (Int_t k = 0; k < 16; k++) {
      e.leftdE[k] = Float_t(ev.left[k]);
      e.rightdE[k] = Float_t(ev.right[k]);
      e.leftdE_adc[k] = Float_t(ev.leftdE_adc[k]);
      e.rightdE_adc[k] = Float_t(ev.rightdE_adc[k]);
    }
    e.strip0dE = Float_t(ev.strip0);
    e.strip17dE = Float_t(ev.strip17);
    e.strip0_adc = Float_t(ev.strip0_adc);
    e.strip17_adc = Float_t(ev.strip17_adc);
    // Mirror IGNORE_SHORT_STRIPS on the raw record too: the decode zeroes the
    // short end, so the raw ADC trace drops the same side for comparability.
    if (Constants::cfg.IGNORE_SHORT_STRIPS)
      for (Int_t s = 1; s <= 16; s++) {
        if ((s % 2) != 0)
          e.rightdE_adc[s - 1] = 0.0f;
        else
          e.leftdE_adc[s - 1] = 0.0f;
      }
    // Both-channel multiplicity: segmented strips where both ends FIRED, read
    // off RAW ADC, not calibrated: short-end gains are 0 (no sim anchor).
    Int_t both = 0;
    for (Int_t k = 0; k < 16; k++)
      if (ev.leftdE_adc[k] > 0 && ev.rightdE_adc[k] > 0)
        both++;
    e.both_mult = both;
    e.seed_ts = seed_ts_in;
    e.reac_mask = mask;
    e.beam_flat = beam;
    res.reservoir.push_back(e);
    if (mask != 0)
      totalGated++;
  }
  res.gated = totalGated;
  res.seen = totalSeen;
  res.normed = totalNormed;
  return res;
}

// Run `n` indexed tasks on `workers` threads, pulling from a shared queue.
void StripSumScatter::RunIndexedParallel(
    Int_t n, Int_t workers, const std::function<void(Int_t)> &task) {
  std::queue<Int_t> work;
  for (Int_t i = 0; i < n; i++)
    work.push(i);
  std::mutex work_mutex;
  std::vector<std::thread> pool;
  for (Int_t w = 0; w < workers; w++) {
    pool.emplace_back([&]() {
      while (true) {
        Int_t i;
        {
          std::lock_guard<std::mutex> lk(work_mutex);
          if (work.empty())
            return;
          i = work.front();
          work.pop();
        }
        task(i);
      }
    });
  }
  for (Int_t w = 0; w < Int_t(pool.size()); w++)
    pool[w].join();
}

void StripSumScatter::AllocateScatters() {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kXLo = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_HI;
  const Int_t kXBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.XBINS;
  const Int_t kYBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.YBINS;
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    TH2F *h =
        new TH2F(Form("scatter_r%d", reac),
                 Form(";norm. #DeltaE strips %d#rightarrow%d [a.u.];norm. "
                      "#DeltaE strips %d#rightarrow%d [a.u.]",
                      kXLo, kXHi, YLoOf(reac), YHiOf(reac)),
                 kXBins, ScatterBuildRange::kXMin, ScatterBuildRange::kXMax,
                 kYBins, ScatterBuildRange::kYMin, ScatterBuildRange::kYMax);
    h->SetDirectory(nullptr);
    h->SetStats(0);
    m_scatter[reac] = h;
  }
}

void StripSumScatter::ReprojectFromReservoir() {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  Long64_t nFilled = 0;
  for (Int_t k = 0; k < Int_t(m_reservoir.size()); k++) {
    const TraceEvt &e = m_reservoir[k];
    if (e.reac_mask == 0)
      continue;
    Double_t total[18];
    e.Totals(total);
    for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
      if (!(e.reac_mask & (1u << ReacIndex(reac))))
        continue;
      Double_t x = 0.0, y = 0.0;
      PlaneXY(total, reac, x, y);
      m_scatter[reac]->Fill(x, y);
      nFilled++;
    }
  }
  std::cout << "strip-sum-scatter: re-projected " << nFilled
            << " tagged entries." << std::endl;
}

void StripSumScatter::FillScatters(const FileSet::GateGroups &groups) {
  const std::vector<Int_t> &runOrder = groups.order;
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kXBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.XBINS;
  const Int_t kYBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.YBINS;

  AllocateScatters();
  // The variant counters, in the order every per-run fill uses.
  {
    const std::vector<std::pair<TString, TagThresholds>> variants =
        ThresholdVariants();
    m_variantNames.clear();
    m_taggedVar.assign(variants.size(),
                       std::vector<Long64_t>(m_tagged.size(), 0));
    for (size_t v = 0; v < variants.size(); v++)
      m_variantNames.push_back(variants[v].first);
  }

  std::vector<GateSpec> activeGates = ActiveGates();

  // Pre-index the chains so worker threads only touch their own group: a
  // run's chunks on SOLARIS, one subfile on CoMPASS.
  const Int_t nRuns = Int_t(runOrder.size());
  std::vector<TChain *> chainVec(nRuns);
  std::vector<TString> labelVec(nRuns);
  for (Int_t i = 0; i < nRuns; i++) {
    chainVec[i] = groups.chain.find(runOrder[i])->second;
    labelVec[i] = groups.label.find(runOrder[i])->second;
  }

  Int_t n_workers =
      TMath::Min(Int_t(std::thread::hardware_concurrency()), nRuns);
  n_workers = TMath::Min(
      n_workers, Constants::cfg.STRIP_SUM_SCATTER_CONFIG.MAX_STRIP_SUM_WORKERS);
  if (n_workers < 1)
    n_workers = 1;
  std::cout << "strip-sum-scatter: " << nRuns << " gate groups on " << n_workers
            << " workers" << std::endl;

  // Phase 1: beam ellipses + series gates, one task per run. The gates within
  // a run stay sequential -- each only sees events passing the prior ones.
  std::vector<SingleRunFitResult> fits(nRuns);
  RunIndexedParallel(nRuns, n_workers, [&](Int_t i) {
    fits[i] = FitRunGates(runOrder[i], labelVec[i], chainVec[i], activeGates);
  });

  // Phase 2: fill, split per file, not per run: a CoMPASS run is hundreds of
  // subfiles in one chain (SOLARIS: one per run); tasks carry their gates.
  struct FillTask {
    Int_t run_idx;
    TString path;
  };
  std::vector<FillTask> tasks;
  {
    std::map<Int_t, Int_t> idx_of_run;
    for (Int_t i = 0; i < nRuns; i++)
      idx_of_run[runOrder[i]] = i;
    std::vector<FileSpec> specs = FileSet::BuildProcessedFileSpecs();
    for (Int_t k = 0; k < Int_t(specs.size()); k++) {
      // The file's gate group: its run on SOLARIS, itself on CoMPASS.
      std::map<TString, Int_t>::const_iterator kit =
          groups.key_of.find(FileSet::EventsName(specs[k]));
      if (kit == groups.key_of.end())
        continue;
      std::map<Int_t, Int_t>::const_iterator it = idx_of_run.find(kit->second);
      if (it == idx_of_run.end() || !fits[it->second].ok)
        continue;
      TString full = IO::GetRootFilesBaseDir() + "/" +
                     FileSet::EventsName(specs[k]) + ".root";
      if (gSystem->AccessPathName(full))
        continue;
      FillTask t;
      t.run_idx = it->second;
      t.path = full;
      tasks.push_back(t);
    }
  }
  Int_t nTasks = Int_t(tasks.size());
  Int_t fill_workers =
      TMath::Min(Int_t(std::thread::hardware_concurrency()), nTasks);
  fill_workers =
      TMath::Min(fill_workers,
                 Constants::cfg.STRIP_SUM_SCATTER_CONFIG.MAX_STRIP_SUM_WORKERS);
  if (fill_workers < 1)
    fill_workers = 1;
  std::cout << "strip-sum-scatter: filling " << nTasks << " files on "
            << fill_workers << " workers" << std::endl;

  // Merge into the totals in TASK ORDER as each prefix finishes, freeing on
  // the spot; the order keeps the pure-beam budget identical to sequential.
  std::vector<SingleRunFillResult> fills(nTasks);
  std::vector<Bool_t> filled(nTasks, kFALSE);
  Int_t next_merge = 0;
  std::mutex merge_mutex;
  Long64_t totalGated = 0, totalSeen = 0;
  const Int_t kBeamReservoirCap =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRACES_PER_CLASS * 10;
  Int_t nBeamKept = 0;
  auto merge = [&](Int_t t) {
    for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
      Int_t ri = ReacIndex(reac);
      if (ri < Int_t(fills[t].scatters.size()) && fills[t].scatters[ri])
        m_scatter[reac]->Add(fills[t].scatters[ri]);
    }
    // Pure-beam events are capped GLOBALLY, not per task (reaction-tagged are
    // never dropped); per-task capping would scale kept beam with file count.
    for (Int_t k = 0; k < Int_t(fills[t].reservoir.size()); k++) {
      const TraceEvt &e = fills[t].reservoir[k];
      if (e.beam_flat) {
        if (nBeamKept >= kBeamReservoirCap)
          continue;
        nBeamKept++;
      }
      m_reservoir.push_back(e);
    }
    totalGated += fills[t].gated;
    totalSeen += fills[t].seen;
    m_nNormed += fills[t].normed;
    for (Int_t k = 0; k < Int_t(fills[t].normed_at.size()); k++)
      m_normedAt[k] += fills[t].normed_at[k];
    for (Int_t k = 0; k < Int_t(fills[t].tagged.size()); k++)
      m_tagged[k] += fills[t].tagged[k];
    if (m_taggedVar.size() == fills[t].tagged_var.size())
      for (size_t v = 0; v < m_taggedVar.size(); v++)
        for (Int_t k = 0; k < Int_t(fills[t].tagged_var[v].size()); k++)
          m_taggedVar[v][k] += fills[t].tagged_var[v][k];
    for (Int_t k = 0; k < Int_t(fills[t].cut_counts.size()); k++)
      m_cutCounts[k] += fills[t].cut_counts[k];
    for (Int_t k = 0; k < Int_t(fills[t].pre_counts.size()); k++)
      m_preCounts[k] += fills[t].pre_counts[k];
    for (Int_t k = 0; k < Int_t(fills[t].scatters.size()); k++)
      delete fills[t].scatters[k];
    fills[t].scatters.clear();
    std::vector<TraceEvt>().swap(fills[t].reservoir);
  };
  RunIndexedParallel(nTasks, fill_workers, [&](Int_t t) {
    Int_t i = tasks[t].run_idx;
    {
      TChain ch("events");
      ch.Add(tasks[t].path);
      if (ch.GetEntries() > 0)
        fills[t] = FillRunScatters(runOrder[i], labelVec[i], &ch, activeGates,
                                   fits[i].series_gates, fits[i].pure_beam);
    }
    std::lock_guard<std::mutex> lk(merge_mutex);
    filled[t] = kTRUE;
    while (next_merge < nTasks && filled[next_merge]) {
      merge(next_merge);
      next_merge++;
    }
  });
  m_nSeen = totalSeen;
  std::cout << "strip-sum-scatter: " << totalGated << " reaction-tagged of "
            << totalSeen << " events (" << m_nNormed
            << " past every pre-tag cut); reservoir " << m_reservoir.size()
            << std::endl;
  WriteCutReport();
}

/// The per-condition counts of the last fill, printed and written to
/// plots/strip_sum_scatter/tag_cuts.txt. Sequential: each condition counts the
/// events that passed everything before it, so a row sums to the events that
/// reached the tag at that strip and the columns read as a funnel. A condition
/// that is off in the config shows "off" rather than 0.
void StripSumScatter::WriteCutReport() const {
  const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;
  const Int_t kReacMin = C.REACTION_STRIP_MIN;
  const Int_t kReacMax = C.REACTION_STRIP_MAX;
  const Bool_t active[kNTagCuts] = {
      kTRUE,
      kTRUE,
      C.REQUIRE_BEAM_UPSTREAM_OF_REAC,
      kTRUE,
      kTRUE,
      C.SMOOTHNESS_NSIGMA > 0.0,
      (C.TAIL_FALL_FROM_STRIP > 0) && C.TAIL_RISE_NSIGMA > 0.0,
      C.TAIL_RERISE_NSIGMA > 0.0,
      C.POST_ABOVE_NSIGMA > 0.0 && C.POST_ABOVE_STRIPS > 0,
      kTRUE,
      C.TAIL_CLIFF_MAX_FRACTION > 0.0};
  const Bool_t pre_active[kNPreCuts] = {kTRUE, kTRUE, kTRUE, kTRUE,
                                        C.BOTH_MULT_MAX >= 0};

  TString out;
  out += Form("strip-sum-scatter tag cuts: %s, %lld events seen\n",
              Paths::DatasetName().Data(), m_nSeen);
  out += "Sequential counts: each condition counts the events that passed "
         "everything before it.\n\n";
  out += "Event-level cuts (before any reaction is asked about):\n";
  for (Int_t c = 1; c < kNPreCuts; c++) {
    if (!pre_active[c])
      out += Form("  %-14s off\n", kPreCutName[c]);
    else
      out += Form("  %-14s %12lld  (%6.2f%% of seen)\n", kPreCutName[c],
                  m_preCounts[c],
                  m_nSeen > 0 ? 100.0 * m_preCounts[c] / m_nSeen : 0.0);
  }
  out += Form("  %-14s %12lld  (%6.2f%% of seen)\n\n", kPreCutName[kPrePass],
              m_preCounts[kPrePass],
              m_nSeen > 0 ? 100.0 * m_preCounts[kPrePass] / m_nSeen : 0.0);

  out += "Tag conditions per reaction strip (events past the event-level "
         "cuts; % of those):\n";
  out += "reac";
  for (Int_t c = 1; c < kNTagCuts; c++)
    out += Form(" | %13s", kTagCutName[c]);
  out += Form(" | %13s\n", kTagCutName[kTagPass]);
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    const Long64_t *row = &m_cutCounts[ReacIndex(reac) * kNTagCuts];
    Long64_t total = 0;
    for (Int_t c = 0; c < kNTagCuts; c++)
      total += row[c];
    out += Form("%4d", reac);
    for (Int_t c = 1; c < kNTagCuts; c++) {
      if (!active[c])
        out += Form(" | %13s", "off");
      else
        out += Form(" | %13s", Form("%lld (%.2f%%)", row[c],
                                    total > 0 ? 100.0 * row[c] / total : 0.0));
    }
    out += Form(" | %13s\n",
                Form("%lld (%.3f%%)", row[kTagPass],
                     total > 0 ? 100.0 * row[kTagPass] / total : 0.0));
  }
  std::cout << out;

  // The cut variation: tagged count per strip under each shifted threshold,
  // beside the nominal, so the sensitivity is readable before the cross
  // section folds it into a systematic.
  if (!m_variantNames.empty()) {
    out += "\nCut variation (tagged events per strip with one threshold "
           "shifted; nominal first):\n";
    out += "reac |       nominal";
    for (size_t v = 0; v < m_variantNames.size(); v++)
      out += Form(" | %9s", m_variantNames[v].Data());
    out += "\n";
    for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
      out += Form("%4d | %13lld", reac, m_tagged[ReacIndex(reac)]);
      for (size_t v = 0; v < m_variantNames.size(); v++)
        out += Form(" | %9lld", m_taggedVar[v][ReacIndex(reac)]);
      out += "\n";
    }
  }

  const TString dir = Paths::ResultsDir() + "/plots/strip_sum_scatter";
  gSystem->mkdir(dir, kTRUE);
  const TString path = dir + "/tag_cuts.txt";
  std::ofstream f(path.Data());
  if (!f) {
    std::cerr << "strip-sum-scatter: cannot write " << path << std::endl;
    return;
  }
  f << out;
  f.close();
  std::cout << "strip-sum-scatter: tag-cut report written to " << path
            << std::endl;
}

void StripSumScatter::PlotScatters() {
  // Display windows recomputed here, not in FillScatters: a cached run skips
  // the fill; constructor zeros made SetRangeUser(0, 0) throw the zoom away.
  YBounds(m_yLo, m_yHi);
  const Int_t kXLo = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_HI;
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;

  std::lock_guard<std::mutex> lock(g_plot_mutex);
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    // Display-only zoom: the histogram is built over the fixed build range, so
    // retuning these windows never forces a refill.
    Int_t ri = ReacIndex(reac);
    m_scatter[reac]->GetXaxis()->SetRangeUser(
        Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_DISPLAY_MIN,
        Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_DISPLAY_MAX);
    m_scatter[reac]->GetYaxis()->SetRangeUser(m_yLo[ri], m_yHi[ri]);
    PlottingUtils::ConfigureAndDraw2DHistogram(m_scatter[reac], c);
    m_scatter[reac]->GetYaxis()->SetTitleOffset(1.3);
    c->SetLeftMargin(0.18);
    PlottingUtils::SaveFigure(c,
                              Form("normsumE_reac%d_s%d_%d_vs_s%d_%d", reac,
                                   YLoOf(reac), YHiOf(reac), kXLo, kXHi),
                              "strip_sum_scatter", PlotSaveOptions::kLINEAR);
    delete c;
  }
}

void StripSumScatter::InteractiveOverlay(Int_t reac) {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kXLo = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_HI;
  const Int_t kTracesPerRegion =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRACES_PER_CLASS;

  if (reac < kReacMin || reac > kReacMax) {
    std::cerr << "strip-sum-scatter: candidate reaction strip " << reac
              << " outside [" << kReacMin << "," << kReacMax
              << "]; skipping interactive overlay." << std::endl;
    return;
  }
  // A cut saved from a previous pass wins over prompting, so a decided run
  // repeats without a person in the loop or DISPLAY: nothing interactive opens.
  TCutG *cutAn = nullptr;
  TCutG *cutAa = nullptr;
  // Saved cuts -- drawn here earlier, or fitted by compute-regions -- win over
  // prompting, so a decided run repeats without a person in the loop.
  if (!Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REGION_CUT_REDRAW) {
    cutAn = LoadRegionCut("region_an", reac);
    cutAa = LoadRegionCut("region_aa", reac);
  }
  TCanvas *cutCanvas = nullptr;

  // TApplication + argv must outlive all GUI canvases below, incl. trace ones:
  // tearing it down early faults ROOT's paint path; so leak it, argv static.
  static Int_t app_argc = 1;
  static char app_arg0[] = "strip-sum-scatter";
  static char *app_argv[] = {app_arg0};

  if (cutAn && cutAa) {
    std::cout << "  [region] loaded saved cuts for reac " << reac << std::endl;
  } else {
    delete cutAn;
    cutAn = nullptr;
    delete cutAa;
    cutAa = nullptr;
    if (!gSystem->Getenv("DISPLAY")) {
      std::cerr << "strip-sum-scatter: no saved region cuts for reac " << reac
                << " and no DISPLAY to draw them; skipping interactive "
                   "region-trace overlay (scatters already saved)."
                << std::endl;
      return;
    }
    // Intentionally not stored and never deleted -- see the note above.
    new TApplication("strip-sum-scatter", &app_argc, app_argv);
    gROOT->SetBatch(kFALSE);

    cutCanvas = new TCanvas("c_strip_sum_regions",
                            "Draw (a,n) then (a,a') regions", 900, 700);
    cutCanvas->SetLogz(kTRUE); // match the saved scatter's z-scale
    m_scatter[reac]->Draw("COLZ");
    cutCanvas->Update();
    cutAn = PromptCut(cutCanvas, "region_an", "(a,n)");
    cutAa = PromptCut(cutCanvas, "region_aa", "(a,a')");
    SaveRegionCuts(reac, cutAn, cutAa);

    cutCanvas->GetListOfPrimitives()->Remove(cutAn);
    cutCanvas->GetListOfPrimitives()->Remove(cutAa);
    gROOT->SetEditorMode();
    gSystem->ProcessEvents();
    // Re-entering batch with TApplication + GUI canvas faults in TPad::PaintBox
    // (TCanvas::Build border paint); deleting canvas or nulling gPad: no help.
    cutCanvas->Clear();
    gSystem->ProcessEvents();
  }

  // ClusterVarHists(reac, cutAa, cutAn, "strip_sum_scatter");

  std::vector<TGraph *> tr_an, tr_aa, tr_beam;
  // Same selected events, raw (un-normalized) ADC -- one entry per normed
  // trace, kept in lock-step so the two overlays show the identical events.
  std::vector<TGraph *> tr_an_adc, tr_aa_adc, tr_beam_adc;
  // Same selected events, Savitzky-Golay smoothed (normed a.u. space), for
  // the smoothing overlay; in lock-step with raw normed traces, same events.
  const Bool_t kSkipSg =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.SKIP_SAVGOL_PLOTS;
  std::vector<TGraph *> tr_an_sg, tr_aa_sg, tr_beam_sg;
  UInt_t bit = (1u << ReacIndex(reac));

  for (Int_t k = 0; k < Int_t(m_reservoir.size()); k++) {
    if (Int_t(tr_an.size()) >= kTracesPerRegion &&
        Int_t(tr_aa.size()) >= kTracesPerRegion &&
        Int_t(tr_beam.size()) >= kTracesPerRegion)
      break;
    const TraceEvt &e = m_reservoir[k];
    Double_t td[18], td_adc[18];
    e.Totals(td);
    e.TotalsAdc(td_adc);
    if (e.beam_flat && Int_t(tr_beam.size()) < kTracesPerRegion) {
      tr_beam.push_back(TraceFromTotal(td));
      tr_beam_adc.push_back(TraceFromTotal(td_adc));
      if (!kSkipSg)
        tr_beam_sg.push_back(SmoothedTraceFromTotal(td));
      continue;
    }
    if (!(e.reac_mask & bit))
      continue;
    Double_t x = 0.0, y = 0.0;
    PlaneXY(td, reac, x, y);
    if (cutAn && Int_t(tr_an.size()) < kTracesPerRegion &&
        cutAn->IsInside(x, y)) {
      tr_an.push_back(TraceFromTotal(td));
      tr_an_adc.push_back(TraceFromTotal(td_adc));
      if (!kSkipSg)
        tr_an_sg.push_back(SmoothedTraceFromTotal(td));
    } else if (cutAa && Int_t(tr_aa.size()) < kTracesPerRegion &&
               cutAa->IsInside(x, y)) {
      tr_aa.push_back(TraceFromTotal(td));
      tr_aa_adc.push_back(TraceFromTotal(td_adc));
      if (!kSkipSg)
        tr_aa_sg.push_back(SmoothedTraceFromTotal(td));
    }
  }

  std::cout << "Sampled traces: beam=" << tr_beam.size()
            << " (a,a')=" << tr_aa.size() << " (a,n)=" << tr_an.size()
            << std::endl;
  DrawRegionTraces(Form("region_traces_reac%d", reac), "strip_sum_scatter",
                   tr_beam, tr_aa, tr_an, 0.6, 1.6, "#DeltaE [a.u.]");
  const Bool_t kMeanTraces =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PLOT_REGION_MEAN_TRACES;
  if (kMeanTraces)
    DrawRegionMeanTraces(Form("region_mean_traces_reac%d", reac),
                         "strip_sum_scatter", tr_beam, tr_aa, tr_an, 0.6, 1.6,
                         "#DeltaE [a.u.]");
  if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PLOT_ADC_TRACES) {
    Double_t adc_y_lo = 0.0, adc_y_hi = 0.0;
    TraceYRange(tr_beam_adc, tr_aa_adc, tr_an_adc, adc_y_lo, adc_y_hi);
    DrawRegionTraces(Form("region_traces_reac%d_adc", reac),
                     "strip_sum_scatter", tr_beam_adc, tr_aa_adc, tr_an_adc,
                     adc_y_lo, adc_y_hi, "#DeltaE [ADC]");
    if (kMeanTraces)
      DrawRegionMeanTraces(Form("region_mean_traces_reac%d_adc", reac),
                           "strip_sum_scatter", tr_beam_adc, tr_aa_adc,
                           tr_an_adc, adc_y_lo, adc_y_hi, "#DeltaE [ADC]");
  }
  if (!kSkipSg) {
    DrawRegionTraces(Form("region_traces_reac%d_sg", reac), "strip_sum_scatter",
                     tr_beam_sg, tr_aa_sg, tr_an_sg, 0.6, 1.6, "#DeltaE [a.u.]",
                     kFALSE);
    if (kMeanTraces)
      DrawRegionMeanTraces(Form("region_mean_traces_reac%d_sg", reac),
                           "strip_sum_scatter", tr_beam_sg, tr_aa_sg, tr_an_sg,
                           0.7, 1.3, "#DeltaE [a.u.]");
  }

  for (Int_t i = 0; i < Int_t(tr_an.size()); i++)
    delete tr_an[i];
  for (Int_t i = 0; i < Int_t(tr_aa.size()); i++)
    delete tr_aa[i];
  for (Int_t i = 0; i < Int_t(tr_beam.size()); i++)
    delete tr_beam[i];
  for (Int_t i = 0; i < Int_t(tr_an_adc.size()); i++)
    delete tr_an_adc[i];
  for (Int_t i = 0; i < Int_t(tr_aa_adc.size()); i++)
    delete tr_aa_adc[i];
  for (Int_t i = 0; i < Int_t(tr_beam_adc.size()); i++)
    delete tr_beam_adc[i];
  for (Int_t i = 0; i < Int_t(tr_an_sg.size()); i++)
    delete tr_an_sg[i];
  for (Int_t i = 0; i < Int_t(tr_aa_sg.size()); i++)
    delete tr_aa_sg[i];
  for (Int_t i = 0; i < Int_t(tr_beam_sg.size()); i++)
    delete tr_beam_sg[i];

  delete cutAn;
  delete cutAa;
}

Bool_t StripSumScatter::Prepare() {
  // Required before the threaded fill touches TChains from worker threads.
  ROOT::EnableThreadSafety();
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                Paths::ResultsDir() + "/plots",
                                Paths::ResultsDir() + "/root_files");
  gROOT->SetBatch(kTRUE);

  // The gate groups: a run's chunks on SOLARIS, one subfile on CoMPASS. Every
  // per-group step below runs one task per group.
  FileSet::GateGroups groups = FileSet::GroupEventsForGating();
  const std::vector<Int_t> &run_order = groups.order;
  std::map<Int_t, TChain *> &chain_by_run = groups.chain;
  if (run_order.empty()) {
    std::cerr << "strip-sum-scatter: no runs found" << std::endl;
    return kFALSE;
  }
  const TString first_label = groups.label[run_order[0]];

  // Every level and step cut is in sigma of the beam, so the beam comes
  // first: the first group's entrance and exit ellipses select its beam
  // events, and each strip's mean and width over them are the reference.
  // The fingerprint stamps the resolved thresholds.
  Double_t strip_mean[18], strip_sigma[18];
  {
    std::vector<GateSpec> no_gates;
    SingleRunFitResult first = FitRunGates(
        run_order[0], first_label, chain_by_run[run_order[0]], no_gates);
    if (!first.pure_beam.ok ||
        !MeasureBeamNoise(chain_by_run[run_order[0]], first.pure_beam,
                          strip_mean, strip_sigma)) {
      std::cerr << "strip-sum-scatter: cannot measure the beam on "
                << first_label << std::endl;
      for (Int_t i = 0; i < Int_t(run_order.size()); i++)
        delete chain_by_run[run_order[i]];
      return kFALSE;
    }
  }
  SetStripNoise(strip_mean, strip_sigma);
  {
    const Int_t kReacMin =
        Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
    const Int_t kReacMax =
        Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
    TString line = Form("strip-sum-scatter: beam mean (sigma) from the gated "
                        "sample of %s:",
                        first_label.Data());
    for (Int_t strip = 1; strip <= 17; strip++)
      line += Form(" s%d %.4f (%.4f)", strip, strip_mean[strip],
                   strip_sigma[strip]);
    std::cout << line << std::endl;
    std::cout << Form("strip-sum-scatter: jump gate %.2f sigma -> %.4f at "
                      "strip %d, %.4f at strip %d",
                      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REAC_JUMP_NSIGMA,
                      JumpMin(kReacMin), kReacMin, JumpMin(kReacMax), kReacMax)
              << std::endl;
  }

  // Build fingerprint and try cache.
  TString fingerprint = BuildFingerprint(groups);
  TString cache_name = CacheName();

  Bool_t loaded = TryLoadCache(cache_name, fingerprint);

  if (!loaded) {
    FillScatters(groups);
    WriteCache(cache_name, fingerprint);
  }

  // Batch plotting (always done).
  PlotScatters();

  // The chains served the beam measurement and the fill; the overlays work
  // from the scatters and the reservoir.
  for (Int_t i = 0; i < Int_t(run_order.size()); i++)
    delete chain_by_run[run_order[i]];
  return kTRUE;
}

void StripSumScatter::Run() {
  if (!Prepare())
    return;

  // Optional sim overlays.
  if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.RERUN_SIM) {
    SimOverlay();
    SimTraceOverlay();
  }

  // Interactive region-trace overlay (requires DISPLAY).
  Int_t reac = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.CANDIDATE_REAC_STRIP;
  InteractiveOverlay(reac);
}

void StripSumScatter::DrawAllTaggedTraces(const TString &subdir) {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kTracesPerRegion =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRACES_PER_CLASS;
  const Bool_t kMeanTraces =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PLOT_REGION_MEAN_TRACES;
  const std::vector<TGraph *> none;

  // The beam sample once: it is the same band under every strip's figure.
  std::vector<TGraph *> tr_beam, tr_beam_adc;
  for (Int_t k = 0; k < Int_t(m_reservoir.size()) &&
                    Int_t(tr_beam.size()) < kTracesPerRegion;
       k++) {
    const TraceEvt &e = m_reservoir[k];
    if (!e.beam_flat)
      continue;
    Double_t td[18], td_adc[18];
    e.Totals(td);
    e.TotalsAdc(td_adc);
    tr_beam.push_back(TraceFromTotal(td));
    tr_beam_adc.push_back(TraceFromTotal(td_adc));
  }

  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    const UInt_t bit = (1u << ReacIndex(reac));
    std::vector<TGraph *> tr_tag, tr_tag_adc;
    // Every tagged event at this strip, uncapped: in the all-tagged mode the
    // count IS these traces, so the figure shows all of them. Only the beam
    // sample above is capped.
    for (Int_t k = 0; k < Int_t(m_reservoir.size()); k++) {
      const TraceEvt &e = m_reservoir[k];
      if (e.beam_flat || !(e.reac_mask & bit))
        continue;
      Double_t td[18], td_adc[18];
      e.Totals(td);
      e.TotalsAdc(td_adc);
      tr_tag.push_back(TraceFromTotal(td));
      tr_tag_adc.push_back(TraceFromTotal(td_adc));
    }
    std::cout << "  reac " << reac << ": " << tr_tag.size()
              << " tagged traces drawn against " << tr_beam.size()
              << " beam traces" << std::endl;
    DrawRegionTraces(Form("all_tagged_traces_reac%d", reac), subdir, tr_beam,
                     none, tr_tag, 0.6, 1.6, "#DeltaE [a.u.]", kTRUE, "Tagged");
    if (kMeanTraces)
      DrawRegionMeanTraces(Form("all_tagged_mean_traces_reac%d", reac), subdir,
                           tr_beam, none, tr_tag, 0.6, 1.6, "#DeltaE [a.u.]");
    if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PLOT_ADC_TRACES) {
      Double_t adc_y_lo = 0.0, adc_y_hi = 0.0;
      TraceYRange(tr_beam_adc, none, tr_tag_adc, adc_y_lo, adc_y_hi);
      DrawRegionTraces(Form("all_tagged_traces_reac%d_adc", reac), subdir,
                       tr_beam_adc, none, tr_tag_adc, adc_y_lo, adc_y_hi,
                       "#DeltaE [ADC]", kTRUE, "Tagged");
      if (kMeanTraces)
        DrawRegionMeanTraces(Form("all_tagged_mean_traces_reac%d_adc", reac),
                             subdir, tr_beam_adc, none, tr_tag_adc, adc_y_lo,
                             adc_y_hi, "#DeltaE [ADC]");
    }
    for (Int_t i = 0; i < Int_t(tr_tag.size()); i++)
      delete tr_tag[i];
    for (Int_t i = 0; i < Int_t(tr_tag_adc.size()); i++)
      delete tr_tag_adc[i];
  }
  for (Int_t i = 0; i < Int_t(tr_beam.size()); i++)
    delete tr_beam[i];
  for (Int_t i = 0; i < Int_t(tr_beam_adc.size()); i++)
    delete tr_beam_adc[i];
}
