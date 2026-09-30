#include "StripSumScatter.hpp"
#include "RegionCuts.hpp"
#include "SelectionDiagram.hpp"
#include <TDecompLU.h>
#include <TMatrixD.h>
#include <TParameter.h>
#include <TRandom3.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>

const char *const kTagCutName[kNTagCuts] = {
    "tagged",    "upstream beam", "jump",       "reac level",
    "tail rise", "re-rise",       "post above", "beam crossing",
    "end strip", "cliff",         "other strip"};
const char *const kPreCutName[kNPreCuts] = {
    "reached tag", "all strips", "beam gate", "pileup",
    "noise",       "smoothness", "both mult"};

namespace {
// A beam-gate axis by name: "s<N>" / "strip N" for a strip, "grid" for the
// grid (GATE_AXIS_GRID). The figure and histogram names keep the strip form.
TString AxisTag(Int_t axis) {
  return axis == GATE_AXIS_GRID ? TString("grid") : TString(Form("s%d", axis));
}
TString AxisName(Int_t axis) {
  return axis == GATE_AXIS_GRID ? TString("grid")
                                : TString(Form("strip %d", axis));
}
// Axis title, as the events-summary figures write it.
TString AxisTitle(Int_t axis) {
  return axis == GATE_AXIS_GRID ? TString("Grid")
                                : TString(Form("Strip%d", axis));
}
// The pure-beam entrance ellipse's axes, per PURE_BEAM_GATE.
void EntranceAxes(Int_t &sx, Int_t &sy) {
  switch (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PURE_BEAM_GATE) {
  case StripSumScatterConfig::PURE_BEAM_GATE_S1_S2:
    sx = 1;
    sy = 2;
    break;
  case StripSumScatterConfig::PURE_BEAM_GATE_S1_GRID:
    sx = 1;
    sy = GATE_AXIS_GRID;
    break;
  default:
    sx = 0;
    sy = 1;
  }
}
} // namespace

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
  T.gate_nsigma = C.GATE_NSIGMA;
  T.pileup_nsigma = C.PILEUP_NSIGMA;
  T.noise_nsigma = C.NOISE_NSIGMA;
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
  T.cross_min_strip = C.POST_CROSS_MIN_STRIP;
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
      // The beam selection: the event-level cuts every count and denominator
      // pass through.
      {"gate", &TagThresholds::gate_nsigma, nom.gate_nsigma > 0.0, ds, kFALSE},
      {"pileup", &TagThresholds::pileup_nsigma, nom.pileup_nsigma > 0.0, ds,
       kFALSE},
      {"noise", &TagThresholds::noise_nsigma, nom.noise_nsigma > 0.0, ds,
       kFALSE},
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
      // A count: the variation is one strip either way (never down to off).
      {"cross", &TagThresholds::cross_min_strip, nom.cross_min_strip > 0.0, 1.0,
       kFALSE},
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

Bool_t StripSumScatter::ConditionActive(TagCut c, const TagThresholds &T) {
  switch (c) {
  case kCutUpstream:
    return T.require_upstream;
  case kCutJump:
  case kCutReacLevel:
  case kCutEndStrip:
    return kTRUE;
  case kCutTailRise:
    return T.tail_rise_nsigma > 0.0 && T.tail_fall_from_strip > 0;
  case kCutRerise:
    return T.tail_rerise_nsigma > 0.0;
  case kCutPostAbove:
    return T.post_above_nsigma > 0.0 && T.post_above_strips > 0;
  case kCutCross:
    return T.cross_min_strip > 0.0;
  case kCutCliff:
    return T.cliff_max > 0.0;
  default:
    return kFALSE;
  }
}

/// Each condition on its own. The tail-shape ones follow the published 87Rb
/// per-strip macros (TracesVisu2.C at A9-A11), each in sigma of the measured
/// noise; all are numerator-only: they describe the residue, not the beam. A
/// strip whose noise was not measured is skipped. A condition that does not
/// apply at this strip (a window past the last strip, the cliff with no
/// strip between the reaction and the end) passes.
Bool_t StripSumScatter::Condition(TagCut c, const EnergyView &ev, Int_t reac,
                                  const TagThresholds &T) {
  const Int_t kLast = Constants::cfg.IGNORE_STRIP_17 ? 16 : 17;
  const Int_t end_strip =
      (Constants::cfg.IGNORE_STRIP_17 ||
       Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REQUIRE_STRIP_16_BELOW_BEAM)
          ? 16
          : 17;
  switch (c) {
  case kCutUpstream:
    // Otherwise a reaction at an earlier strip can pass this strip's jump
    // gate on a noise fluctuation and be counted here as well.
    return BeamUpstreamOf(ev, reac, T);
  case kCutJump:
    return ev.Total(reac) - ev.Total(reac - 1) >
           T.jump_nsigma * StripSigma(reac);
  case kCutReacLevel:
    return ev.Total(reac) > StripMean(reac) + T.jump_nsigma * StripSigma(reac);
  case kCutTailRise: {
    // Falling tail: no rise beyond the noise from the start strip on.
    const Int_t start = TMath::Max(reac + 1, T.tail_fall_from_strip + 1);
    for (Int_t s = start; s <= kLast; s++)
      if (StripSigma(s) > 0.0 &&
          ev.Total(s) - ev.Total(s - 1) > T.tail_rise_nsigma * StripSigma(s))
        return kFALSE;
    return kTRUE;
  }
  case kCutRerise: {
    // No return: back at the beam after the peak, then above it again.
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
        return kFALSE;
      if (ev.Total(s) <= StripMean(s) + T.tail_return_nsigma * StripSigma(s))
        returned = kTRUE;
    }
    return kTRUE;
  }
  case kCutPostAbove:
    // Persistence: the excess must hold for post_above_strips strips
    // after the reaction, each more than post_above_nsigma of its spread.
    for (Int_t s = reac + 1; s <= TMath::Min(reac + T.post_above_strips, kLast);
         s++)
      if (StripSigma(s) > 0.0 &&
          !(ev.Total(s) > StripMean(s) + T.post_above_nsigma * StripSigma(s)))
        return kFALSE;
    return kTRUE;
  case kCutCross: {
    // The residue's range: the trace may not read below the beam before
    // strip cross_min_strip (an absolute strip); the window is reac+1..S.
    const Int_t upto = TMath::Min(Int_t(T.cross_min_strip + 0.5), kLast);
    for (Int_t s = reac + 1; s <= upto; s++)
      if (StripSigma(s) > 0.0 && ev.Total(s) < StripMean(s))
        return kFALSE;
    return kTRUE;
  }
  case kCutEndStrip:
    // The end strip below the beam by n sigma of its own spread: a residue
    // that has stopped (or is stopping) reads low there, the beam does not.
    return ev.Total(end_strip) <
           StripMean(end_strip) - T.end_strip_nsigma * StripSigma(end_strip);
  case kCutCliff: {
    // A stop is gradual: the last step may take only part of the fall from
    // the peak; the elastic class keeps its excess to end_strip-1.
    if (!(end_strip - 1 > reac))
      return kTRUE;
    Double_t peak = ev.Total(reac);
    for (Int_t s = reac + 1; s < end_strip; s++)
      peak = TMath::Max(peak, ev.Total(s));
    const Double_t fall = peak - ev.Total(end_strip);
    const Double_t last = ev.Total(end_strip - 1) - ev.Total(end_strip);
    return !(fall > 0.0 && last / fall > T.cliff_max);
  }
  default:
    return kTRUE;
  }
}

TagCut StripSumScatter::TailReason(const EnergyView &ev, Int_t reac,
                                   const TagThresholds &T) {
  for (Int_t c = kCutTailRise; c <= kCutPostAbove; c++)
    if (ConditionActive(TagCut(c), T) && !Condition(TagCut(c), ev, reac, T))
      return TagCut(c);
  return kTagPass;
}

// Sequential: the first active condition that fails, in enum order, so the
// cut report reads as a funnel.
TagCut StripSumScatter::RejectReason(const EnergyView &ev, Int_t reac,
                                     const TagThresholds &T) {
  for (Int_t c = kCutUpstream; c <= kCutCliff; c++)
    if (ConditionActive(TagCut(c), T) && !Condition(TagCut(c), ev, reac, T))
      return TagCut(c);
  return kTagPass;
}

TagCut StripSumScatter::RejectReason(const EnergyView &ev, Int_t reac) {
  return RejectReason(ev, reac, NominalThresholds());
}

Int_t StripSumScatter::ResolveTag(const EnergyView &ev,
                                  const std::vector<Bool_t> &pass) {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  Int_t best = -1;
  Double_t best_jump = -1.0e9;
  for (Int_t k = 0; k < Int_t(pass.size()); k++) {
    if (!pass[k])
      continue;
    const Int_t reac = kReacMin + k;
    if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TAG_RESOLVE ==
        StripSumScatterConfig::TAG_RESOLVE_FIRST_STRIP)
      return reac;
    const Double_t sig = StripSigma(reac);
    const Double_t jump =
        (ev.Total(reac) - ev.Total(reac - 1)) / (sig > 0.0 ? sig : 1.0);
    if (jump > best_jump) {
      best_jump = jump;
      best = reac;
    }
  }
  return best;
}

namespace {

// One step of the selection description: the condition, what an event must
// satisfy to go on, and the cut count it counts to.
void AddSelectionStep(std::vector<SelectionStep> &out,
                      SelectionStep::Stage stage, const char *id,
                      const char *name, const TString &detail, Bool_t on,
                      Int_t pre_cut, Int_t tag_cut) {
  SelectionStep st;
  st.stage = stage;
  st.id = id;
  st.name = name;
  st.detail = detail;
  st.on = on;
  st.pre_cut = pre_cut;
  st.tag_cut = tag_cut;
  out.push_back(st);
}

// "fewer than 1" reads better as "none".
TString FewerStrips(Int_t n) {
  return n <= 1 ? TString("none") : TString(Form("fewer than %d", n));
}

} // namespace

// The selection in words, one step per condition in the order RejectReason
// and the fill apply them; disabled steps stay in the list with `on` false.
std::vector<SelectionStep> StripSumScatter::DescribeSelection() {
  const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;
  const Bool_t ign0 = Constants::cfg.IGNORE_STRIP_0;
  const Bool_t req0 = !ign0 && Constants::cfg.REQUIRE_STRIP_0;
  const Bool_t ign17 = Constants::cfg.IGNORE_STRIP_17;
  const Int_t last = ign17 ? 16 : 17;
  const Int_t end_strip = (ign17 || C.REQUIRE_STRIP_16_BELOW_BEAM) ? 16 : 17;
  std::vector<SelectionStep> out;

  Int_t ent_x = 0, ent_y = 1;
  EntranceAxes(ent_x, ent_y);
  const TString entrance = AxisName(ent_x) + " vs " + AxisName(ent_y);
  // Input: what the numbers below are measured against.
  AddSelectionStep(
      out, SelectionStep::kInput, "events", "calibrated events",
      Form("E(s) = strip total in beam units (beam = 1)%s; beam gates fitted "
           "per %s",
           Constants::ActiveIgnoreShortStrips() ? ", long ends only" : "",
           Constants::ActiveUseSolarisData() ? "run" : "subfile"),
      kTRUE, -1, -1);
  AddSelectionStep(
      out, SelectionStep::kInput, "beam_ref", "beam reference",
      Form("beam mean and σ per strip from the first %s's pure-beam sample: "
           "every strip fired, inside the entrance (%s) and exit (%s) "
           "ellipses at %.1f σ, 3σ-clipped",
           Constants::ActiveUseSolarisData() ? "run" : "subfile",
           entrance.Data(), ign17 ? "strips 15 vs 16" : "strips 16 vs 17",
           C.PURE_BEAM_NSIGMA),
      kTRUE, -1, -1);

  // Event level, in the order of the fill loop.
  AddSelectionStep(
      out, SelectionStep::kEventLevel, "all_strips", kPreCutName[kPreAllStrips],
      Form("every strip %d-%d read above zero", req0 ? 0 : 1, last), kTRUE,
      kPreAllStrips, -1);
  const TString gates =
      C.GATE_NSIGMA > 0.0
          ? TString(Form("strip %d within %.1f σ of its fitted beam peak",
                         C.GATE_STRIP, C.GATE_NSIGMA))
          : TString("no beam gate");
  AddSelectionStep(out, SelectionStep::kEventLevel, "gate",
                   kPreCutName[kPreGate], gates, kTRUE, kPreGate, -1);
  AddSelectionStep(
      out, SelectionStep::kEventLevel, "pileup", kPreCutName[kPrePileup],
      Form("%s of strips 1-16 at or above beam + %.1f σ",
           FewerStrips(C.PILEUP_MIN_STRIPS).Data(), C.PILEUP_NSIGMA),
      C.PILEUP_NSIGMA > 0.0, kPrePileup, -1);
  AddSelectionStep(out, SelectionStep::kEventLevel, "noise",
                   kPreCutName[kPreNoise],
                   Form("%s of strips 1-16 at or below beam − %.1f σ",
                        FewerStrips(C.NOISE_MIN_STRIPS).Data(), C.NOISE_NSIGMA),
                   C.NOISE_NSIGMA > 0.0, kPreNoise, -1);
  AddSelectionStep(
      out, SelectionStep::kEventLevel, "smooth", kPreCutName[kPreSmooth],
      Form("every step between strips 1-%d within %.1f σ of the later "
           "strip's spread, the single largest rise excepted",
           last, C.SMOOTHNESS_NSIGMA),
      C.SMOOTHNESS_NSIGMA > 0.0, kPreSmooth, -1);
  AddSelectionStep(
      out, SelectionStep::kEventLevel, "both_mult", kPreCutName[kPreBothMult],
      Form("at most %d of strips 1-%d with both ends fired (raw ADC)",
           C.BOTH_MULT_MAX, TMath::Min(16, C.BOTH_MULT_COUNT_TO)),
      C.BOTH_MULT_MAX >= 0, kPreBothMult, -1);

  // Per reaction strip, in the order of RejectReason. The loop marker first:
  // it carries the strip range and is not a cut.
  AddSelectionStep(
      out, SelectionStep::kPerStrip, "reac_loop",
      Form("for each reaction strip reac = %d..%d", C.REACTION_STRIP_MIN,
           C.REACTION_STRIP_MAX),
      "every strip is tested on its own; the first failing condition ends "
      "its test",
      kTRUE, -1, -1);
  AddSelectionStep(out, SelectionStep::kPerStrip, "upstream",
                   kTagCutName[kCutUpstream],
                   Form("every strip 1..reac−1 within %.1f σ of the beam",
                        C.BEAM_UPSTREAM_NSIGMA),
                   C.REQUIRE_BEAM_UPSTREAM_OF_REAC, -1, kCutUpstream);
  AddSelectionStep(
      out, SelectionStep::kPerStrip, "jump", kTagCutName[kCutJump],
      Form("E(reac) − E(reac−1) > %.1f σ(reac)", C.REAC_JUMP_NSIGMA), kTRUE, -1,
      kCutJump);
  AddSelectionStep(out, SelectionStep::kPerStrip, "reac_level",
                   kTagCutName[kCutReacLevel],
                   Form("E(reac) > beam + %.1f σ(reac)", C.REAC_JUMP_NSIGMA),
                   kTRUE, -1, kCutReacLevel);
  AddSelectionStep(
      out, SelectionStep::kPerStrip, "tail_rise", kTagCutName[kCutTailRise],
      Form("no step up above %.1f σ from strip max(reac, %d)+1 to %d",
           C.TAIL_RISE_NSIGMA, C.TAIL_FALL_FROM_STRIP, last),
      C.TAIL_FALL_FROM_STRIP > 0 && C.TAIL_RISE_NSIGMA > 0.0, -1, kCutTailRise);
  AddSelectionStep(
      out, SelectionStep::kPerStrip, "rerise", kTagCutName[kCutRerise],
      Form("once back within %.1f σ of the beam after the peak, never above "
           "beam + %.1f σ again",
           C.TAIL_RETURN_NSIGMA, C.TAIL_RERISE_NSIGMA),
      C.TAIL_RERISE_NSIGMA > 0.0, -1, kCutRerise);
  AddSelectionStep(
      out, SelectionStep::kPerStrip, "post_above", kTagCutName[kCutPostAbove],
      Form("strips reac+1 to reac+%d all above beam + %.1f σ",
           C.POST_ABOVE_STRIPS, C.POST_ABOVE_NSIGMA),
      C.POST_ABOVE_NSIGMA > 0.0 && C.POST_ABOVE_STRIPS > 0, -1, kCutPostAbove);
  AddSelectionStep(
      out, SelectionStep::kPerStrip, "cross", kTagCutName[kCutCross],
      Form("strips reac+1 through %d all at or above the beam mean: the "
           "residue reaches strip %d before crossing the beam",
           C.POST_CROSS_MIN_STRIP, C.POST_CROSS_MIN_STRIP),
      C.POST_CROSS_MIN_STRIP > 0, -1, kCutCross);
  AddSelectionStep(
      out, SelectionStep::kPerStrip, "last_strip", kTagCutName[kCutEndStrip],
      Form("strip %d below beam − %.1f σ", end_strip, C.END_STRIP_NSIGMA),
      kTRUE, -1, kCutEndStrip);
  AddSelectionStep(
      out, SelectionStep::kPerStrip, "cliff", kTagCutName[kCutCliff],
      Form("the last step, strip %d to %d, at most %.0f%% of the fall from "
           "the peak",
           end_strip - 1, end_strip, 100.0 * C.TAIL_CLIFF_MAX_FRACTION),
      C.TAIL_CLIFF_MAX_FRACTION > 0.0, -1, kCutCliff);
  AddSelectionStep(
      out, SelectionStep::kPerStrip, "one_tag", kTagCutName[kCutOtherTag],
      C.TAG_RESOLVE == StripSumScatterConfig::TAG_RESOLVE_FIRST_STRIP
          ? "one tag per event: the first strip that passes takes it"
          : "one tag per event: the passing strip with the largest jump takes "
            "it",
      kTRUE, -1, kCutOtherTag);

  // What becomes of a tag.
  AddSelectionStep(out, SelectionStep::kOutcome, "tagged", "tagged",
                   "a reaction at reac; one strip per event", kTRUE, -1,
                   kTagPass);
  AddSelectionStep(
      out, SelectionStep::kOutcome, "n_beam", "beam count",
      "events past the event-level cuts with beam upstream of reac and no "
      "tag before it: the beam incident on reac, the denominator",
      kTRUE, -1, -1);
  if (C.AN_REGION_MODE == StripSumScatterConfig::AN_REGION_ALL_TAGGED) {
    TString sys = "no cut variation";
    if (C.CUT_VARIATION)
      sys = Form("systematic: each beam-selection and identification "
                 "condition shifted ±%.1f σ (cliff ±%.2f), one at a time, "
                 "largest change of σ per condition in quadrature, with the "
                 "gas pressure uncertainty",
                 C.CUT_VARIATION_NSIGMA_STEP, C.CUT_VARIATION_CLIFF_STEP);
    AddSelectionStep(
        out, SelectionStep::kOutcome, "xs", "cross section, all tagged",
        "σ(reac) = tagged / (beam count × target atoms per strip); " + sys,
        kTRUE, -1, -1);
  } else {
    AddSelectionStep(
        out, SelectionStep::kOutcome, "xs", "cross section, mixture fit",
        Form("scatter x = ΣE(%d..%d), y = ΣE(reac+1..reac+%d)%s; bivariate "
             "Gaussian mixture per strip, (a,n) = reaction component within "
             "%.1f σ, count as attributed by the fit",
             C.X_LO, C.X_HI, C.POST_TRIGGER_SUM_STRIPS,
             C.Y_RATIO_TO_UPSTREAM ? " over the upstream mean" : "",
             C.AN_REGION_NSIGMA),
        kTRUE, -1, -1);
  }
  return out;
}

// Sigma-clipped mean and width, started from the median and MAD, not the
// RMS: never sheds a second population (0.41 a.u. vs 0.06 on 37Cl 97).
static Double_t ClippedWidth(const std::vector<Double_t> &values,
                             Double_t *mean_out = nullptr) {
  const Int_t kClipPasses = 3;
  const Double_t kClipNSigma = 3.0;
  // The variance a Gaussian keeps inside the clip (0.973 at 3 sigma), so the
  // width is a Gaussian's sigma, as BeamFitUtils::ClippedMoments' is.
  const Double_t kKeptVar =
      1.0 - 2.0 * kClipNSigma * TMath::Gaus(kClipNSigma, 0.0, 1.0, kTRUE) /
                TMath::Erf(kClipNSigma / TMath::Sqrt2());
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
    width = TMath::Sqrt(TMath::Max(0.0, sum2 / Double_t(kept) - mean * mean) /
                        kKeptVar);
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
  chain->SetBranchStatus("Grid", 1);
  // Absent in files built before the branch existed; harmless to enable.
  if (chain->GetBranch("SeedTs"))
    chain->SetBranchStatus("SeedTs", 1);
}

Bool_t StripSumScatter::AllStripsFired(const EnergyView &ev) {
  if (!Constants::cfg.IGNORE_STRIP_0 && Constants::cfg.REQUIRE_STRIP_0 &&
      !(ev.strip0 > 0.0))
    return kFALSE;
  if (!Constants::cfg.IGNORE_STRIP_17 && !(ev.strip17 > 0.0))
    return kFALSE;
  for (Int_t s = 1; s <= 16; s++)
    if (!(ev.Total(s) > 0.0))
      return kFALSE;
  return kTRUE;
}

Bool_t StripSumScatter::IsPureBeam(const EnergyView &ev,
                                   const BeamEllipses &be) {
  if (!be.ok)
    return kFALSE;
  if (!AllStripsFired(ev))
    return kFALSE;
  // Must pass BOTH the entrance AND exit ellipses, at the beam sample's own
  // level (PURE_BEAM_NSIGMA), not the event-level gate's.
  const Double_t n = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PURE_BEAM_NSIGMA;
  switch (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PURE_BEAM_GATE) {
  case StripSumScatterConfig::PURE_BEAM_GATE_S1_S2:
    if (!PassesEllipse(be.s1_s2, ev, 1, 2, n))
      return kFALSE;
    break;
  case StripSumScatterConfig::PURE_BEAM_GATE_S1_GRID:
    if (!PassesEllipse(be.s1_grid, ev, 1, GATE_AXIS_GRID, n))
      return kFALSE;
    break;
  default:
    if (!PassesEllipse(be.s0_s1, ev, 0, 1, n))
      return kFALSE;
  }
  if (be.use_s15_s16) {
    if (!PassesEllipse(be.s15_s16, ev, 15, 16, n))
      return kFALSE;
  } else {
    if (!PassesEllipse(be.s16_s17, ev, 16, 17, n))
      return kFALSE;
  }
  return kTRUE;
}

Bool_t StripSumScatter::IsPileup(const EnergyView &ev) {
  return IsPileup(ev, Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PILEUP_NSIGMA);
}

Bool_t StripSumScatter::IsPileup(const EnergyView &ev, Double_t kNSigma) {
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
  return IsNoise(ev, Constants::cfg.STRIP_SUM_SCATTER_CONFIG.NOISE_NSIGMA);
}

Bool_t StripSumScatter::IsNoise(const EnergyView &ev, Double_t kNSigma) {
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
// The largest step with its single largest rise exempted: a reaction is
// one rise then gradual change; strip 0 stays out (its own scale).
Double_t StripSumScatter::MaxStepNSigma(const EnergyView &ev) {
  const Int_t kLast = Constants::cfg.IGNORE_STRIP_17 ? 16 : 17;
  Double_t worst = 0.0, best_rise = 0.0, second_rise = 0.0;
  for (Int_t s = 2; s <= kLast; s++) {
    if (!(StripSigma(s) > 0.0))
      continue;
    const Double_t d = (ev.Total(s) - ev.Total(s - 1)) / StripSigma(s);
    if (d > 0.0) {
      // Rises are ranked: the largest is exempt, the runner-up counts.
      if (d > best_rise) {
        second_rise = best_rise;
        best_rise = d;
      } else if (d > second_rise) {
        second_rise = d;
      }
    } else if (-d > worst) {
      worst = -d;
    }
  }
  return TMath::Max(worst, second_rise);
}

Bool_t StripSumScatter::IsSmooth(const EnergyView &ev, Double_t nsigma) {
  return !(nsigma > 0.0) || MaxStepNSigma(ev) <= nsigma;
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

TString StripSumScatter::CacheName() {
  // A tagged epoch keeps its own cache: its scatters and reservoir are of
  // different events files than the untagged eras'.
  TString name = "StripSumScatter_cache";
  if (Constants::ActiveFileTag().Length() > 0)
    name += "_" + Constants::ActiveFileTag();
  name += ".root";
  return name;
}

Bool_t StripSumScatter::PassesGate(const BeamGate1D &gate, const EnergyView &ev,
                                   Int_t strip, Double_t nsigma) {
  if (!(nsigma > 0.0))
    return kTRUE;
  const Double_t v = ev.Axis(strip);
  if (!(v > 0.0) || !gate.ok)
    return kFALSE;
  return TMath::Abs(v - gate.mu) < nsigma * gate.sigma;
}

Bool_t StripSumScatter::PassesEllipse(const BeamFit2D &gate,
                                      const EnergyView &ev, Int_t sx, Int_t sy,
                                      Double_t nsigma) {
  Double_t g0 = ev.Axis(sx);
  Double_t g1 = ev.Axis(sy);
  if (!(g0 > 0.0 && g1 > 0.0))
    return kFALSE;
  return BeamFitUtils::InEllipseXY(gate, g0, g1, nsigma);
}

BeamFit2D StripSumScatter::FindBeamGate(TChain *chain, Int_t sx, Int_t sy,
                                        const TString &tag,
                                        const TString &subdir,
                                        Double_t nsigma) {
  // The figure histogram is coarse; the fit is on a fine one around the seed.
  const Int_t kGateBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_BINS;
  const Double_t kGateMin = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_MIN;
  const Double_t kGateMax = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_MAX;
  const Int_t kSeedHalfBins =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.SEED_HALF_BINS;
  const Long64_t kSampleMaxPoints =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.SAMPLE_MAX_POINTS;

  BeamFit2D out;
  EnergyView ev;
  ev.Attach(chain);
  EnableEventBranches(chain);
  TH2F *h =
      new TH2F(Form("h2_beamgate_%s_%s_%s", AxisTag(sx).Data(),
                    AxisTag(sy).Data(), tag.Data()),
               Form(";%s #DeltaE [a.u.];%s #DeltaE [a.u.]",
                    AxisTitle(sx).Data(), AxisTitle(sy).Data()),
               kGateBins, kGateMin, kGateMax, kGateBins, kGateMin, kGateMax);
  h->SetDirectory(nullptr);
  Long64_t n = chain->GetEntries();
  Long64_t stride = FileSet::SampleStride(n, kSampleMaxPoints);
  // The same events unbinned, for the clipped moments below.
  std::vector<std::pair<Float_t, Float_t>> pts;
  pts.reserve(std::min(n, kSampleMaxPoints > 0 ? kSampleMaxPoints : n));
  for (Long64_t j = 0; j < n; j += stride) {
    chain->GetEntry(j);
    ev.Decode();
    Double_t x = ev.Axis(sx);
    Double_t y = ev.Axis(sy);
    if (x > 0.0 && y > 0.0) {
      h->Fill(x, y);
      pts.push_back(std::make_pair(Float_t(x), Float_t(y)));
    }
  }
  if (h->GetEntries() < 100) {
    delete h;
    return out;
  }
  Moments2D m = BeamFitUtils::SeedSpotMoments(h, kSeedHalfBins);
  if (m.weight <= 0) {
    delete h;
    return out;
  }
  // The width the cuts are in is a fitted sigma: the spot events, finely
  // binned around the seed, fitted with a correlated Gaussian on a pedestal.
  BeamFitUtils::SpotFit spot = BeamFitUtils::FitSpotFromPoints(h, pts, m);
  out = spot.fit;
  {
    std::lock_guard<std::mutex> lk(g_log_mutex);
    Constants::Detail() << "  beam spot " << AxisTag(sx) << "/" << AxisTag(sy)
                        << " " << tag << ": "
                        << (spot.fit_used ? "fit"
                                          : "CLIPPED MOMENTS (fit failed)")
                        << " sigma=(" << out.sigma_x << "," << out.sigma_y
                        << ") rho=" << out.rho << " chi2/ndf=" << spot.chi2_ndf
                        << std::endl;
  }

  // The per-group beam-gate figures are the only thing written under
  // strip_sum_scatter/<group>; groups outside the plot sample write none.
  if (!Constants::SavePlots()) {
    delete h;
    return out;
  }
  std::lock_guard<std::mutex> lock(g_plot_mutex);
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureAndDraw2DHistogram(h, c);
  // The ellipse PassesEllipse applies.
  BeamFitUtils::DrawEllipse(out, nsigma, kRed + 1);
  // Same folder and naming as the strip gates: beam_gate_s0_s1,
  // beam_gate_s1_grid, ... under strip_sum_scatter/<group>.
  PlottingUtils::SaveFigure(
      c, Form("beam_gate_%s_%s", AxisTag(sx).Data(), AxisTag(sy).Data()),
      subdir, PlotSaveOptions::kLINEAR);
  delete c;

  delete h;
  return out;
}

BeamGate1D StripSumScatter::FindStripGate(TChain *chain, Int_t strip,
                                          const TString &tag,
                                          const TString &subdir,
                                          Double_t nsigma, Double_t center) {
  const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;
  const Double_t kSpotWindow = 4.0;
  const Double_t kFitWindow = 2.0;
  const Int_t kSpotBins = 160;
  const Int_t kMinEvents = 100;
  BeamGate1D out;
  EnergyView ev;
  ev.Attach(chain);
  EnableEventBranches(chain);
  TH1F *h = new TH1F(Form("h1_gate_s%d_%s", strip, tag.Data()),
                     Form(";%s #DeltaE [a.u.];Events", AxisTitle(strip).Data()),
                     C.GATE_BINS, C.GATE_MIN, C.GATE_MAX);
  h->SetDirectory(nullptr);
  const Long64_t n = chain->GetEntries();
  const Long64_t stride = FileSet::SampleStride(n, C.SAMPLE_MAX_POINTS);
  std::vector<Double_t> values;
  values.reserve(C.SAMPLE_MAX_POINTS > 0 ? std::min(n, C.SAMPLE_MAX_POINTS)
                                         : n);
  for (Long64_t j = 0; j < n; j += stride) {
    chain->GetEntry(j);
    ev.Decode();
    const Double_t v = ev.Axis(strip);
    if (!(v > 0.0))
      continue;
    h->Fill(v);
    values.push_back(v);
  }
  if (Int_t(values.size()) < kMinEvents) {
    delete h;
    return out;
  }
  Double_t mean = 0.0;
  const Double_t width = ClippedWidth(values, &mean);
  if (!(width > 0.0)) {
    delete h;
    return out;
  }
  if (center > 0.0)
    mean = center;
  // The clip seeds the fit and is the fallback; the fit is on a fine
  // histogram around the seed, the inner kFitWindow fitted and the rest judged.
  TH1F *hf = new TH1F(Form("%s_fit", h->GetName()), "", kSpotBins,
                      mean - kSpotWindow * width, mean + kSpotWindow * width);
  hf->SetDirectory(nullptr);
  for (Int_t k = 0; k < Int_t(values.size()); k++)
    hf->Fill(values[k]);
  Double_t chi2_ndf = -1.0;
  Double_t fit_mu = mean;
  Double_t fit_sigma = width;
  const Bool_t fit_ok = BeamFitUtils::FitPeak(
      hf, mean, width, kFitWindow, center > 0.0, fit_mu, fit_sigma, &chi2_ndf);
  delete hf;
  out.mu = fit_ok ? fit_mu : mean;
  out.sigma = fit_ok ? fit_sigma : width;
  out.ok = kTRUE;
  {
    std::lock_guard<std::mutex> lk(g_log_mutex);
    Constants::Detail() << "  strip gate s" << strip << " " << tag << ": "
                        << (fit_ok ? "fit" : "CLIPPED WIDTH (fit failed)")
                        << " mu=" << out.mu << " sigma=" << out.sigma
                        << " chi2/ndf=" << chi2_ndf << std::endl;
  }
  if (Constants::SavePlots()) {
    std::lock_guard<std::mutex> lock(g_plot_mutex);
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kTRUE);
    PlottingUtils::ConfigureAndDrawHistogram(h, kBlue + 1);
    const Double_t top = h->GetMaximum();
    TLine *lo = new TLine(out.mu - nsigma * out.sigma, 0.0,
                          out.mu - nsigma * out.sigma, top);
    TLine *hi = new TLine(out.mu + nsigma * out.sigma, 0.0,
                          out.mu + nsigma * out.sigma, top);
    lo->SetLineColor(kRed + 1);
    hi->SetLineColor(kRed + 1);
    lo->SetLineWidth(2);
    hi->SetLineWidth(2);
    lo->Draw();
    hi->Draw();
    PlottingUtils::SaveFigure(c, Form("beam_gate_s%d", strip), subdir);
    delete lo;
    delete hi;
    delete c;
  }
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
  std::vector<TraceClass> classes;
  TraceClass c_aa;
  c_aa.traces = aa;
  c_aa.label = "(#alpha,#alpha')";
  c_aa.color = SimClassColor("aa");
  classes.push_back(c_aa);
  TraceClass c_an;
  c_an.traces = an;
  c_an.label = an_label ? an_label : "(#alpha,n)";
  c_an.color = SimClassColor("an");
  classes.push_back(c_an);
  DrawRegionTraces(save_name, subdir, beam, classes, y_min, y_max, y_title,
                   beam_sigma_measured);
}

Int_t StripSumScatter::SimClassColor(const TString &base) {
  if (base == "aa")
    return kAzure + 2;
  if (base == "an")
    return kRed + 1;
  if (base == "ap")
    return kGreen + 2;
  if (base == "a2n")
    return kMagenta + 2;
  if (base == "ag")
    return kOrange + 7;
  return kGray + 2;
}

// The per-point mean of a trace set with the +-1 RMS band as the y error:
// the shared core of the beam and region mean trace plots. Caller draws it.
static TGraphErrors *TraceMeanBand(const std::vector<TGraph *> &traces) {
  const Int_t npts = traces[0]->GetN();
  std::vector<Double_t> mean(npts, 0.0), m2(npts, 0.0);
  for (Int_t t = 0; t < Int_t(traces.size()); t++) {
    const Double_t *yv = traces[t]->GetY();
    for (Int_t p = 0; p < npts; p++) {
      mean[p] += yv[p];
      m2[p] += yv[p] * yv[p];
    }
  }
  const Double_t nt = Double_t(traces.size());
  const Double_t *xv = traces[0]->GetX();
  TGraphErrors *ge = new TGraphErrors(npts);
  for (Int_t p = 0; p < npts; p++) {
    mean[p] /= nt;
    const Double_t var = m2[p] / nt - mean[p] * mean[p];
    ge->SetPoint(p, xv[p], mean[p]);
    ge->SetPointError(p, 0.0, var > 0.0 ? TMath::Sqrt(var) : 0.0);
  }
  return ge;
}

void StripSumScatter::DrawRegionTraces(const TString &save_name,
                                       const TString &subdir,
                                       const std::vector<TGraph *> &beam,
                                       const std::vector<TraceClass> &classes,
                                       Double_t y_min, Double_t y_max,
                                       const char *y_title,
                                       Bool_t beam_sigma_measured) {
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
    beam_band = TraceMeanBand(beam);
    // Where the strip noise was measured it reads better than the trace
    // spread: swap it in over the RMS.
    if (beam_sigma_measured) {
      for (Int_t p = 0; p < beam_band->GetN(); p++) {
        const Double_t measured =
            StripSigma(Int_t(TMath::Nint(beam_band->GetX()[p]))) *
            beam_band->GetY()[p];
        if (measured > 0.0)
          beam_band->SetPointError(p, 0.0, measured);
      }
    }
    // Dark mean line over a light band.
    beam_band->SetLineColor(kGray + 3);
    beam_band->SetLineWidth(3);
    beam_band->SetFillColorAlpha(kGray + 1, 0.35);
    beam_band->Draw("3 SAME");
  }
  for (size_t k = 0; k < classes.size(); k++)
    DrawTraceSet(classes[k].traces, classes[k].color);
  if (beam_band)
    beam_band->Draw("LX SAME"); // mean line above the traces

  // A legend proxy per class; an empty class is left out: the all-tagged
  // overlay has one class of traces against the beam.
  std::vector<TGraph *> proxies;
  TLegend *leg = PlottingUtils::AddLegend(0.71, 0.88, 0.70, 0.86);
  if (beam_band)
    leg->AddEntry(beam_band, "Beam #pm 1#sigma", "lf");
  for (size_t k = 0; k < classes.size(); k++) {
    if (classes[k].traces.empty())
      continue;
    TGraph *p = new TGraph(1);
    p->SetPoint(0, -1e9, -1e9);
    p->SetLineColor(classes[k].color);
    p->SetLineWidth(3);
    proxies.push_back(p);
    leg->AddEntry(p, classes[k].label, "l");
  }
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
    TGraphErrors *ge = TraceMeanBand(tr);
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

// Build a TGraph trace from a Savitzky-Golay-smoothed set of per-strip totals.
TGraph *StripSumScatter::SmoothedTraceFromTotal(const Double_t *total) {
  Double_t sgd[18];
  SavitzkyGolay(total, sgd);
  return EventsSummary::BuildTraceFromTotals(sgd);
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

  const Double_t kGateNSigma =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA;
  const Double_t kGateCenter =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_CENTER;
  const Int_t kGateBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_BINS;
  const Double_t kGateMin = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_MIN;
  const Double_t kGateMax = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_MAX;

  const Int_t kXBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.XBINS;
  const Int_t kYBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.YBINS;
  const Int_t kGateStrip = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_STRIP;

  // Two parts split by the bar. Before: what decides tagging and keeping; a
  // change there refills. After: plane-only, re-projected from the reservoir.
  TString s =
      Form("v45 reac[%d,%d] short=%d bmult[%d,%d] pileup=%.2fsig,%d "
           "noise=%.2fsig,%d "
           "jump=%.2fsig smooth=%.2fsig "
           "end=%.2fsig gate[s%d,%.2f,c%.4f,%d,%.3f,%.3f]",
           kReacMin, kReacMax, Int_t(Constants::ActiveIgnoreShortStrips()),
           Constants::cfg.STRIP_SUM_SCATTER_CONFIG.BOTH_MULT_MAX,
           Constants::cfg.STRIP_SUM_SCATTER_CONFIG.BOTH_MULT_COUNT_TO,
           Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PILEUP_NSIGMA,
           Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PILEUP_MIN_STRIPS,
           Constants::cfg.STRIP_SUM_SCATTER_CONFIG.NOISE_NSIGMA,
           Constants::cfg.STRIP_SUM_SCATTER_CONFIG.NOISE_MIN_STRIPS,
           kReacJumpNSigma, kSmoothNSigma, kEndStripMax, kGateStrip,
           kGateNSigma, kGateCenter, kGateBins, kGateMin, kGateMax);
  // Which unsegmented strips the all-strips cut requires: a change refills.
  s += Form(
      " req[s0=%d,s17=%d]",
      Int_t(!Constants::cfg.IGNORE_STRIP_0 && Constants::cfg.REQUIRE_STRIP_0),
      Int_t(!Constants::cfg.IGNORE_STRIP_17));
  // The pure-beam entrance ellipse decides the beam sample and its widths.
  s += Form(" one=%d",
            Int_t(Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TAG_RESOLVE));
  s += Form(" pure=%d,%.2fsig",
            Int_t(Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PURE_BEAM_GATE),
            Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PURE_BEAM_NSIGMA);
  // The tail-shape conditions, in sigma, and the measured noise every
  // sigma-unit cut resolves through: a drift in the noise refills.
  {
    const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;
    s += Form(" tail[fall=%d,%.2fsig above=%d,%.2fsig rerise=%.2fsig,%.2fsig "
              "cross=%d cliff=%.2f]",
              C.TAIL_FALL_FROM_STRIP, C.TAIL_RISE_NSIGMA, C.POST_ABOVE_STRIPS,
              C.POST_ABOVE_NSIGMA, C.TAIL_RETURN_NSIGMA, C.TAIL_RERISE_NSIGMA,
              C.POST_CROSS_MIN_STRIP, C.TAIL_CLIFF_MAX_FRACTION);
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
  if (base == "aa")
    return "(#alpha,#alpha')";
  if (base == "an")
    return "(#alpha,n)";
  if (base == "ap")
    return "(#alpha,p)";
  if (base == "a2n")
    return "(#alpha,2n)";
  if (base == "ag")
    return "(#alpha,#gamma)";
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
  // Reference the simulated beam file (the same sim as the overlay
  // populations) so the sim beam lands exactly on 1; else SIM_BEAM_FILE.
  TString file;
  std::vector<RemixSim::SimFileSpec> specs = RemixSim::BuildFileSpecs();
  for (Int_t i = 0; i < Int_t(specs.size()); i++) {
    if (RemixSim::TagWithoutStrip(specs[i].tag) == "beam") {
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
  if (Constants::ActiveIgnoreShortStrips())
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

// The segmented strips where both ends FIRED, read off RAW ADC (sees the
// short end even when IGNORE_SHORT_STRIPS zeroes it), over strips 1..hi.
static Int_t CountBothEnds(const EnergyView &ev, Int_t hi) {
  Int_t nboth = 0;
  for (Int_t s = 1; s <= hi; s++)
    if (ev.leftdE_adc[s - 1] > 0 && ev.rightdE_adc[s - 1] > 0)
      nboth++;
  return nboth;
}

void StripSumScatter::SimTagReport() {
  const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;
  const Int_t kReacMin = C.REACTION_STRIP_MIN;
  const Int_t kReacMax = C.REACTION_STRIP_MAX;
  const Int_t nReac = kReacMax - kReacMin + 1;
  std::vector<RemixSim::SimFileSpec> specs = RemixSim::BuildFileSpecs();
  if (specs.empty())
    return;
  // One file per class and strip (BuildFileSpecs stops on a second); the
  // beam files as a class of their own (strip -1).
  struct SimClass {
    TString base;
    Int_t strip;
    TString file;
  };
  std::vector<SimClass> classes;
  for (Int_t i = 0; i < Int_t(specs.size()); i++) {
    SimClass c;
    c.base = RemixSim::TagWithoutStrip(specs[i].tag);
    c.strip = RemixSim::ReactionStripOf(specs[i].tag);
    c.file = RemixSim::SimRootPath(specs[i]);
    if (c.base != "beam" && (c.strip < kReacMin || c.strip > kReacMax))
      continue;
    classes.push_back(c);
  }
  Double_t gain[18];
  if (!SimBeamGains(gain))
    for (Int_t s = 0; s < 18; s++)
      gain[s] = 1.0;

  std::ostringstream out;
  out << Form("strip-sum-scatter: the data selection over the simulated "
              "populations (%s), against the measured beam spread; the beam "
              "gate is not applied (fitted on data), the rest is the fill's "
              "code.",
              Paths::DatasetName().Data())
      << std::endl;
  out << "Event level: sequential (all strips, pileup, noise, smoothness, "
         "both mult). Own strip: the first failing tag condition, sequential; "
         "tagged = passed them all. One tag: where the event's single tag "
         "went under the one-tag rule (own, an earlier strip, a later strip, "
         "none)"
      << std::endl;
  out << std::endl;
  for (size_t k = 0; k < classes.size(); k++) {
    const SimClass &c = classes[k];
    TFile *fl = IO::OpenForReading(c.file);
    if (!fl || fl->IsZombie()) {
      delete fl;
      out << Form("%s: cannot open %s", c.base.Data(), c.file.Data())
          << std::endl;
      continue;
    }
    TTree *t = static_cast<TTree *>(fl->Get("events_MeV"));
    RemixSim::Event e;
    if (!t || !e.Attach(t)) {
      out << Form("%s: no events_MeV in %s", c.base.Data(), c.file.Data())
          << std::endl;
      fl->Close();
      delete fl;
      continue;
    }
    const Long64_t n = t->GetEntries();
    std::vector<Long64_t> pre(kNPreCuts, 0);
    std::vector<Long64_t> own(kNTagCuts, 0);
    Long64_t tag_own = 0, tag_earlier = 0, tag_later = 0, tag_none = 0;
    std::vector<Long64_t> tag_at(nReac, 0);
    EnergyView ev;
    std::vector<Bool_t> pass(nReac, kFALSE);
    for (Long64_t j = 0; j < n; j++) {
      t->GetEntry(j);
      // The sim's deposit per end, on the data's beam scale (SimTotal's
      // rule), into the view the selection reads.
      for (Int_t s = 1; s <= 16; s++) {
        ev.left[s - 1] = gain[s] * e.Left(s);
        ev.right[s - 1] = gain[s] * e.Right(s);
        if (Constants::ActiveIgnoreShortStrips()) {
          const Double_t tot = gain[s] * e.Total(s);
          ev.left[s - 1] = (s % 2) != 0 ? tot : 0.0;
          ev.right[s - 1] = (s % 2) != 0 ? 0.0 : tot;
        }
        ev.leftdE_adc[s - 1] = ev.left[s - 1] > 0.0 ? 1 : 0;
        ev.rightdE_adc[s - 1] = ev.right[s - 1] > 0.0 ? 1 : 0;
      }
      ev.strip0 = gain[0] * e.strip0;
      ev.strip17 = gain[17] * e.strip17;
      // Event level, sequential as in the fill.
      if (!AllStripsFired(ev)) {
        pre[kPreAllStrips]++;
        continue;
      }
      if (IsPileup(ev)) {
        pre[kPrePileup]++;
        continue;
      }
      if (IsNoise(ev)) {
        pre[kPreNoise]++;
        continue;
      }
      if (!IsSmooth(ev, C.SMOOTHNESS_NSIGMA)) {
        pre[kPreSmooth]++;
        continue;
      }
      if (C.BOTH_MULT_MAX >= 0) {
        if (CountBothEnds(ev, TMath::Min(16, C.BOTH_MULT_COUNT_TO)) >
            C.BOTH_MULT_MAX) {
          pre[kPreBothMult]++;
          continue;
        }
      }
      pre[kPrePass]++;
      for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
        const TagCut why = RejectReason(ev, reac);
        pass[ReacIndex(reac)] = why == kTagPass;
        if (reac == c.strip)
          own[why]++;
      }
      const Int_t won = ResolveTag(ev, pass);
      if (won < 0)
        tag_none++;
      else {
        tag_at[ReacIndex(won)]++;
        if (won == c.strip)
          tag_own++;
        else if (won < c.strip)
          tag_earlier++;
        else
          tag_later++;
      }
    }
    fl->Close();
    delete fl;

    out << Form("%s%s: %lld events (%s)", PrettyLabel(c.base).Data(),
                c.strip >= 0 ? Form(" at strip %d", c.strip) : "", n,
                gSystem->BaseName(c.file))
        << std::endl;
    TString line = "  event level:";
    for (Int_t p = 1; p < kNPreCuts; p++)
      if (p != kPreGate)
        line += Form(" %s %lld,", kPreCutName[p], pre[p]);
    line += Form(" reached the tag %lld (%.1f%%)", pre[kPrePass],
                 n > 0 ? 100.0 * pre[kPrePass] / n : 0.0);
    out << line << std::endl;
    if (c.strip >= 0) {
      line = Form("  own strip %d:", c.strip);
      for (Int_t w = 1; w < kNTagCuts; w++)
        if (w != kCutOtherTag)
          line += Form(" %s %lld,", kTagCutName[w], own[w]);
      line +=
          Form(" tagged %lld (%.1f%% of events, %.1f%% of those reaching "
               "the tag)",
               own[kTagPass], n > 0 ? 100.0 * own[kTagPass] / n : 0.0,
               pre[kPrePass] > 0 ? 100.0 * own[kTagPass] / pre[kPrePass] : 0.0);
      out << line << std::endl;
      out << Form("  one tag: own %lld, earlier strip %lld, later strip %lld, "
                  "none %lld",
                  tag_own, tag_earlier, tag_later, tag_none)
          << std::endl;
    }
    line = "  tagged at:";
    for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
      line += Form(" s%d %lld", reac, tag_at[ReacIndex(reac)]);
    line += Form("  (none %lld)", tag_none);
    out << line << std::endl;
    out << std::endl;
  }
  const std::string text = out.str();
  std::cout << text;
  const TString dir = Paths::ResultsDir() + "/plots/strip_sum_scatter";
  gSystem->mkdir(dir, kTRUE);
  std::ofstream fo((dir + "/sim_tag_report.txt").Data());
  fo << text;
}

namespace {

// One candidate plane, in the shape PlaneXY produces: x sums strips xl..xh,
// y sums yl..yh and may be divided by the event's own upstream mean.
struct StudyPlane {
  Int_t xl, xh, yl, yh;
  Bool_t ratio;
};

// A projected population for one strip and plane: label, colour, marker
// style, and the points as separate x and y lists. The real data sets
// #emphasize so it draws over the simulated clouds rather than blending in.
struct StudySample {
  TString label;
  Int_t color = kBlack;
  Int_t marker = 20; ///< ROOT marker style; the data samples override it.
  Bool_t emphasize = kFALSE;
  std::vector<Double_t> x;
  std::vector<Double_t> y;
};

// A plane's score: the LDA distances of (a,p) from (a,n) and (a,a'), and the
// (a,p) capture at (a,n) leakages of 1, 2 and 5 percent (with (a,a') held
// under 2 percent throughout).
struct StudyRow {
  StudyPlane p;
  Double_t d_an = 0.0, d_aa = 0.0;
  Double_t cap1 = 0.0, cap2 = 0.0, cap5 = 0.0;
  Bool_t current = kFALSE;
};

// The plane coordinates from the cumulative sums. The upstream divisor
// follows PlaneXY: the mean of strips 1..reac-1, divided only when positive.
void StudyXY(const StudyEvt &e, const StudyPlane &p, Int_t reac, Double_t &x,
             Double_t &y) {
  x = e.cs[p.xh] - e.cs[p.xl - 1];
  y = e.cs[p.yh] - e.cs[p.yl - 1];
  if (p.ratio && reac > 1) {
    const Double_t up = e.cs[reac - 1] / Double_t(reac - 1);
    if (up > 0.0)
      y /= up;
  }
}

// The linear separability of two 2-D samples: the Mahalanobis distance
// between their means under the pooled covariance. 0 when either sample is
// too small or the covariance degenerate.
Double_t StudyLda(const std::vector<Double_t> &x1,
                  const std::vector<Double_t> &y1,
                  const std::vector<Double_t> &x2,
                  const std::vector<Double_t> &y2) {
  const Int_t n1 = Int_t(x1.size()), n2 = Int_t(x2.size());
  if (n1 < 30 || n2 < 30)
    return 0.0;
  Double_t m1x = 0.0, m1y = 0.0, m2x = 0.0, m2y = 0.0;
  for (Int_t i = 0; i < n1; i++) {
    m1x += x1[i];
    m1y += y1[i];
  }
  for (Int_t i = 0; i < n2; i++) {
    m2x += x2[i];
    m2y += y2[i];
  }
  m1x /= Double_t(n1);
  m1y /= Double_t(n1);
  m2x /= Double_t(n2);
  m2y /= Double_t(n2);
  Double_t sxx = 0.0, syy = 0.0, sxy = 0.0;
  for (Int_t i = 0; i < n1; i++) {
    const Double_t dx = x1[i] - m1x, dy = y1[i] - m1y;
    sxx += dx * dx;
    syy += dy * dy;
    sxy += dx * dy;
  }
  for (Int_t i = 0; i < n2; i++) {
    const Double_t dx = x2[i] - m2x, dy = y2[i] - m2y;
    sxx += dx * dx;
    syy += dy * dy;
    sxy += dx * dy;
  }
  sxx /= Double_t(n1 + n2 - 2);
  syy /= Double_t(n1 + n2 - 2);
  sxy /= Double_t(n1 + n2 - 2);
  const Double_t det = sxx * syy - sxy * sxy;
  if (!(std::fabs(det) > 1.0e-30))
    return 0.0;
  const Double_t dx = m1x - m2x, dy = m1y - m2y;
  const Double_t d2 =
      (syy * dx * dx - 2.0 * sxy * dx * dy + sxx * dy * dy) / det;
  return d2 > 0.0 ? std::sqrt(d2) : 0.0;
}

// The greedy purity frontier: fill the three samples on a common grid, order
// the cells by signal fraction (Laplace-smoothed), and walk them, reporting
// the largest signal capture reached while both background leakages stayed
// within their budgets. Along that walk the counts only grow, so the walk
// stops the first time either budget is crossed.
Double_t
StudyCapture(const std::vector<Double_t> &sx, const std::vector<Double_t> &sy,
             const std::vector<Double_t> &b1x, const std::vector<Double_t> &b1y,
             const std::vector<Double_t> &b2x, const std::vector<Double_t> &b2y,
             Double_t leak1, Double_t leak2) {
  const Int_t nTot = Int_t(sx.size()) + Int_t(b1x.size()) + Int_t(b2x.size());
  if (Int_t(sx.size()) < 30 || Int_t(b1x.size()) < 30 || nTot < 90)
    return 0.0;
  Double_t xlo = 1.0e30, xhi = -1.0e30, ylo = 1.0e30, yhi = -1.0e30;
  const std::vector<Double_t> *xs[3] = {&sx, &b1x, &b2x};
  const std::vector<Double_t> *ys[3] = {&sy, &b1y, &b2y};
  for (Int_t k = 0; k < 3; k++)
    for (size_t i = 0; i < xs[k]->size(); i++) {
      xlo = TMath::Min(xlo, (*xs[k])[i]);
      xhi = TMath::Max(xhi, (*xs[k])[i]);
      ylo = TMath::Min(ylo, (*ys[k])[i]);
      yhi = TMath::Max(yhi, (*ys[k])[i]);
    }
  if (!(xhi > xlo) || !(yhi > ylo))
    return 0.0;
  const Int_t kB = 72;
  std::vector<Double_t> ns(kB * kB, 0.0), nb1(kB * kB, 0.0), nb2(kB * kB, 0.0);
  for (Int_t k = 0; k < 3; k++) {
    std::vector<Double_t> &cnt = k == 0 ? ns : (k == 1 ? nb1 : nb2);
    for (size_t i = 0; i < xs[k]->size(); i++) {
      Int_t bx = Int_t(((*xs[k])[i] - xlo) / (xhi - xlo) * Double_t(kB));
      Int_t by = Int_t(((*ys[k])[i] - ylo) / (yhi - ylo) * Double_t(kB));
      bx = TMath::Min(TMath::Max(bx, 0), kB - 1);
      by = TMath::Min(TMath::Max(by, 0), kB - 1);
      cnt[by * kB + bx] += 1.0;
    }
  }
  std::vector<std::pair<Double_t, Int_t>> order;
  for (Int_t c = 0; c < kB * kB; c++) {
    const Double_t tot = ns[c] + nb1[c] + nb2[c];
    if (tot <= 0.0)
      continue;
    order.push_back(std::make_pair((ns[c] + 0.1) / (tot + 0.4), c));
  }
  std::sort(
      order.begin(), order.end(),
      [](const std::pair<Double_t, Int_t> &a,
         const std::pair<Double_t, Int_t> &b) { return a.first > b.first; });
  const Double_t nS = Double_t(sx.size());
  const Double_t nB1 = Double_t(b1x.size());
  const Double_t nB2 = Double_t(b2x.size());
  Double_t cap = 0.0, g1 = 0.0, g2 = 0.0;
  for (size_t o = 0; o < order.size(); o++) {
    const Int_t c = order[o].second;
    const Double_t next1 = g1 + nb1[c] / nB1;
    const Double_t next2 = g2 + nb2[c] / nB2;
    if (next1 > leak1 || next2 > leak2)
      break;
    g1 = next1;
    g2 = next2;
    cap += ns[c] / nS;
  }
  return cap;
}

// Report sort: by the capture at 2 percent (a,n) leakage.
Bool_t StudyRowLess(const StudyRow &a, const StudyRow &b) {
  return a.cap2 > b.cap2;
}

// One study figure: the four projected populations as translucent point
// clouds on a frame spanning them all, with axes named for the plane.
void StudyFigure(const TString &name, const TString &subtitle,
                 const std::vector<StudySample> &samples) {
  Double_t xlo = 1.0e30, xhi = -1.0e30, ylo = 1.0e30, yhi = -1.0e30;
  for (size_t k = 0; k < samples.size(); k++)
    for (size_t i = 0; i < samples[k].x.size(); i++) {
      xlo = TMath::Min(xlo, samples[k].x[i]);
      xhi = TMath::Max(xhi, samples[k].x[i]);
      ylo = TMath::Min(ylo, samples[k].y[i]);
      yhi = TMath::Max(yhi, samples[k].y[i]);
    }
  if (!(xhi > xlo) || !(yhi > ylo))
    return;
  const Double_t mx = 0.05 * (xhi - xlo), my = 0.05 * (yhi - ylo);
  TH2F frame("study_frame", subtitle, 10, xlo - mx, xhi + mx, 10, ylo - my,
             yhi + my);
  frame.SetStats(0);
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  c->SetLeftMargin(0.18);
  frame.Draw();
  // Bottom-right: the beam sits bottom-left and the reactions top, leaving
  // this corner clear. Tall enough for both data entries too.
  TLegend *leg = PlottingUtils::AddLegend(0.60, 0.88, 0.09, 0.50);
  for (size_t k = 0; k < samples.size(); k++) {
    const StudySample &s = samples[k];
    if (s.x.empty())
      continue;
    TGraph *g = new TGraph();
    const Long64_t stride = FileSet::SampleStride(Int_t(s.x.size()), 4000);
    Int_t n = 0;
    for (size_t i = 0; i < s.x.size(); i += Long64_t(stride))
      g->SetPoint(n++, s.x[i], s.y[i]);
    g->Set(n);
    g->SetMarkerStyle(s.marker);
    g->SetMarkerSize(s.emphasize ? 0.55 : 0.3);
    g->SetMarkerColorAlpha(s.color, s.emphasize ? 0.8 : 0.35);
    g->SetLineColor(s.color);
    g->Draw("P SAME");
    leg->AddEntry(g, s.label, "p");
  }
  leg->Draw();
  PlottingUtils::SaveFigure(c, name, "sim_scatter", PlotSaveOptions::kLINEAR);
  delete leg;
  delete c;
}

// The one sim beam file (the first "beam" spec, else SIM_BEAM_FILE),
// opened; null when it cannot be. Callers close and delete it.
TFile *StudyBeamFile() {
  TString file;
  std::vector<RemixSim::SimFileSpec> specs = RemixSim::BuildFileSpecs();
  for (Int_t i = 0; i < Int_t(specs.size()); i++) {
    if (RemixSim::TagWithoutStrip(specs[i].tag) == "beam") {
      file = RemixSim::SimRootPath(specs[i]);
      break;
    }
  }
  if (file.Length() == 0)
    file =
        Paths::DatasetDir() + "/sim_root_files/" + Constants::cfg.SIM_BEAM_FILE;
  TFile *f = IO::OpenForReading(file);
  if (!f || f->IsZombie()) {
    if (f)
      delete f;
    return nullptr;
  }
  return f;
}

// Mean and rms of the gain-normalized per-strip sim beam total, over the
// same >0 events SimBeamGains averages -- so mean[s] is 1 wherever a gain
// was set, and rms[s] is the sim's own relative resolution on that strip.
Bool_t StudySimBeamMoments(const Double_t *gain, Double_t *mean,
                           Double_t *rms) {
  for (Int_t s = 0; s < 18; s++) {
    mean[s] = 0.0;
    rms[s] = 0.0;
  }
  TFile *f = StudyBeamFile();
  if (!f)
    return kFALSE;
  TTree *t = static_cast<TTree *>(f->Get("events_MeV"));
  RemixSim::Event e;
  if (!t || !e.Attach(t)) {
    f->Close();
    delete f;
    return kFALSE;
  }
  const Long64_t n = t->GetEntries();
  const Long64_t stride = FileSet::SampleStride(
      n, Constants::cfg.STRIP_SUM_SCATTER_CONFIG.SAMPLE_MAX_POINTS);
  Double_t sum[18] = {0}, sum2[18] = {0};
  Long64_t cnt[18] = {0};
  for (Long64_t j = 0; j < n; j += stride) {
    t->GetEntry(j);
    for (Int_t s = 0; s < 18; s++) {
      const Double_t tot = gain[s] * e.Total(s);
      if (tot > 0.0) {
        sum[s] += tot;
        sum2[s] += tot * tot;
        cnt[s]++;
      }
    }
  }
  f->Close();
  delete f;
  for (Int_t s = 0; s < 18; s++) {
    if (cnt[s] < 2 || !(sum[s] > 0.0))
      continue;
    mean[s] = sum[s] / Double_t(cnt[s]);
    const Double_t v = sum2[s] / Double_t(cnt[s]) - mean[s] * mean[s];
    rms[s] = v > 0.0 ? std::sqrt(v) : 0.0;
  }
  return kTRUE;
}

// The relative rms (rms / mean) of a sample of event totals.
Double_t StudyRelRms(const std::vector<Double_t> &v) {
  if (v.size() < 2)
    return 0.0;
  Double_t sum = 0.0, sum2 = 0.0;
  for (size_t i = 0; i < v.size(); i++) {
    sum += v[i];
    sum2 += v[i] * v[i];
  }
  const Double_t m = sum / Double_t(v.size());
  const Double_t s2 = sum2 / Double_t(v.size()) - m * m;
  return (s2 > 0.0 && m > 0.0) ? std::sqrt(s2) / m : 0.0;
}

// Relative rms of the sim beam's full total (strips 1..16) under one
// smearing setting: each strip its multiplicative jitter, the event its
// common flux factor. A separate rng keeps the calibration reproducible.
Double_t StudyBeamTotalRms(const Double_t *gain, const Double_t *jit,
                           Double_t sigmaFlux, ULong64_t seed) {
  TFile *f = StudyBeamFile();
  if (!f)
    return 0.0;
  TTree *t = static_cast<TTree *>(f->Get("events_MeV"));
  RemixSim::Event e;
  if (!t || !e.Attach(t)) {
    f->Close();
    delete f;
    return 0.0;
  }
  TRandom3 rng(seed);
  const Long64_t n = t->GetEntries();
  const Long64_t stride = FileSet::SampleStride(
      n, Constants::cfg.STRIP_SUM_SCATTER_CONFIG.SAMPLE_MAX_POINTS);
  std::vector<Double_t> tots;
  for (Long64_t j = 0; j < n; j += stride) {
    t->GetEntry(j);
    const Double_t g = sigmaFlux > 0.0 ? 1.0 + rng.Gaus(0.0, sigmaFlux) : 1.0;
    Double_t tot = 0.0;
    for (Int_t s = 1; s <= 16; s++) {
      Double_t f_s = g;
      if (jit[s] > 0.0)
        f_s *= 1.0 + rng.Gaus(0.0, jit[s]);
      tot += f_s * gain[s] * e.Total(s);
    }
    if (tot > 0.0)
      tots.push_back(tot);
  }
  f->Close();
  delete f;
  return StudyRelRms(tots);
}

// The full total (strips 1..16) of every pure-beam data event, the sample
// the sim beam total is matched against.
std::vector<Double_t>
StudyReservoirBeamTotals(const std::vector<TraceEvt> &res) {
  std::vector<Double_t> tots;
  for (size_t k = 0; k < res.size(); k++) {
    const TraceEvt &e = res[k];
    if (!e.beam_flat)
      continue;
    Double_t tot = 0.0;
    for (Int_t s = 1; s <= 16; s++)
      tot += e.Total(s);
    if (tot > 0.0)
      tots.push_back(tot);
  }
  return tots;
}

// The data tagged at one strip, as study events, straight from the
// reservoir: the fit's data sample.
std::vector<StudyEvt> StudyDataEvts(Int_t reac,
                                    const std::vector<TraceEvt> &res) {
  const Int_t ri =
      reac - Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  std::vector<StudyEvt> out;
  for (size_t k = 0; k < res.size(); k++) {
    const TraceEvt &e = res[k];
    if (!(e.reac_mask & (1u << ri)))
      continue;
    Double_t tot[18];
    e.Totals(tot);
    StudyEvt se;
    se.cls = -1;
    se.src = reac;
    se.tagged = reac;
    se.cs[0] = 0.0;
    for (Int_t s = 1; s <= 17; s++)
      se.cs[s] = se.cs[s - 1] + tot[s];
    out.push_back(se);
  }
  return out;
}

// The smearing shared by the plane study and the template fit. A positive
// jitter is one value for every strip; a negative one auto-calibrates:
// per-strip jitter from the strip widths, then a per-event flux jitter from
// the beam total's spread, taken out of the per-strip budget so each strip
// stays matched. The applied jitters land in jit and sigmaFlux.
void StudySmearingCalibration(const Double_t *gain, Double_t jitter,
                              const std::vector<Double_t> &dataBeamTotals,
                              Double_t *jit, Double_t &sigmaFlux,
                              std::ostringstream &calib) {
  for (Int_t s = 0; s < 18; s++)
    jit[s] = jitter > 0.0 ? jitter : 0.0;
  sigmaFlux = 0.0;
  if (!(jitter < 0.0))
    return;
  Double_t simMu[18], simRms[18];
  if (!StudySimBeamMoments(gain, simMu, simRms)) {
    std::cerr << "strip-sum-scatter: auto jitter cannot read the sim beam "
                 "file; the run proceeds unsmeared."
              << std::endl;
    return;
  }
  Double_t jit0[18];
  calib << "  Resolution match: multiplicative smearing added to the "
           "simulated traces (in quadrature on the sim's own width): a "
           "per-strip jitter so each strip's relative beam rms reaches the "
           "measured one, plus one per-event flux jitter so the beam "
           "total's spread reaches the measured one."
        << std::endl;
  calib << Form("  %6s  %11s  %11s  %8s", "strip", "sim rel-rms",
                "data rel-rms", "jitter")
        << std::endl;
  for (Int_t s = 0; s < 18; s++) {
    jit0[s] = 0.0;
    const Double_t dataMu = StripSumScatter::StripMean(s),
                   dataSg = StripSumScatter::StripSigma(s);
    if (!(dataMu > 0.0) || !(dataSg > 0.0) || !(simMu[s] > 0.0))
      continue;
    const Double_t wd = dataSg / dataMu, ws = simRms[s] / simMu[s];
    const Double_t diff = wd * wd - ws * ws;
    jit0[s] = diff > 0.0 ? std::sqrt(diff) : 0.0;
  }
  const Double_t dataTot = StudyRelRms(dataBeamTotals);
  const Double_t simTot0 = StudyBeamTotalRms(gain, jit0, 0.0, 20260930);
  const Double_t flux2 = dataTot * dataTot - simTot0 * simTot0;
  sigmaFlux = flux2 > 0.0 ? std::sqrt(flux2) : 0.0;
  for (Int_t s = 0; s < 18; s++) {
    const Double_t d = jit0[s] * jit0[s] - sigmaFlux * sigmaFlux;
    jit[s] = d > 0.0 ? std::sqrt(d) : 0.0;
  }
  const Double_t simTotFin = StudyBeamTotalRms(gain, jit, sigmaFlux, 20260931);
  for (Int_t s = 0; s < 18; s++) {
    const Double_t dataMu = StripSumScatter::StripMean(s),
                   dataSg = StripSumScatter::StripSigma(s);
    if (!(dataMu > 0.0) || !(dataSg > 0.0) || !(simMu[s] > 0.0))
      continue;
    calib << Form("  %6d  %11.4f  %11.4f  %8.4f", s, simRms[s] / simMu[s],
                  dataSg / dataMu, jit[s])
          << std::endl;
  }
  calib << Form("  beam total rel-rms: data %.4f, sim %.4f before and "
                "%.4f after, with %.4f per-event flux jitter (in addition "
                "to the per-strip jitter above).",
                dataTot, simTot0, simTotFin, sigmaFlux)
        << std::endl;
  calib << std::endl;
}

// --- Template fit: the extended binned Poisson likelihood of unit-summed
// templates, maximized by EM. The yields update multiplicatively and the
// likelihood is concave in them, so EM lands on the maximum likelihood.
struct StudyFitResult {
  std::vector<Double_t> N;   ///< Yields at the fit
  std::vector<Double_t> cov; ///< Covariance, K*K row-major
  Double_t dev = 0.0;        ///< -2 ln L at the fit
  Int_t iters = 0;           ///< EM iterations run
  Bool_t converged = kFALSE; ///< Yields stopped moving
  Bool_t singular = kFALSE;  ///< The covariance would not invert
};

StudyFitResult StudyFitEM(Int_t nb, const std::vector<Double_t> &d,
                          const std::vector<std::vector<Double_t>> &P) {
  const Int_t K = Int_t(P.size());
  std::vector<Double_t> N(K, 0.0), mu(nb, 0.0);
  Double_t nTot = 0.0;
  for (Int_t b = 0; b < nb; b++)
    nTot += d[b];
  for (Int_t k = 0; k < K; k++)
    N[k] = nTot / Double_t(K);
  StudyFitResult r;
  for (Int_t it = 0; it < 50000; it++) {
    r.iters = it + 1;
    for (Int_t b = 0; b < nb; b++) {
      Double_t m = 0.0;
      for (Int_t k = 0; k < K; k++)
        m += N[k] * P[k][b];
      mu[b] = m > 0.0 ? m : 1.0e-9;
    }
    std::vector<Double_t> Nn(K, 0.0);
    for (Int_t k = 0; k < K; k++) {
      Double_t u = 0.0;
      for (Int_t b = 0; b < nb; b++)
        u += d[b] * P[k][b] / mu[b];
      Nn[k] = N[k] * u;
    }
    Double_t maxrel = 0.0;
    for (Int_t k = 0; k < K; k++)
      maxrel = TMath::Max(maxrel,
                          std::fabs(Nn[k] - N[k]) / TMath::Max(Nn[k], 1.0e-6));
    N = Nn;
    r.converged = maxrel < 1.0e-7;
    if (r.converged)
      break;
  }
  r.N = N;
  r.dev = 0.0;
  for (Int_t b = 0; b < nb; b++)
    r.dev += 2.0 * (mu[b] - d[b] +
                    (d[b] > 0.0 ? d[b] * std::log(d[b] / mu[b]) : 0.0));
  // The Hessian of -ln L in the yields, inverted for the covariance.
  TMatrixD H(K, K);
  for (Int_t b = 0; b < nb; b++)
    for (Int_t i = 0; i < K; i++)
      for (Int_t j = i; j < K; j++) {
        const Double_t h = P[i][b] * P[j][b] / mu[b];
        H(i, j) += h;
        if (j != i)
          H(j, i) += h;
      }
  r.cov.assign(K * K, 0.0);
  TDecompLU de(H);
  r.singular = !de.Decompose();
  if (!r.singular) {
    const TMatrixD C = de.Invert();
    for (Int_t i = 0; i < K; i++)
      for (Int_t j = 0; j < K; j++)
        r.cov[i * K + j] = C(i, j);
  }
  return r;
}

// The q-quantile of a sample, linearly interpolated.
Double_t StudyQuantile(std::vector<Double_t> v, Double_t q) {
  if (v.empty())
    return 0.0;
  std::sort(v.begin(), v.end());
  const Double_t idx = q * Double_t(Int_t(v.size()) - 1);
  const Int_t lo = Int_t(idx);
  const Int_t hi = TMath::Min(lo + 1, Int_t(v.size()) - 1);
  const Double_t f = idx - Double_t(lo);
  return v[lo] * (1.0 - f) + v[hi] * f;
}

// The mean and standard deviation of a sample.
void StudyMoments(const std::vector<Double_t> &v, Double_t &mean,
                  Double_t &sd) {
  mean = 0.0;
  sd = 0.0;
  if (v.size() < 2)
    return;
  Double_t s = 0.0, s2 = 0.0;
  for (size_t i = 0; i < v.size(); i++) {
    s += v[i];
    s2 += v[i] * v[i];
  }
  mean = s / Double_t(v.size());
  const Double_t v2 = s2 / Double_t(v.size()) - mean * mean;
  sd = v2 > 0.0 ? std::sqrt(v2) : 0.0;
}

// One fit figure: the 1-D y distribution with the fitted components and
// the data, and the pulls (d - mu) / sqrt(mu) below it.
void StudyFitFigure(const TString &name, const TString &title, Int_t nb,
                    const Double_t *ye, const std::vector<Double_t> &d,
                    const std::vector<Double_t> &mu,
                    const std::vector<TString> &labels,
                    const std::vector<Int_t> &colors,
                    const std::vector<Double_t> &N,
                    const std::vector<Double_t> &E,
                    const std::vector<std::vector<Double_t>> &P) {
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  c->Divide(1, 2, 0.0, 0.0, 0.0);
  c->cd(1);
  // The frame is filled with the fit total so its auto y range covers the
  // components; an empty frame would clamp the axis and clip every template.
  TH1F frame("fit_frame", title, nb, ye);
  frame.SetStats(0);
  frame.GetYaxis()->SetTitle("events / bin");
  std::vector<TH1F *> comps(N.size());
  for (size_t k = 0; k < N.size(); k++) {
    TH1F *h = new TH1F(Form("fit_c%d", Int_t(k)), labels[k].Data(), nb, ye);
    for (Int_t b = 0; b < nb; b++)
      h->SetBinContent(b + 1, P[k][b] * N[k]);
    h->SetLineColor(colors[k]);
    h->SetLineWidth(1.5);
    comps[k] = h;
  }
  TH1F *tot = new TH1F("fit_tot", "fit", nb, ye);
  for (Int_t b = 0; b < nb; b++) {
    Double_t s = 0.0;
    for (size_t k = 0; k < N.size(); k++)
      s += comps[k]->GetBinContent(b + 1);
    tot->SetBinContent(b + 1, s);
  }
  Double_t ymax = 0.0;
  for (Int_t b = 0; b < nb; b++) {
    if (tot->GetBinContent(b + 1) > ymax)
      ymax = tot->GetBinContent(b + 1);
    if (d[b] > ymax)
      ymax = d[b];
  }
  frame.SetMaximum(1.15 * ymax + 1.0);
  frame.SetMinimum(0.0);
  frame.Draw();
  for (size_t k = 0; k < comps.size(); k++)
    comps[k]->Draw("hist SAME");
  tot->SetLineColor(kBlack);
  tot->SetLineWidth(2);
  tot->Draw("hist SAME");
  TGraphErrors *gd = new TGraphErrors(nb);
  for (Int_t b = 0; b < nb; b++) {
    gd->SetPoint(b, 0.5 * (ye[b] + ye[b + 1]), d[b]);
    gd->SetPointError(b, 0.0, std::sqrt(d[b] > 0.0 ? d[b] : 1.0));
  }
  gd->SetMarkerStyle(20);
  gd->SetMarkerSize(0.5);
  gd->SetMarkerColor(kBlack);
  gd->Draw("P SAME");
  TLegend *leg = PlottingUtils::AddLegend(0.12, 0.52, 0.50, 0.90);
  leg->AddEntry(tot, "fit", "l");
  for (size_t k = 0; k < N.size(); k++)
    leg->AddEntry(comps[k],
                  Form("%s  %.0f #plusmn %.0f", labels[k].Data(), N[k], E[k]),
                  "l");
  leg->AddEntry(gd, "data", "p");
  leg->Draw();
  // Bottom pad: the pulls.
  c->cd(2);
  TGraphErrors *gp = new TGraphErrors(nb);
  for (Int_t b = 0; b < nb; b++) {
    const Double_t dd = d[b], m = mu[b];
    const Double_t pull =
        m > 1.0e-9 ? (dd - m) / std::sqrt(m) : (dd > 0.0 ? 6.0 : 0.0);
    gp->SetPoint(b, 0.5 * (ye[b] + ye[b + 1]), pull);
    gp->SetPointError(b, 0.0, 1.0);
  }
  gp->SetMarkerStyle(20);
  gp->SetMarkerSize(0.4);
  gp->SetMarkerColor(kBlack);
  gp->Draw("AP");
  gp->GetYaxis()->SetRangeUser(-6.0, 6.0);
  gp->GetYaxis()->SetTitle("pull");
  gp->GetXaxis()->SetTitle("y");
  PlottingUtils::SaveFigure(c, name, "sim_scatter", PlotSaveOptions::kLINEAR);
  delete gp;
  delete gd;
  delete tot;
  for (size_t k = 0; k < comps.size(); k++)
    delete comps[k];
  delete leg;
  delete c;
}

// Fill the figure samples: the four simulated populations, then the real
// beam and the real data tagged at this strip, so both draw over the clouds.
void StudyMakeSamples(Int_t reac, const std::vector<const StudyEvt *> *pop,
                      const std::vector<StudyEvt> &dataBeam,
                      const std::vector<StudyEvt> &data, const StudyPlane &p,
                      std::vector<StudySample> &samples) {
  samples.assign(6, StudySample());
  samples[0].label = "Beam";
  samples[0].color = kGray + 2;
  samples[1].label = "(#alpha,#alpha')";
  samples[1].color = kAzure + 2;
  samples[2].label = "(#alpha,n)";
  samples[2].color = kRed + 1;
  samples[3].label = "(#alpha,p)";
  samples[3].color = kGreen + 2;
  samples[4].label = "Data beam";
  samples[4].color = kBlack;
  samples[4].marker = 26;
  samples[4].emphasize = kTRUE;
  samples[5].label = Form("Data tagged at %d", reac);
  samples[5].color = kBlack;
  samples[5].marker = 22;
  samples[5].emphasize = kTRUE;
  for (Int_t k = 0; k < 4; k++)
    for (size_t i = 0; i < pop[k].size(); i++) {
      Double_t x = 0.0, y = 0.0;
      StudyXY(*pop[k][i], p, reac, x, y);
      samples[k].x.push_back(x);
      samples[k].y.push_back(y);
    }
  for (size_t i = 0; i < dataBeam.size(); i++) {
    Double_t x = 0.0, y = 0.0;
    StudyXY(dataBeam[i], p, reac, x, y);
    samples[4].x.push_back(x);
    samples[4].y.push_back(y);
  }
  for (size_t i = 0; i < data.size(); i++) {
    Double_t x = 0.0, y = 0.0;
    StudyXY(data[i], p, reac, x, y);
    samples[5].x.push_back(x);
    samples[5].y.push_back(y);
  }
}

} // namespace

std::vector<StudyEvt> StripSumScatter::StudyCollect(const Double_t *gain,
                                                    const Double_t *jit,
                                                    Double_t sigmaFlux,
                                                    ULong64_t seed) {
  const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;
  const Int_t kReacMin = C.REACTION_STRIP_MIN;
  const Int_t kReacMax = C.REACTION_STRIP_MAX;
  const Int_t nReac = kReacMax - kReacMin + 1;
  const TagThresholds T = NominalThresholds();
  // The beam gate was fitted on data (Prepare's first group); a failed fit
  // or a grid-axis gate cannot be applied to the sim, which has no grid.
  const Bool_t applyGate = T.gate_nsigma > 0.0 && m_simStudyGate.ok &&
                           C.GATE_STRIP != GATE_AXIS_GRID;
  std::vector<RemixSim::SimFileSpec> specs = RemixSim::BuildFileSpecs();
  if (specs.empty()) {
    std::cerr << "strip-sum-scatter: no sim control files; no simulated "
                 "events to collect."
              << std::endl;
    return std::vector<StudyEvt>();
  }
  TRandom3 rng(seed);
  std::vector<StudyEvt> evts;
  for (size_t k = 0; k < specs.size(); k++) {
    const TString base = RemixSim::TagWithoutStrip(specs[k].tag);
    Int_t cls = -1;
    if (base == "beam")
      cls = 0;
    else if (base == "aa")
      cls = 1;
    else if (base == "an")
      cls = 2;
    else if (base == "ap")
      cls = 3;
    if (cls < 0)
      continue;
    const Int_t src = cls == 0 ? -1 : RemixSim::ReactionStripOf(specs[k].tag);
    if (cls > 0 && (src < kReacMin || src > kReacMax))
      continue;
    TFile *fl = IO::OpenForReading(RemixSim::SimRootPath(specs[k]));
    if (!fl || fl->IsZombie()) {
      delete fl;
      std::cerr << "  cannot open the sim file for " << specs[k].tag << ": "
                << RemixSim::SimRootPath(specs[k]) << std::endl;
      continue;
    }
    TTree *t = static_cast<TTree *>(fl->Get("events_MeV"));
    RemixSim::Event e;
    if (!t || !e.Attach(t)) {
      std::cerr << "  no events_MeV for " << specs[k].tag << std::endl;
      if (fl)
        fl->Close();
      delete fl;
      continue;
    }
    const Long64_t n = t->GetEntries();
    Long64_t kept = 0, tagged = 0;
    EnergyView ev;
    std::vector<Bool_t> pass(nReac, kFALSE);
    for (Long64_t j = 0; j < n; j++) {
      t->GetEntry(j);
      // The sim's deposit per end on the data's beam scale, exactly as
      // SimTagReport builds it, with the smearing: one factor per event
      // (the flux) times one per strip, so the strip total carries both.
      const Double_t g = sigmaFlux > 0.0 ? 1.0 + rng.Gaus(0.0, sigmaFlux) : 1.0;
      for (Int_t s = 1; s <= 16; s++) {
        Double_t f = g;
        if (jit[s] > 0.0)
          f *= 1.0 + rng.Gaus(0.0, jit[s]);
        const Double_t tot = f * gain[s] * e.Total(s);
        if (Constants::ActiveIgnoreShortStrips()) {
          ev.left[s - 1] = (s % 2) != 0 ? tot : 0.0;
          ev.right[s - 1] = (s % 2) != 0 ? 0.0 : tot;
        } else {
          ev.left[s - 1] = Float_t(f * gain[s] * e.Left(s));
          ev.right[s - 1] = Float_t(f * gain[s] * e.Right(s));
        }
        ev.leftdE_adc[s - 1] = ev.left[s - 1] > 0.0 ? 1 : 0;
        ev.rightdE_adc[s - 1] = ev.right[s - 1] > 0.0 ? 1 : 0;
      }
      ev.strip0 =
          Float_t(g * (jit[0] > 0.0 ? 1.0 + rng.Gaus(0.0, jit[0]) : 1.0) *
                  gain[0] * e.strip0);
      ev.strip17 =
          Float_t(g * (jit[17] > 0.0 ? 1.0 + rng.Gaus(0.0, jit[17]) : 1.0) *
                  gain[17] * e.strip17);
      // The event-level funnel, in the fill's order.
      if (!AllStripsFired(ev))
        continue;
      if (applyGate &&
          !PassesGate(m_simStudyGate, ev, C.GATE_STRIP, T.gate_nsigma))
        continue;
      if (IsPileup(ev, T.pileup_nsigma))
        continue;
      if (IsNoise(ev, T.noise_nsigma))
        continue;
      if (T.smooth_nsigma > 0.0 && MaxStepNSigma(ev) > T.smooth_nsigma)
        continue;
      if (C.BOTH_MULT_MAX >= 0 &&
          CountBothEnds(ev, TMath::Min(16, C.BOTH_MULT_COUNT_TO)) >
              C.BOTH_MULT_MAX)
        continue;
      // The tag at every strip and the one-tag rule, nominal thresholds.
      for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
        pass[ReacIndex(reac)] = RejectReason(ev, reac) == kTagPass;
      StudyEvt se;
      se.cls = cls;
      se.src = src;
      se.tagged = ResolveTag(ev, pass);
      se.cs[0] = 0.0;
      for (Int_t s = 1; s <= 17; s++)
        se.cs[s] = se.cs[s - 1] + ev.Total(s);
      evts.push_back(se);
      kept++;
      if (se.tagged >= 0)
        tagged++;
    }
    fl->Close();
    delete fl;
    std::cout << Form("  %s at %s: %lld events, %lld past the event-level "
                      "cuts, %lld tagged somewhere",
                      base.Data(), src >= 0 ? Form("strip %d", src) : "beam", n,
                      kept, tagged)
              << std::endl;
  }
  return evts;
}

void StripSumScatter::PlaneStudy() {
  const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;
  const Int_t kReacMin = C.REACTION_STRIP_MIN;
  const Int_t kReacMax = C.REACTION_STRIP_MAX;
  const Int_t nReac = kReacMax - kReacMin + 1;
  const TagThresholds T = NominalThresholds();
  std::cout << "strip-sum-scatter: plane study, scoring candidate scatter "
               "planes on the simulated populations."
            << std::endl;

  // The beam gate was fitted on data (Prepare's first group); a failed fit
  // or a grid-axis gate cannot be applied to the sim, which has no grid.
  const Bool_t applyGate = T.gate_nsigma > 0.0 && m_simStudyGate.ok &&
                           C.GATE_STRIP != GATE_AXIS_GRID;
  if (T.gate_nsigma > 0.0 && !applyGate)
    std::cout << "  note: the beam gate is skipped in the study (the sim has "
                 "no grid, or the gate was not fitted)."
              << std::endl;

  Double_t gain[18];
  if (!SimBeamGains(gain))
    for (Int_t s = 0; s < 18; s++)
      gain[s] = 1.0;
  // Smearing of the simulated traces: PLANE_STUDY_EXTRA_JITTER > 0 is one
  // Gaussian per strip; < 0 auto-calibrates the per-strip jitter and a
  // per-event flux jitter to the measured beam (the calibration table goes
  // into the report).
  const Double_t jitter = C.PLANE_STUDY_EXTRA_JITTER;
  const Bool_t autoJit = jitter < 0.0;
  Double_t jit[18];
  Double_t sigmaFlux = 0.0;
  std::ostringstream calib;
  StudySmearingCalibration(gain, jitter, StudyReservoirBeamTotals(m_reservoir),
                           jit, sigmaFlux, calib);

  // --- Collect: every sim population through the fill's full selection.
  std::vector<StudyEvt> evts = StudyCollect(gain, jit, sigmaFlux, 20260929);
  if (evts.empty())
    return;

  // Pure-beam data events from the reservoir, projected exactly like the
  // sim beam, so the figures can compare the measured and simulated widths.
  std::vector<StudyEvt> dataBeamEvts;
  for (size_t k = 0; k < m_reservoir.size(); k++) {
    const TraceEvt &e = m_reservoir[k];
    if (!e.beam_flat)
      continue;
    Double_t tot[18];
    e.Totals(tot);
    StudyEvt se;
    se.cls = -1;
    se.src = -1;
    se.tagged = -1;
    se.cs[0] = 0.0;
    for (Int_t s = 1; s <= 17; s++)
      se.cs[s] = se.cs[s - 1] + tot[s];
    dataBeamEvts.push_back(se);
  }

  // --- Score: candidate planes per reaction strip.
  std::ostringstream out;
  out << "strip-sum-scatter: plane study (" << Paths::DatasetName().Data()
      << "): simulated populations pushed through the fill's full nominal "
         "selection (event-level cuts"
      << (applyGate ? ", beam gate," : ", no beam gate,")
      << " tag, one-tag rule), projected onto candidate planes x = sum of "
         "strips 1..xh, y = sum of reac+1..yh (/event upstream mean when "
         "ratio); smearing: "
      << (autoJit ? "auto-calibrated per strip and per event (flux) to the "
                    "measured beam"
                  : Form("extra jitter %.1f%%", 100.0 * jitter))
      << "." << std::endl;
  if (autoJit)
    out << calib.str();
  out << "Samples per strip: (a,p), (a,n) and (a,a') GENERATED at and TAGGED "
         "at that strip. d(ap|an), d(ap|aa): LDA distance between the clouds. "
         "cap@k: (a,p) capture at (a,n) leakage k (with (a,a') leakage kept "
         "under 2%), from a cell-greedy purity frontier: an optimistic "
         "ceiling, not a drawable cut. Migration is deliberately excluded: "
         "events tagged away from their birth strip are another strip's "
         "problem (see sim_tag_report.txt)."
      << std::endl;
  out << "Figures also carry the pure-beam data from the reservoir ("
      << dataBeamEvts.size() << " events) next to the simulated beam."
      << std::endl;
  out << std::endl;

  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    // The four populations for this strip: beam for the figures only, aa/an
    // as backgrounds, ap as signal.
    std::vector<const StudyEvt *> pop[4];
    for (size_t i = 0; i < evts.size(); i++) {
      const StudyEvt &e = evts[i];
      if (e.cls == 0)
        pop[0].push_back(&e);
      else if (e.src == reac && e.tagged == reac) {
        if (e.cls == 1)
          pop[1].push_back(&e);
        else if (e.cls == 2)
          pop[2].push_back(&e);
        else if (e.cls == 3)
          pop[3].push_back(&e);
      }
    }
    // The real data tagged at this strip, straight from the reservoir: the
    // figures carry it so the simulated clouds can be read against it.
    std::vector<StudyEvt> dataEvts;
    for (size_t k = 0; k < m_reservoir.size(); k++) {
      const TraceEvt &e = m_reservoir[k];
      if (!(e.reac_mask & (1u << ReacIndex(reac))))
        continue;
      Double_t tot[18];
      e.Totals(tot);
      StudyEvt se;
      se.cls = -1;
      se.src = reac;
      se.tagged = reac;
      se.cs[0] = 0.0;
      for (Int_t s = 1; s <= 17; s++)
        se.cs[s] = se.cs[s - 1] + tot[s];
      dataEvts.push_back(se);
    }
    out << Form("(%s) at strip %d: n(ap)=%zu n(an)=%zu n(aa)=%zu n(beam)=%zu "
                "n(data)=%zu",
                "a,p", reac, pop[3].size(), pop[2].size(), pop[1].size(),
                pop[0].size(), dataEvts.size())
        << std::endl;
    if (pop[3].size() < 30 || pop[2].size() < 30) {
      out << "  too few (a,p) or (a,n) tagged here to score any plane."
          << std::endl
          << std::endl;
      continue;
    }

    // The candidate grid: xh and yh sweep every window, each with and
    // without the upstream-mean ratio. The current analysis plane is marked.
    std::vector<StudyPlane> planes;
    for (Int_t xh = reac; xh <= 17; xh++)
      for (Int_t yh = reac + 1; yh <= 17; yh++) {
        StudyPlane p;
        p.xl = 1;
        p.xh = xh;
        p.yl = reac + 1;
        p.yh = yh;
        p.ratio = kTRUE;
        planes.push_back(p);
        p.ratio = kFALSE;
        planes.push_back(p);
      }

    std::vector<StudyRow> rows;
    for (size_t q = 0; q < planes.size(); q++) {
      const StudyPlane &p = planes[q];
      std::vector<Double_t> px[4], py[4];
      for (Int_t k = 0; k < 4; k++) {
        px[k].reserve(pop[k].size());
        py[k].reserve(pop[k].size());
      }
      for (Int_t k = 1; k < 4; k++)
        for (size_t i = 0; i < pop[k].size(); i++) {
          Double_t x = 0.0, y = 0.0;
          StudyXY(*pop[k][i], p, reac, x, y);
          px[k].push_back(x);
          py[k].push_back(y);
        }
      StudyRow row;
      row.p = p;
      row.d_an = StudyLda(px[3], py[3], px[2], py[2]);
      row.d_aa = StudyLda(px[3], py[3], px[1], py[1]);
      row.cap1 =
          StudyCapture(px[3], py[3], px[2], py[2], px[1], py[1], 0.01, 0.02);
      row.cap2 =
          StudyCapture(px[3], py[3], px[2], py[2], px[1], py[1], 0.02, 0.02);
      row.cap5 =
          StudyCapture(px[3], py[3], px[2], py[2], px[1], py[1], 0.05, 0.02);
      row.current = p.xl == C.X_LO && p.xh == C.X_HI && p.yl == YLoOf(reac) &&
                    p.yh == YHiOf(reac) && p.ratio == C.Y_RATIO_TO_UPSTREAM;
      rows.push_back(row);
    }

    // Find the baseline row first (the grid's own current flag).
    const StudyRow *cur = nullptr;
    for (size_t q = 0; q < rows.size(); q++)
      if (rows[q].current)
        cur = &rows[q];
    out << Form("  current plane  x 1->%-2d  y %d->%-2d %-6s d(an)=%6.2f "
                "d(aa)=%6.2f  cap@1%%=%5.1f%%  cap@2%%=%5.1f%%  "
                "cap@5%%=%5.1f%%",
                C.X_HI, YLoOf(reac), YHiOf(reac),
                C.Y_RATIO_TO_UPSTREAM ? "ratio" : "plain",
                cur ? cur->d_an : 0.0, cur ? cur->d_aa : 0.0,
                cur ? 100.0 * cur->cap1 : 0.0, cur ? 100.0 * cur->cap2 : 0.0,
                cur ? 100.0 * cur->cap5 : 0.0)
        << std::endl;

    std::sort(rows.begin(), rows.end(), StudyRowLess);
    out << "  ranked planes (top 12 by cap@2%):" << std::endl;
    const Int_t kTop = 12;
    for (Int_t q = 0; q < Int_t(rows.size()) && q < kTop; q++)
      out << Form("   %2d   x 1->%-2d  y %d->%-2d %-6s d(an)=%6.2f d(aa)=%6.2f "
                  " cap@1%%=%5.1f%%  cap@2%%=%5.1f%%  cap@5%%=%5.1f%%%s",
                  q + 1, rows[q].p.xh, rows[q].p.yl, rows[q].p.yh,
                  rows[q].p.ratio ? "ratio" : "plain", rows[q].d_an,
                  rows[q].d_aa, 100.0 * rows[q].cap1, 100.0 * rows[q].cap2,
                  100.0 * rows[q].cap5, rows[q].current ? "   <- current" : "")
          << std::endl;
    out << std::endl;

    // Figures: the top three and the current plane, simulated populations,
    // the pure-beam and tagged data from the reservoir.
    const char *figTag[3] = {"top1", "top2", "top3"};
    std::vector<StudySample> samples;
    for (Int_t q = 0; q < 3 && q < Int_t(rows.size()); q++) {
      StudyMakeSamples(reac, pop, dataBeamEvts, dataEvts, rows[q].p, samples);
      const TString title = Form(
          "strip %d: x #Sigma s1#rightarrows%d, y #Sigma s%d#rightarrows%d%s",
          reac, rows[q].p.xh, rows[q].p.yl, rows[q].p.yh,
          rows[q].p.ratio ? " / upstream" : "");
      StudyFigure(Form("plane_study_reac%d_%s", reac, figTag[q]), title,
                  samples);
    }
    if (cur) {
      StudyPlane cp;
      cp.xl = C.X_LO;
      cp.xh = C.X_HI;
      cp.yl = YLoOf(reac);
      cp.yh = YHiOf(reac);
      cp.ratio = C.Y_RATIO_TO_UPSTREAM;
      StudyMakeSamples(reac, pop, dataBeamEvts, dataEvts, cp, samples);
      StudyFigure(Form("plane_study_reac%d_current", reac),
                  Form("strip %d: the current analysis plane, selection as "
                       "the fill applies it",
                       reac),
                  samples);
    }
  }

  const std::string text = out.str();
  std::cout << text;
  const TString dir = Paths::ResultsDir() + "/plots/sim_scatter";
  gSystem->mkdir(dir, kTRUE);
  std::ofstream fo((dir + "/plane_study_report.txt").Data());
  fo << text;
  std::cout << "strip-sum-scatter: plane study written to " << dir
            << "/plane_study_report.txt" << std::endl;
}

void StripSumScatter::TemplateFit() {
  const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;
  const Int_t kReacMin = C.REACTION_STRIP_MIN;
  const Int_t kReacMax = C.REACTION_STRIP_MAX;
  std::cout << "strip-sum-scatter: template fit, the data tagged at each "
               "strip as a sum of the simulated populations."
            << std::endl;

  Double_t gain[18];
  if (!SimBeamGains(gain))
    for (Int_t s = 0; s < 18; s++)
      gain[s] = 1.0;
  // The same smearing calibration as the plane study, so the templates
  // carry the measured beam resolution and flux.
  Double_t jit[18];
  Double_t sigmaFlux = 0.0;
  std::ostringstream calib;
  StudySmearingCalibration(gain, C.PLANE_STUDY_EXTRA_JITTER,
                           StudyReservoirBeamTotals(m_reservoir), jit,
                           sigmaFlux, calib);

  std::vector<StudyEvt> evts = StudyCollect(gain, jit, sigmaFlux, 20260930);
  if (evts.empty()) {
    std::cerr << "strip-sum-scatter: no simulated events; skipping the "
                 "template fit."
              << std::endl;
    return;
  }

  // The data yield profile per strip, the migration template's weights:
  // its ap/an/aa mixture follows the physical strip yields, not the sim's
  // uniform generation.
  Int_t nDataAt[10] = {0};
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
    nDataAt[ReacIndex(reac)] = Int_t(StudyDataEvts(reac, m_reservoir).size());

  const Int_t nb = 40;
  std::ostringstream out;
  out << "strip-sum-scatter: template fit (" << Paths::DatasetName().Data()
      << "): the data tagged at each strip, binned " << nb
      << " deep on the y axis of the current analysis plane, fitted as a "
         "sum of the simulated components (each with a free yield) by the "
         "extended binned Poisson likelihood, maximized by EM. The "
         "simulated traces carry the same smearing as the plane study. "
         "The fit is 1-D in y on purpose: the (a,p) proton is charged and "
         "its Q value is 0.7 MeV less negative, so (a,p) sits at higher y "
         "than (a,n), and a 1-D fit is far better conditioned than the 2-D "
         "one, whose clouds are only ~1.7 sigma apart."
      << std::endl;
  if (C.PLANE_STUDY_EXTRA_JITTER < 0.0)
    out << calib.str();
  out << "Components: (a,p), (a,n), (a,a') GENERATED at and TAGGED at the "
         "strip; the migration (born at another strip, tagged here) with "
         "the data yield profile as its composition weight; the beam. The "
         "migration's own composition is the sim's, so its yield is a "
         "shape, not a per-class count."
      << std::endl;
  out << std::endl;

  // Accumulate the per-strip fit rows for template_fit.root.
  struct FitRow {
    Int_t reac;
    Double_t n_data;
    Double_t n_ap, n_an, n_aa, n_mig, n_beam;
    Double_t err_ap, err_an, err_aa, err_mig, err_beam;
    Double_t f_ap, err_fap;
    Double_t chi2, dof;
  };
  std::vector<FitRow> rows;

  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    std::vector<StudyEvt> data = StudyDataEvts(reac, m_reservoir);
    out << Form("strip %d: n_data=%zu", reac, data.size()) << std::endl;
    if (data.size() < 30) {
      out << "  too few data to fit." << std::endl << std::endl;
      continue;
    }

    // The current analysis plane for this strip (only its y axis is used).
    StudyPlane p;
    p.xl = C.X_LO;
    p.xh = C.X_HI;
    p.yl = YLoOf(reac);
    p.yh = YHiOf(reac);
    p.ratio = C.Y_RATIO_TO_UPSTREAM;

    // Component y lists and weights. Index 0 ap, 1 an, 2 aa, 3 migration,
    // 4 beam.
    std::vector<Double_t> cy[5], cw[5];
    for (size_t i = 0; i < evts.size(); i++) {
      const StudyEvt &e = evts[i];
      if (e.tagged != reac)
        continue;
      Double_t x = 0.0, y = 0.0;
      StudyXY(e, p, reac, x, y);
      if (e.cls > 0 && e.src == reac) {
        // cls 3 ap -> 0, 2 an -> 1, 1 aa -> 2, the component order.
        const Int_t idx = (e.cls == 3) ? 0 : (e.cls == 2) ? 1 : 2;
        cy[idx].push_back(y);
        cw[idx].push_back(1.0);
      } else if (e.cls > 0 && e.src != reac) {
        const Int_t sri = ReacIndex(e.src);
        cy[3].push_back(y);
        cw[3].push_back(nDataAt[sri] > 0 ? Double_t(nDataAt[sri]) : 0.0);
      } else if (e.cls == 0) {
        cy[4].push_back(y);
        cw[4].push_back(1.0);
      }
    }
    // The data y values.
    std::vector<Double_t> dy;
    for (size_t i = 0; i < data.size(); i++) {
      Double_t x = 0.0, y = 0.0;
      StudyXY(data[i], p, reac, x, y);
      dy.push_back(y);
    }

    // Drop empty components, keeping the index for the report.
    Int_t compIdx[5] = {0, 1, 2, 3, 4};
    Int_t nComp = 0;
    for (Int_t a = 0; a < 5; a++) {
      const Int_t k = compIdx[a];
      Double_t wsum = 0.0;
      for (size_t i = 0; i < cw[k].size(); i++)
        wsum += cw[k][i];
      if (cw[k].size() >= 20 && wsum > 0.0) {
        compIdx[nComp++] = k;
      }
    }
    if (nComp < 2) {
      out << "  too few simulated components to fit." << std::endl << std::endl;
      continue;
    }

    // Bin range from the 2nd-98th percentile of the data y.
    const Double_t ylo = StudyQuantile(dy, 0.02), yhi = StudyQuantile(dy, 0.98);
    if (!(yhi > ylo)) {
      out << "  degenerate data range; skipping." << std::endl << std::endl;
      continue;
    }
    const Double_t ybin = (yhi - ylo) / Double_t(nb);

    // Data binned counts.
    std::vector<Double_t> d(nb, 0.0);
    for (size_t i = 0; i < dy.size(); i++) {
      const Int_t iy = Int_t((dy[i] - ylo) / ybin);
      if (iy < 0 || iy >= nb)
        continue;
      d[iy] += 1.0;
    }

    // Templates, unit-summed, weighted.
    std::vector<std::vector<Double_t>> P;
    std::vector<TString> labels;
    std::vector<Int_t> colors;
    for (Int_t a = 0; a < nComp; a++) {
      const Int_t k = compIdx[a];
      std::vector<Double_t> pk(nb, 0.0);
      Double_t wsum = 0.0;
      for (size_t i = 0; i < cy[k].size(); i++) {
        const Int_t iy = Int_t((cy[k][i] - ylo) / ybin);
        if (iy < 0 || iy >= nb)
          continue;
        pk[iy] += cw[k][i];
        wsum += cw[k][i];
      }
      if (wsum <= 0.0)
        continue;
      for (Int_t b = 0; b < nb; b++)
        pk[b] /= wsum;
      P.push_back(pk);
      const Int_t idx = k;
      if (idx == 0) {
        labels.push_back("(#alpha,p)");
        colors.push_back(kGreen + 2);
      } else if (idx == 1) {
        labels.push_back("(#alpha,n)");
        colors.push_back(kRed + 1);
      } else if (idx == 2) {
        labels.push_back("(#alpha,#alpha')");
        colors.push_back(kAzure + 2);
      } else if (idx == 3) {
        labels.push_back("migration");
        colors.push_back(kOrange + 2);
      } else {
        labels.push_back("beam");
        colors.push_back(kGray + 2);
      }
    }
    if (P.size() < 2) {
      out << "  too few in-bin components to fit." << std::endl << std::endl;
      continue;
    }

    // DIAGNOSTIC: the y moments of the data and each template, so the
    // templates can be compared to the data directly. The beam (data and
    // sim) is included to separate a scale problem from a reaction-shape
    // problem: the beam is flat, so its y is ~1 with no reaction content.
    {
      const char *cname[5] = {"(a,p)", "(a,n)", "(a,a')", "migr", "beam"};
      Double_t m = 0.0, sd = 0.0;
      StudyMoments(dy, m, sd);
      out << Form("  y moments (mean +- sd): data n=%zu  %.4f +- %.4f",
                  dy.size(), m, sd)
          << std::endl;
      std::vector<Double_t> dataBeamY;
      for (size_t k = 0; k < m_reservoir.size(); k++) {
        const TraceEvt &e = m_reservoir[k];
        if (!e.beam_flat)
          continue;
        Double_t tot[18];
        e.Totals(tot);
        StudyEvt se;
        se.cls = -1;
        se.src = reac;
        se.tagged = -1;
        se.cs[0] = 0.0;
        for (Int_t s = 1; s <= 17; s++)
          se.cs[s] = se.cs[s - 1] + tot[s];
        Double_t bx = 0.0, by = 0.0;
        StudyXY(se, p, reac, bx, by);
        dataBeamY.push_back(by);
      }
      Double_t bm = 0.0, bsd = 0.0;
      StudyMoments(dataBeamY, bm, bsd);
      out << Form("   %-7s n=%5zu  %.4f +- %.4f  (data beam)", "beam",
                  dataBeamY.size(), bm, bsd)
          << std::endl;
      std::vector<Double_t> simBeamY;
      for (size_t i = 0; i < evts.size(); i++) {
        if (evts[i].cls != 0)
          continue;
        Double_t bx = 0.0, by = 0.0;
        StudyXY(evts[i], p, reac, bx, by);
        simBeamY.push_back(by);
      }
      Double_t sbm = 0.0, sbd = 0.0;
      StudyMoments(simBeamY, sbm, sbd);
      out << Form("   %-7s n=%5zu  %.4f +- %.4f  (sim beam)", "beam",
                  simBeamY.size(), sbm, sbd)
          << std::endl;
      for (Int_t a = 0; a < nComp; a++) {
        const Int_t k = compIdx[a];
        Double_t cm = 0.0, csd = 0.0;
        StudyMoments(cy[k], cm, csd);
        out << Form("   %-7s n=%5zu  %.4f +- %.4f", cname[k], cy[k].size(), cm,
                    csd)
            << std::endl;
      }
      // The binned data y, so the shape (a single shifted peak vs a reaction
      // peak plus a low-y beam tail) can be read directly.
      out << Form("  data y bins [%.2f..%.2f]:", ylo, yhi);
      for (Int_t b = 0; b < nb; b++) {
        if (b % 8 == 7)
          out << std::endl;
        out << Form(" %4.0f", d[b]);
      }
      out << std::endl;
    }

    const Int_t K = Int_t(P.size());
    const StudyFitResult r = StudyFitEM(nb, d, P);
    Double_t ntot = 0.0;
    for (Int_t k = 0; k < K; k++)
      ntot += r.N[k];

    // The fitted yields back-mapped to the component order. P keeps the
    // components in compIdx order minus any dropped for being empty
    // in-bin, so rebuild that map the same way the templates were.
    Double_t N[5] = {0.0}, E[5] = {0.0};
    std::vector<Int_t> pToComp;
    for (Int_t a = 0; a < nComp; a++) {
      const Int_t k = compIdx[a];
      Double_t wsum = 0.0;
      for (size_t i = 0; i < cw[k].size(); i++) {
        const Int_t iy = Int_t((cy[k][i] - ylo) / ybin);
        if (iy >= 0 && iy < nb)
          wsum += cw[k][i];
      }
      if (wsum > 0.0)
        pToComp.push_back(k);
    }
    for (Int_t b = 0; b < K && b < Int_t(pToComp.size()); b++) {
      N[pToComp[b]] = r.N[b];
      E[pToComp[b]] =
          r.cov[b * K + b] > 0.0 ? std::sqrt(r.cov[b * K + b]) : 0.0;
    }

    // f(a,p) and its error, by covariance propagation on N_ap / total.
    Double_t f_ap = ntot > 0.0 ? N[0] / ntot : 0.0;
    Double_t err_fap = 0.0;
    Int_t apos = -1;
    for (Int_t b = 0; b < Int_t(pToComp.size()); b++)
      if (pToComp[b] == 0)
        apos = b;
    if (apos >= 0 && ntot > 0.0) {
      // f = N_a / sum N; df/dN_a = 1 - f, df/dN_j = -f for j != a.
      Double_t v = 0.0;
      for (Int_t i = 0; i < K; i++) {
        const Double_t gi = (i == apos) ? (1.0 - f_ap) : -f_ap;
        for (Int_t j = 0; j < K; j++) {
          const Double_t gj = (j == apos) ? (1.0 - f_ap) : -f_ap;
          v += gi * gj * r.cov[i * K + j];
        }
      }
      err_fap = v > 0.0 ? std::sqrt(v) : 0.0;
    }

    const Double_t chi2dof = (nb - K) > 0 ? r.dev / Double_t(nb - K) : 0.0;
    out << Form("  fitted (EM converged=%s, %d iters):",
                r.converged ? "yes" : "no", r.iters)
        << std::endl;
    const char *cname[5] = {"(a,p)", "(a,n)", "(a,a')", "migr", "beam"};
    for (Int_t b = 0; b < Int_t(pToComp.size()); b++) {
      const Int_t k = pToComp[b];
      out << Form("    %-7s  %8.1f +- %6.1f", cname[k], N[k], E[k])
          << std::endl;
    }
    out << Form("  total  %8.1f  (data %zu)", ntot, data.size()) << std::endl;
    out << Form("  chi2/dof  %6.2f  (%d bins, %d components)%s", chi2dof, nb, K,
                r.singular ? "  [covariance singular]" : "")
        << std::endl;
    out << Form("  f(a,p) = %6.3f +- %5.3f   (fraction of tagged that is "
                "(a,p))",
                f_ap, err_fap)
        << std::endl;
    out << std::endl;

    FitRow row;
    row.reac = reac;
    row.n_data = Double_t(data.size());
    row.n_ap = N[0];
    row.n_an = N[1];
    row.n_aa = N[2];
    row.n_mig = N[3];
    row.n_beam = N[4];
    row.err_ap = E[0];
    row.err_an = E[1];
    row.err_aa = E[2];
    row.err_mig = E[3];
    row.err_beam = E[4];
    row.f_ap = f_ap;
    row.err_fap = err_fap;
    row.chi2 = r.dev;
    row.dof = Double_t(nb - K);
    rows.push_back(row);

    // The figure.
    std::vector<Double_t> mu(nb, 0.0);
    for (Int_t b = 0; b < nb; b++)
      for (Int_t k = 0; k < K; k++)
        mu[b] += r.N[k] * P[k][b];
    std::vector<Double_t> ye(nb + 1);
    for (Int_t i = 0; i <= nb; i++)
      ye[i] = ylo + i * ybin;
    std::vector<Double_t> Evec(K, 0.0);
    for (Int_t k = 0; k < K; k++)
      Evec[k] = r.cov[k * K + k] > 0.0 ? std::sqrt(r.cov[k * K + k]) : 0.0;
    StudyFitFigure(
        Form("template_fit_reac%d", reac),
        Form("strip %d: 1-D template fit, y #Sigma s%d#rightarrows%d%s", reac,
             p.yl, p.yh, p.ratio ? " / upstream" : ""),
        nb, ye.data(), d, mu, labels, colors, r.N, Evec, P);
  }

  const std::string text = out.str();
  std::cout << text;
  const TString dir = Paths::ResultsDir() + "/plots/sim_scatter";
  gSystem->mkdir(dir, kTRUE);
  std::ofstream fo((dir + "/template_fit_report.txt").Data());
  fo << text;

  // The per-strip fit rows, for a cross-section to read the (a,p) scale.
  TFile *rf = new TFile((dir + "/template_fit.root").Data(), "RECREATE");
  TTree *tr = new TTree("template_fit", "Template fit per strip");
  Int_t reac_i;
  Double_t n_data, n_ap, n_an, n_aa, n_mig, n_beam;
  Double_t err_ap, err_an, err_aa, err_mig, err_beam;
  Double_t f_ap, err_fap, chi2, dof;
  tr->Branch("reac", &reac_i, "reac/I");
  tr->Branch("n_data", &n_data, "n_data/D");
  tr->Branch("n_ap", &n_ap, "n_ap/D");
  tr->Branch("n_an", &n_an, "n_an/D");
  tr->Branch("n_aa", &n_aa, "n_aa/D");
  tr->Branch("n_mig", &n_mig, "n_mig/D");
  tr->Branch("n_beam", &n_beam, "n_beam/D");
  tr->Branch("err_ap", &err_ap, "err_ap/D");
  tr->Branch("err_an", &err_an, "err_an/D");
  tr->Branch("err_aa", &err_aa, "err_aa/D");
  tr->Branch("err_mig", &err_mig, "err_mig/D");
  tr->Branch("err_beam", &err_beam, "err_beam/D");
  tr->Branch("f_ap", &f_ap, "f_ap/D");
  tr->Branch("err_fap", &err_fap, "err_fap/D");
  tr->Branch("chi2", &chi2, "chi2/D");
  tr->Branch("dof", &dof, "dof/D");
  for (size_t i = 0; i < rows.size(); i++) {
    reac_i = rows[i].reac;
    n_data = rows[i].n_data;
    n_ap = rows[i].n_ap;
    n_an = rows[i].n_an;
    n_aa = rows[i].n_aa;
    n_mig = rows[i].n_mig;
    n_beam = rows[i].n_beam;
    err_ap = rows[i].err_ap;
    err_an = rows[i].err_an;
    err_aa = rows[i].err_aa;
    err_mig = rows[i].err_mig;
    err_beam = rows[i].err_beam;
    f_ap = rows[i].f_ap;
    err_fap = rows[i].err_fap;
    chi2 = rows[i].chi2;
    dof = rows[i].dof;
    tr->Fill();
  }
  tr->Write();
  delete tr;
  delete rf;

  std::cout << "strip-sum-scatter: template fit written to " << dir
            << "/template_fit_report.txt and template_fit.root" << std::endl;
}

/// Per reaction strip, overlay TRACES_PER_CLASS sampled per-strip traces of
/// every sim population the control directory holds, in the experimental
/// DrawRegionTraces style (beam grey, (a,a') azure, (a,n) red, (a,p) green,
/// SimClassColor for the rest). The beam reference is the same for every
/// strip. Sampled fresh each run (40 traces/file is trivial).
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
  // One file per class and strip (BuildFileSpecs stops on a second).
  // Classes in the order they are first seen (the specs are sorted by tag).
  std::vector<TString> bases;
  std::map<TString, std::map<Int_t, TString>> files; // base, strip -> file
  std::vector<TString> beam_files;
  for (Int_t i = 0; i < Int_t(specs.size()); i++) {
    const TString base = RemixSim::TagWithoutStrip(specs[i].tag);
    const TString file = RemixSim::SimRootPath(specs[i]);
    const Int_t strip = RemixSim::ReactionStripOf(specs[i].tag);
    if (base == "beam") {
      beam_files.push_back(file);
    } else if (strip >= kReacMin && strip <= kReacMax) {
      if (!files.count(base))
        bases.push_back(base);
      files[base][strip] = file;
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
    std::vector<TraceClass> classes;
    Bool_t any = kFALSE;
    for (size_t b = 0; b < bases.size(); b++) {
      TraceClass cl;
      cl.label = PrettyLabel(bases[b]);
      cl.color = SimClassColor(bases[b]);
      if (files[bases[b]].count(r))
        cl.traces = SimPopTraces(files[bases[b]][r], gain, kTracesPerRegion);
      any = any || !cl.traces.empty();
      classes.push_back(cl);
    }
    if (!any)
      continue;
    DrawRegionTraces(Form("sim_region_traces_reac%d", r), "sim_scatter",
                     beam_traces, classes, 0.6, 1.6, "#DeltaE [a.u.]");
    for (size_t b = 0; b < classes.size(); b++)
      for (Int_t i = 0; i < Int_t(classes[b].traces.size()); i++)
        delete classes[b].traces[i];
  }
  for (Int_t i = 0; i < Int_t(beam_traces.size()); i++)
    delete beam_traces[i];
}

/// Fingerprint of the sim inputs + window geometry: each sim file's size+mtime
/// (cheap, no open) plus the reaction-strip range and x window. Regenerating
/// the sim (new mtimes) or changing the windows invalidates the cached overlay.
TString StripSumScatter::SimFingerprint(
    const std::vector<RemixSim::SimFileSpec> &specs) {
  // v2: per-strip norm via SimBeamGains(); the gain file is stamped by the
  // per-spec loop; v3: norm hardcoded to 1 a.u. (NORM_MUSIC_MEV gone).
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kXLo = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_LO;
  const Int_t kXHi = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.X_HI;

  // v4: one population per class and strip, so caches that held every
  // population twice are rebuilt.
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
    // One population per class and strip (BuildFileSpecs stops on a second).
    std::map<std::pair<TString, Int_t>, std::pair<SimPop, Bool_t>> chosen;
    for (Int_t i = 0; i < Int_t(specs.size()); i++) {
      const TString base = RemixSim::TagWithoutStrip(specs[i].tag);
      const Int_t strip = RemixSim::ReactionStripOf(specs[i].tag);
      std::pair<TString, Int_t> key(base, strip);
      SimPop p;
      p.file = RemixSim::SimRootPath(specs[i]);
      p.label = PrettyLabel(specs[i].tag);
      chosen[key] = std::make_pair(p, kTRUE);
    }
    std::map<Int_t, std::vector<SimPop>> reacted;
    std::vector<SimPop> refs;
    for (std::map<std::pair<TString, Int_t>,
                  std::pair<SimPop, Bool_t>>::const_iterator it =
             chosen.begin();
         it != chosen.end(); ++it) {
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
      // grey, (a,a') azure, (a,n) red, (a,p) green -- keyed off the
      // population label, same values as SimClassColor.
      TString lab = g->GetTitle();
      Int_t color = kBlack;
      if (lab == "Beam")
        color = kGray + 2;
      else if (lab == "(#alpha,#alpha')")
        color = kAzure + 2;
      else if (lab == "(#alpha,n)")
        color = kRed + 1;
      else if (lab == "(#alpha,p)")
        color = kGreen + 2;
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
  m_normedVar.clear();
  if (TNamed *vn = dynamic_cast<TNamed *>(cf->Get("cut_variants"))) {
    TString names = vn->GetTitle();
    Int_t from = 0;
    TString tok;
    while (names.Tokenize(tok, from, ",")) {
      if (tok.IsNull())
        continue;
      m_variantNames.push_back(tok);
      std::vector<Long64_t> counts(m_tagged.size(), 0),
          normed(m_tagged.size(), 0);
      for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
        if (TParameter<Long64_t> *p = dynamic_cast<TParameter<Long64_t> *>(
                cf->Get(Form("n_tagged_r%d_%s", reac, tok.Data()))))
          counts[ReacIndex(reac)] = p->GetVal();
        if (TParameter<Long64_t> *p = dynamic_cast<TParameter<Long64_t> *>(
                cf->Get(Form("n_normed_r%d_%s", reac, tok.Data()))))
          normed[ReacIndex(reac)] = p->GetVal();
      }
      m_taggedVar.push_back(counts);
      m_normedVar.push_back(normed);
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
  // The cut-variation counts: the variant names as one list, then per
  // variant and strip the tagged count and the denominator.
  {
    TString names;
    for (size_t v = 0; v < m_variantNames.size(); v++)
      names += (v ? "," : "") + m_variantNames[v];
    TNamed("cut_variants", names.Data()).Write();
    for (size_t v = 0; v < m_variantNames.size(); v++)
      for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
        TParameter<Long64_t>(
            Form("n_tagged_r%d_%s", reac, m_variantNames[v].Data()),
            m_taggedVar[v][ReacIndex(reac)])
            .Write();
        TParameter<Long64_t>(
            Form("n_normed_r%d_%s", reac, m_variantNames[v].Data()),
            m_normedVar[v][ReacIndex(reac)])
            .Write();
      }
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
SingleRunFitResult StripSumScatter::FitRunGates(Int_t key, const TString &label,
                                                TChain *chain) {
  (void)key;
  const TString run = label; // the group's name in every message
  SingleRunFitResult res;
  if (!chain || chain->GetEntries() == 0)
    return res;

  // --- Beam classification ellipses (no gating) ---
  BeamEllipses be;
  be.ok = kFALSE;
  {
    const TString tag = label;
    const TString subdir = "strip_sum_scatter/" + label;
    Int_t ent_sx = 0, ent_sy = 1;
    EntranceAxes(ent_sx, ent_sy);
    const TString ent_tag = AxisTag(ent_sx) + "/" + AxisTag(ent_sy);
    const Double_t kPure =
        Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PURE_BEAM_NSIGMA;
    BeamFit2D ent_ell = FindBeamGate(chain, ent_sx, ent_sy, tag, subdir, kPure);
    if (ent_sy == GATE_AXIS_GRID)
      be.s1_grid = ent_ell;
    else if (ent_sx == 0)
      be.s0_s1 = ent_ell;
    else
      be.s1_s2 = ent_ell;
    if (ent_ell.ok) {
      std::lock_guard<std::mutex> lk(g_log_mutex);
      Constants::Detail() << "  " << run << " beam ellipse " << ent_tag
                          << ": mu=(" << ent_ell.mu_x << "," << ent_ell.mu_y
                          << ")" << std::endl;
    } else {
      std::lock_guard<std::mutex> lk(g_log_mutex);
      std::cerr << "  " << run << " beam ellipse " << ent_tag
                << " failed; skipping run" << std::endl;
      return res;
    }
    if (Constants::cfg.IGNORE_STRIP_17) {
      be.use_s15_s16 = kTRUE;
      be.s15_s16 = FindBeamGate(chain, 15, 16, tag, subdir, kPure);
      if (be.s15_s16.ok) {
        std::lock_guard<std::mutex> lk(g_log_mutex);
        Constants::Detail()
            << "  " << run << " beam ellipse s15/s16: mu=(" << be.s15_s16.mu_x
            << "," << be.s15_s16.mu_y << ")" << std::endl;
      } else {
        std::lock_guard<std::mutex> lk(g_log_mutex);
        std::cerr << "  " << run << " beam ellipse s15/s16 failed; skipping run"
                  << std::endl;
        return res;
      }
    } else {
      be.use_s15_s16 = kFALSE;
      be.s16_s17 = FindBeamGate(chain, 16, 17, tag, subdir, kPure);
      if (be.s16_s17.ok) {
        std::lock_guard<std::mutex> lk(g_log_mutex);
        Constants::Detail()
            << "  " << run << " beam ellipse s16/s17: mu=(" << be.s16_s17.mu_x
            << "," << be.s16_s17.mu_y << ")" << std::endl;
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

  // --- The strip gate ---
  const StripSumScatterConfig &C = Constants::cfg.STRIP_SUM_SCATTER_CONFIG;
  res.gate =
      FindStripGate(chain, C.GATE_STRIP, label, "strip_sum_scatter/" + label,
                    C.GATE_NSIGMA, C.GATE_CENTER);
  if (!res.gate.ok) {
    std::lock_guard<std::mutex> lk(g_log_mutex);
    std::cerr << "  " << run << " strip gate s" << C.GATE_STRIP
              << " failed; skipping run" << std::endl;
  }
  res.ok = res.gate.ok;
  return res;
}

// Which variants move the beam selection: those re-evaluate the event
// level themselves; the rest share the nominal decision.
static void VariantMovesEventLevel(
    const std::vector<std::pair<TString, TagThresholds>> &variants,
    const TagThresholds &nominal, std::vector<Bool_t> &varies_event) {
  varies_event.assign(variants.size(), kFALSE);
  for (size_t v = 0; v < variants.size(); v++) {
    const TagThresholds &T = variants[v].second;
    varies_event[v] = T.gate_nsigma != nominal.gate_nsigma ||
                      T.pileup_nsigma != nominal.pileup_nsigma ||
                      T.noise_nsigma != nominal.noise_nsigma;
  }
}

// The trace record: calibrated and raw ADC strips, the strip-0/-17 ends,
// the mirrored short side, the both-channel multiplicity, seed, mask, beam.
static void FillTraceEvt(TraceEvt &e, const EnergyView &ev, ULong64_t seed_ts,
                         UInt_t mask, Bool_t beam) {
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
  if (Constants::ActiveIgnoreShortStrips())
    for (Int_t s = 1; s <= 16; s++) {
      if ((s % 2) != 0)
        e.rightdE_adc[s - 1] = 0.0f;
      else
        e.leftdE_adc[s - 1] = 0.0f;
    }
  e.both_mult = CountBothEnds(ev, 16);
  e.seed_ts = seed_ts;
  e.reac_mask = mask;
  e.beam_flat = beam;
}

StripSumScatter::EventLevelVerdict
StripSumScatter::JudgeEventLevel(const BeamGate1D &gate, const EnergyView &ev,
                                 Int_t gate_strip, Double_t step_z,
                                 const TagThresholds &T) {
  EventLevelVerdict v;
  v.gate = PassesGate(gate, ev, gate_strip, T.gate_nsigma);
  v.pileup = IsPileup(ev, T.pileup_nsigma);
  v.noise = IsNoise(ev, T.noise_nsigma);
  v.smooth_ok = !(T.smooth_nsigma > 0.0) || step_z <= T.smooth_nsigma;
  v.event = v.gate && !v.pileup && !v.noise;
  // Sequential rejection: the first cut that fails counts.
  if (!v.gate)
    v.pre_cut = kPreGate;
  else if (v.pileup)
    v.pre_cut = kPrePileup;
  else if (v.noise)
    v.pre_cut = kPreNoise;
  else if (!v.smooth_ok)
    v.pre_cut = kPreSmooth;
  return v;
}

Bool_t StripSumScatter::AnyVariantKeeps(
    const BeamGate1D &gate, const EnergyView &ev, Int_t gate_strip,
    const std::vector<std::pair<TString, TagThresholds>> &variants,
    const std::vector<Bool_t> &varies_event, Bool_t nom_event,
    std::vector<Bool_t> &var_event) {
  // Every variant's event-level decision, before the smoothness it tests
  // itself.
  var_event.assign(variants.size(), nom_event);
  for (size_t v = 0; v < variants.size(); v++)
    if (varies_event[v]) {
      const TagThresholds &T = variants[v].second;
      var_event[v] = PassesGate(gate, ev, gate_strip, T.gate_nsigma) &&
                     !IsPileup(ev, T.pileup_nsigma) &&
                     !IsNoise(ev, T.noise_nsigma);
    }
  Bool_t any_event = nom_event;
  for (size_t v = 0; !any_event && v < variants.size(); v++)
    any_event = var_event[v];
  return any_event;
}

void StripSumScatter::TagUnderVariants(
    const EnergyView &ev, Double_t step_z, Int_t kReacMin, Int_t kReacMax,
    const std::vector<std::pair<TString, TagThresholds>> &variants,
    const std::vector<Bool_t> &var_event, std::vector<Bool_t> &pass,
    SingleRunFillResult &res) {
  for (size_t v = 0; v < variants.size(); v++) {
    const TagThresholds &T = variants[v].second;
    if (!var_event[v] || (T.smooth_nsigma > 0.0 && step_z > T.smooth_nsigma))
      continue;
    for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
      pass[ReacIndex(reac)] = RejectReason(ev, reac, T) == kTagPass;
    const Int_t won = ResolveTag(ev, pass);
    if (won >= 0)
      res.tagged_var[v][ReacIndex(won)]++;
    for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
      if ((won < 0 || reac <= won) && BeamUpstreamOf(ev, reac, T))
        res.normed_var[v][ReacIndex(reac)]++;
  }
}

Int_t StripSumScatter::TagNominal(const EnergyView &ev,
                                  const Double_t totals[18], Int_t kReacMin,
                                  Int_t kReacMax, Int_t nReacStrips,
                                  std::vector<Bool_t> &pass,
                                  SingleRunFillResult &res, UInt_t &mask) {
  std::vector<TagCut> why(nReacStrips, kTagPass);
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    why[ReacIndex(reac)] = RejectReason(ev, reac);
    pass[ReacIndex(reac)] = why[ReacIndex(reac)] == kTagPass;
  }
  const Int_t won = ResolveTag(ev, pass);
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    TagCut w = why[ReacIndex(reac)];
    if (w == kTagPass && reac != won)
      w = kCutOtherTag;
    res.cut_counts[ReacIndex(reac) * kNTagCuts + w]++;
  }
  if (won >= 0) {
    mask |= (1u << ReacIndex(won));
    res.tagged[ReacIndex(won)]++;
    Double_t x = 0.0, y = 0.0;
    PlaneXY(totals, won, x, y);
    res.scatters[ReacIndex(won)]->Fill(x, y);
  }
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++)
    if ((won < 0 || reac <= won) && BeamUpstreamOf(ev, reac))
      res.normed_at[ReacIndex(reac)]++;
  return won;
}

void StripSumScatter::AllocateRunScatters(Int_t run, Int_t kReacMin,
                                          Int_t kReacMax, Int_t kXBins,
                                          Int_t kYBins,
                                          std::vector<TH2F *> &scatters) {
  // Private, directory-less clones over the same fixed build range as the
  // merged ones.
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    TH2F *h =
        new TH2F(Form("scatter_r%d_run%d", reac, run), "", kXBins,
                 ScatterBuildRange::kXMin, ScatterBuildRange::kXMax, kYBins,
                 ScatterBuildRange::kYMin, ScatterBuildRange::kYMax);
    h->SetDirectory(nullptr);
    scatters[ReacIndex(reac)] = h;
  }
}

/// One run's scatter fill. Each call builds PRIVATE scatter histograms and its
/// own reservoir slice, so the runs never touch shared state; the caller merges
/// them in run order, which makes the threaded result identical to sequential.
SingleRunFillResult
StripSumScatter::FillRunScatters(Int_t key, const TString &label, TChain *chain,
                                 const BeamGate1D &gate,
                                 const BeamEllipses &runBeam) {
  const Int_t run = key; // unique per task: names the private histograms
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kXBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.XBINS;
  const Int_t kYBins = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.YBINS;
  const Int_t kGateStrip = Constants::cfg.STRIP_SUM_SCATTER_CONFIG.GATE_STRIP;
  const Int_t kBeamReservoirCap =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRACES_PER_CLASS * 10;
  const Int_t nReacStrips = kReacMax - kReacMin + 1;

  SingleRunFillResult res;
  res.scatters.assign(nReacStrips, nullptr);
  if (!chain)
    return res;
  AllocateRunScatters(run, kReacMin, kReacMax, kXBins, kYBins, res.scatters);
  Long64_t totalGated = 0, totalSeen = 0, totalNormed = 0;
  res.tagged.assign(nReacStrips, 0);
  res.normed_at.assign(nReacStrips, 0);
  res.cut_counts.assign(nReacStrips * kNTagCuts, 0);
  res.pre_counts.assign(kNPreCuts, 0);
  // The cut-variation set, fixed for the run: each event past the event-level
  // cuts is tagged once per variant as well.
  const std::vector<std::pair<TString, TagThresholds>> variants =
      ThresholdVariants();
  const TagThresholds nominal = NominalThresholds();
  res.tagged_var.assign(variants.size(), std::vector<Long64_t>(nReacStrips, 0));
  res.normed_var.assign(variants.size(), std::vector<Long64_t>(nReacStrips, 0));
  std::vector<Bool_t> variesEvent;
  VariantMovesEventLevel(variants, nominal, variesEvent);
  Int_t nBeamKept = 0;
  EnergyView ev;
  ev.Attach(chain);
  EnableEventBranches(chain);
  ULong64_t seed_ts_in = 0;
  if (chain->GetBranch("SeedTs"))
    chain->SetBranchAddress("SeedTs", &seed_ts_in);
  Long64_t n = chain->GetEntries();
  {
    std::lock_guard<std::mutex> lk(g_log_mutex);
    Constants::Detail() << "  " << label << ": filling " << nReacStrips
                        << " reaction-strip scatters over " << n << " events..."
                        << std::endl;
  }

  for (Long64_t j = 0; j < n; j++) {
    chain->GetEntry(j);
    ev.Decode();
    totalSeen++;

    // An event with a strip that did not fire is incomplete whatever strip
    // is asked about, so it goes before anything per strip.
    if (!AllStripsFired(ev)) {
      res.pre_counts[kPreAllStrips]++;
      continue;
    }
    // The beam selection, decided at the nominal level for the count and
    // at each variant's level for the cut variation; sequential rejection.
    const Double_t step_z = MaxStepNSigma(ev);
    const EventLevelVerdict nom =
        JudgeEventLevel(gate, ev, kGateStrip, step_z, nominal);
    if (nom.pre_cut >= 0)
      res.pre_counts[nom.pre_cut]++;
    std::vector<Bool_t> var_event;
    if (!AnyVariantKeeps(gate, ev, kGateStrip, variants, variesEvent, nom.event,
                         var_event))
      continue;
    // Both-ends multiplicity: counted on raw ADC so it sees the short end
    // even when IGNORE_SHORT_STRIPS zeroes it in the decode.
    if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.BOTH_MULT_MAX >= 0) {
      const Int_t hi = TMath::Min(
          16, Constants::cfg.STRIP_SUM_SCATTER_CONFIG.BOTH_MULT_COUNT_TO);
      if (CountBothEnds(ev, hi) >
          Constants::cfg.STRIP_SUM_SCATTER_CONFIG.BOTH_MULT_MAX) {
        // Sequential: counted here only if the rest let it through. Not
        // varied: it drops the event for every variant too.
        if (nom.event && nom.smooth_ok)
          res.pre_counts[kPreBothMult]++;
        continue;
      }
    }
    // One tag per event: every strip's verdict, then the attributed strip
    // (ResolveTag); a passing strip that lost it counts as "other strip".

    // Per-strip denominator: what met a reaction's conditions at `reac`
    // and had not reacted before; N(reac+1) = N(reac) - reactions at reac.
    UInt_t mask = 0;
    Double_t totals[18];
    ev.Totals(totals);
    std::vector<Bool_t> pass(nReacStrips, kFALSE);
    TagUnderVariants(ev, step_z, kReacMin, kReacMax, variants, var_event, pass,
                     res);
    if (nom.event && nom.smooth_ok) {
      res.pre_counts[kPrePass]++;
      // The last pre-reaction count: here, not at `seen`, is what makes the
      // tag-count ratio a cross section: same gate + quality efficiencies.
      totalNormed++;
      TagNominal(ev, totals, kReacMin, kReacMax, nReacStrips, pass, res, mask);
    }

    // Keep reaction-passing events for traces; cap pure-beam (mutually
    // exclusive, no reaction jump); a per-task bound: the merge re-caps it.
    if (!nom.event || !nom.smooth_ok)
      continue;
    Bool_t beam = (mask == 0) && IsPureBeam(ev, runBeam);
    if (mask == 0 && !(beam && nBeamKept < kBeamReservoirCap))
      continue;
    if (beam)
      nBeamKept++;

    TraceEvt e;
    FillTraceEvt(e, ev, seed_ts_in, mask, beam);
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

// Phase 1: one run's beam ellipses and strip gate. The first groups of the
// epoch draw their gate figures; the rest are silent.
void StripSumScatter::FitGateRun(Int_t i, std::vector<SingleRunFitResult> &fits,
                                 const std::vector<Int_t> &runOrder,
                                 const std::vector<TString> &labelVec,
                                 const std::vector<TChain *> &chainVec) {
  Constants::SetPlotsThisFile(Constants::InPlotSample(i));
  fits[i] = FitRunGates(runOrder[i], labelVec[i], chainVec[i]);
}

// Merge one file's fill into the totals; the pure-beam events are capped
// GLOBALLY, not per task, so the kept beam does not scale with file count.
void StripSumScatter::MergeRunFill(Int_t t,
                                   std::vector<SingleRunFillResult> &fills,
                                   Int_t &nBeamKept, Long64_t &totalGated,
                                   Long64_t &totalSeen) {
  const Int_t kReacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t kReacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  const Int_t kBeamReservoirCap =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TRACES_PER_CLASS * 10;
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    Int_t ri = ReacIndex(reac);
    if (ri < Int_t(fills[t].scatters.size()) && fills[t].scatters[ri])
      m_scatter[reac]->Add(fills[t].scatters[ri]);
  }
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
      for (Int_t k = 0; k < Int_t(fills[t].tagged_var[v].size()); k++) {
        m_taggedVar[v][k] += fills[t].tagged_var[v][k];
        m_normedVar[v][k] += fills[t].normed_var[v][k];
      }
  for (Int_t k = 0; k < Int_t(fills[t].cut_counts.size()); k++)
    m_cutCounts[k] += fills[t].cut_counts[k];
  for (Int_t k = 0; k < Int_t(fills[t].pre_counts.size()); k++)
    m_preCounts[k] += fills[t].pre_counts[k];
  for (Int_t k = 0; k < Int_t(fills[t].scatters.size()); k++)
    delete fills[t].scatters[k];
  fills[t].scatters.clear();
  std::vector<TraceEvt>().swap(fills[t].reservoir);
}

// Phase 2: one file of the fill; finished prefixes merge in task order,
// freeing on the spot; the budget stays identical to sequential.
void StripSumScatter::FillScatterFile(
    Int_t t, const std::vector<FillTask> &tasks,
    const std::vector<Int_t> &runOrder, const std::vector<TString> &labelVec,
    const std::vector<SingleRunFitResult> &fits, FillPhase &phase) {
  const Int_t i = tasks[t].run_idx;
  // The first files of the epoch log their fill; the rest are silent.
  Constants::SetPlotsThisFile(Constants::InPlotSample(t));
  {
    TChain ch("events");
    ch.Add(tasks[t].path);
    if (ch.GetEntries() > 0)
      phase.fills[t] = FillRunScatters(runOrder[i], labelVec[i], &ch,
                                       fits[i].gate, fits[i].pure_beam);
  }
  std::lock_guard<std::mutex> lk(phase.merge_mutex);
  phase.filled[t] = kTRUE;
  const Int_t nTasks = Int_t(tasks.size());
  while (phase.next_merge < nTasks && phase.filled[phase.next_merge]) {
    MergeRunFill(phase.next_merge, phase.fills, phase.nBeamKept,
                 phase.totalGated, phase.totalSeen);
    phase.next_merge++;
  }
}

void StripSumScatter::FillScatters(const FileSet::GateGroups &groups) {
  const std::vector<Int_t> &runOrder = groups.order;

  AllocateScatters();
  // The variant counters, in the order every per-run fill uses.
  {
    const std::vector<std::pair<TString, TagThresholds>> variants =
        ThresholdVariants();
    m_variantNames.clear();
    m_taggedVar.assign(variants.size(),
                       std::vector<Long64_t>(m_tagged.size(), 0));
    m_normedVar.assign(variants.size(),
                       std::vector<Long64_t>(m_tagged.size(), 0));
    for (size_t v = 0; v < variants.size(); v++)
      m_variantNames.push_back(variants[v].first);
  }

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

  // Phase 1: beam ellipses and the strip gate, one task per run.
  std::vector<SingleRunFitResult> fits(nRuns);
  RunIndexedParallel(nRuns, n_workers, [&](Int_t i) {
    FitGateRun(i, fits, runOrder, labelVec, chainVec);
  });

  // Phase 2: fill, split per file, not per run: a CoMPASS run is hundreds of
  // subfiles in one chain (SOLARIS: one per run); tasks carry their gates.
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

  // The workers fill their file and merge the finished prefixes in task
  // order, so the pure-beam budget stays identical to sequential.
  FillPhase phase;
  phase.fills.resize(nTasks);
  phase.filled.assign(nTasks, kFALSE);
  RunIndexedParallel(nTasks, fill_workers, [&](Int_t t) {
    FillScatterFile(t, tasks, runOrder, labelVec, fits, phase);
  });
  m_nSeen = phase.totalSeen;
  std::cout << "strip-sum-scatter: " << phase.totalGated
            << " reaction-tagged of " << phase.totalSeen << " events ("
            << m_nNormed << " past every pre-tag cut); reservoir "
            << m_reservoir.size() << std::endl;
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
  // Which conditions are on, from the same list the diagram draws.
  Bool_t active[kNTagCuts], pre_active[kNPreCuts];
  for (Int_t c = 0; c < kNTagCuts; c++)
    active[c] = kTRUE;
  for (Int_t c = 0; c < kNPreCuts; c++)
    pre_active[c] = kTRUE;
  const std::vector<SelectionStep> steps = DescribeSelection();
  for (size_t i = 0; i < steps.size(); i++) {
    if (steps[i].tag_cut >= 0)
      active[steps[i].tag_cut] = steps[i].on;
    if (steps[i].pre_cut >= 0)
      pre_active[steps[i].pre_cut] = steps[i].on;
  }

  std::ostringstream out;
  out << Form("strip-sum-scatter tag cuts: %s, %lld events seen",
              Paths::DatasetName().Data(), m_nSeen)
      << std::endl;
  out << "Sequential counts: each condition counts the events that passed "
         "everything before it."
      << std::endl;
  out << std::endl;
  out << "Event-level cuts (before any reaction is asked about):" << std::endl;
  for (Int_t c = 1; c < kNPreCuts; c++) {
    if (!pre_active[c])
      out << Form("  %-14s off", kPreCutName[c]) << std::endl;
    else
      out << Form("  %-14s %12lld  (%6.2f%% of seen)", kPreCutName[c],
                  m_preCounts[c],
                  m_nSeen > 0 ? 100.0 * m_preCounts[c] / m_nSeen : 0.0)
          << std::endl;
  }
  out << Form("  %-14s %12lld  (%6.2f%% of seen)", kPreCutName[kPrePass],
              m_preCounts[kPrePass],
              m_nSeen > 0 ? 100.0 * m_preCounts[kPrePass] / m_nSeen : 0.0)
      << std::endl;
  out << std::endl;

  out << "Tag conditions per reaction strip (events past the event-level "
         "cuts; % of those):"
      << std::endl;
  out << "reac";
  for (Int_t c = 1; c < kNTagCuts; c++)
    out << Form(" | %13s", kTagCutName[c]);
  out << Form(" | %13s", kTagCutName[kTagPass]) << std::endl;
  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    const Long64_t *row = &m_cutCounts[ReacIndex(reac) * kNTagCuts];
    Long64_t total = 0;
    for (Int_t c = 0; c < kNTagCuts; c++)
      total += row[c];
    out << Form("%4d", reac);
    for (Int_t c = 1; c < kNTagCuts; c++) {
      if (!active[c])
        out << Form(" | %13s", "off");
      else
        out << Form(" | %13s", Form("%lld (%.2f%%)", row[c],
                                    total > 0 ? 100.0 * row[c] / total : 0.0));
    }
    out << Form(" | %13s",
                Form("%lld (%.3f%%)", row[kTagPass],
                     total > 0 ? 100.0 * row[kTagPass] / total : 0.0))
        << std::endl;
  }
  std::cout << out.str();

  // The cut variation: tagged counts per strip under each shifted
  // threshold, beside the nominal, so the sensitivity is readable.
  if (!m_variantNames.empty()) {
    out << std::endl;
    out << "Cut variation (tagged events per strip with one threshold "
           "shifted; nominal first):"
        << std::endl;
    out << "reac |       nominal";
    for (size_t v = 0; v < m_variantNames.size(); v++)
      out << Form(" | %9s", m_variantNames[v].Data());
    out << std::endl;
    for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
      out << Form("%4d | %13lld", reac, m_tagged[ReacIndex(reac)]);
      for (size_t v = 0; v < m_variantNames.size(); v++)
        out << Form(" | %9lld", m_taggedVar[v][ReacIndex(reac)]);
      out << std::endl;
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
  f << out.str();
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

  // The selection as configured, as a block diagram: needs no data, so it
  // is on disk whatever happens below.
  SelectionDiagram::Write(
      DescribeSelection(), Paths::ResultsDir() + "/plots/strip_sum_scatter",
      Paths::DatasetName() +
          (Constants::GetActiveEpoch() ? " " + Constants::GetActiveEpoch()->name
                                       : TString("")) +
          " event selection");

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
  // first: the first group's ellipses select its beam events, the reference.
  Double_t strip_mean[18], strip_sigma[18];
  {
    Constants::SetPlotsThisFile(Constants::InPlotSample(0));
    SingleRunFitResult first =
        FitRunGates(run_order[0], first_label, chain_by_run[run_order[0]]);
    if (!first.pure_beam.ok ||
        !MeasureBeamNoise(chain_by_run[run_order[0]], first.pure_beam,
                          strip_mean, strip_sigma)) {
      std::cerr << "strip-sum-scatter: cannot measure the beam on "
                << first_label << std::endl;
      for (Int_t i = 0; i < Int_t(run_order.size()); i++)
        delete chain_by_run[run_order[i]];
      return kFALSE;
    }
    m_simStudyGate = first.gate;
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

  // Optional sim overlays, and the selection run over the sim.
  if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.RERUN_SIM) {
    SimOverlay();
    SimTraceOverlay();
    SimTagReport();
  }

  // Optional (a,p) plane study, scored on the simulated populations.
  if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.PLANE_STUDY)
    PlaneStudy();

  // Optional template fit of the tagged data against the sim components.
  if (Constants::cfg.STRIP_SUM_SCATTER_CONFIG.TEMPLATE_FIT)
    TemplateFit();

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
  // Same events Savitzky-Golay smoothed, as the strip-sum-scatter overlay
  // draws them, under the same option.
  const Bool_t kSkipSg =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.SKIP_SAVGOL_PLOTS;
  const std::vector<TGraph *> none;

  // The beam sample once: it is the same band under every strip's figure.
  std::vector<TGraph *> tr_beam, tr_beam_adc, tr_beam_sg;
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
    if (!kSkipSg)
      tr_beam_sg.push_back(SmoothedTraceFromTotal(td));
  }

  for (Int_t reac = kReacMin; reac <= kReacMax; reac++) {
    const UInt_t bit = (1u << ReacIndex(reac));
    std::vector<TGraph *> tr_tag, tr_tag_adc, tr_tag_sg;
    // Every tagged event at this strip, uncapped: in the all-tagged mode
    // the count IS these traces, so the figure shows all of them.
    for (Int_t k = 0; k < Int_t(m_reservoir.size()); k++) {
      const TraceEvt &e = m_reservoir[k];
      if (e.beam_flat || !(e.reac_mask & bit))
        continue;
      Double_t td[18], td_adc[18];
      e.Totals(td);
      e.TotalsAdc(td_adc);
      tr_tag.push_back(TraceFromTotal(td));
      tr_tag_adc.push_back(TraceFromTotal(td_adc));
      if (!kSkipSg)
        tr_tag_sg.push_back(SmoothedTraceFromTotal(td));
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
    if (!kSkipSg) {
      DrawRegionTraces(Form("all_tagged_traces_reac%d_sg", reac), subdir,
                       tr_beam_sg, none, tr_tag_sg, 0.6, 1.6, "#DeltaE [a.u.]",
                       kFALSE, "Tagged");
      if (kMeanTraces)
        DrawRegionMeanTraces(Form("all_tagged_mean_traces_reac%d_sg", reac),
                             subdir, tr_beam_sg, none, tr_tag_sg, 0.7, 1.3,
                             "#DeltaE [a.u.]");
    }
    for (Int_t i = 0; i < Int_t(tr_tag.size()); i++)
      delete tr_tag[i];
    for (Int_t i = 0; i < Int_t(tr_tag_adc.size()); i++)
      delete tr_tag_adc[i];
    for (Int_t i = 0; i < Int_t(tr_tag_sg.size()); i++)
      delete tr_tag_sg[i];
  }
  for (Int_t i = 0; i < Int_t(tr_beam.size()); i++)
    delete tr_beam[i];
  for (Int_t i = 0; i < Int_t(tr_beam_adc.size()); i++)
    delete tr_beam_adc[i];
  for (Int_t i = 0; i < Int_t(tr_beam_sg.size()); i++)
    delete tr_beam_sg[i];
}
