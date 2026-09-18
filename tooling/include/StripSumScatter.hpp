#ifndef STRIP_SUM_SCATTER_HPP
#define STRIP_SUM_SCATTER_HPP

#include "BeamFit2D.hpp"
#include "Constants.hpp"
#include "EventsSummary.hpp"
#include "FileSet.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "Normalization.hpp"
#include "PlottingUtils.hpp"
#include "RemixSim.hpp"
#include <Rtypes.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TCutG.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TKey.h>
#include <TLegend.h>
#include <TMath.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

/**
 * @file StripSumScatter.hpp
 * @brief The reaction search: strip-sum scatters, beam gating and tagging.
 *
 * Builds, per reaction strip, a two-dimensional scatter in which beam-like and
 * reaction events separate, then tags reactions in it. The scatters and a
 * reservoir of tagged events are cached, so retuning a display window or a
 * region cut costs a reprojection rather than another pass over the data.
 */

/**
 * @brief The conditions of the reaction tag, in the order they are applied.
 *
 * `StripSumScatter::RejectReason` returns the first one an event fails at a
 * reaction strip, so the per-condition counts are sequential: each counts
 * events that passed everything before it. `kTagPass` is a tag.
 */
enum TagCut {
  kTagPass = 0,
  kCutAllStrips, ///< A strip did not fire.
  kCutUpstream,  ///< A strip before the reaction was not beam-like.
  kCutJump,      ///< Jump at the reaction strip below the gate.
  kCutReacLevel, ///< Reaction-strip deposit below 1 + the gate.
  kCutSmooth,    ///< A post-reaction step above SMOOTHNESS_NSIGMA.
  kCutTailRise,  ///< A rise in the tail (TAIL_RISE_NSIGMA).
  kCutRerise,    ///< Back at the beam, then above it again (TAIL_RERISE_*).
  kCutPostAbove, ///< The excess did not persist (POST_ABOVE_*).
  kCutEndStrip,  ///< The end strip not END_STRIP_NSIGMA below the beam.
  kCutCliff,     ///< The fall happened in the last step (TAIL_CLIFF_*).
  kNTagCuts
};
/// @brief Short labels for TagCut, indexed by it.
extern const char *const kTagCutName[kNTagCuts];

/**
 * @brief Every threshold the tag applies, as one value set.
 *
 * The nominal set is read off `StripSumScatterConfig`
 * (StripSumScatter::NominalThresholds); the cut-variation systematic shifts
 * one field at a time (StripSumScatter::ThresholdVariants) and re-runs the
 * same tag, so the count's sensitivity to each threshold is measured by the
 * code that applies it.
 */
struct TagThresholds {
  Bool_t require_upstream = kTRUE;
  Double_t upstream_nsigma = 0.0;
  Double_t jump_nsigma = 0.0;
  Double_t smooth_nsigma = 0.0;
  Int_t tail_fall_from_strip = 0;
  Double_t tail_rise_nsigma = 0.0;
  Double_t tail_return_nsigma = 0.0;
  Double_t tail_rerise_nsigma = 0.0;
  Double_t post_above_nsigma = 0.0;
  Int_t post_above_strips = 0;
  Double_t end_strip_nsigma = 0.0;
  Double_t cliff_max = 0.0;
};

/**
 * @brief The event-level cuts applied before any reaction is asked about, in
 * order; sequential like TagCut. `kPrePass` reached the tag.
 */
enum PreCut {
  kPrePass = 0,
  kPreGate,     ///< Failed a beam gate.
  kPrePileup,   ///< Pileup.
  kPreNoise,    ///< Noise.
  kPreBothMult, ///< Both-ends multiplicity (BOTH_MULT_MAX).
  kNPreCuts
};
/// @brief Short labels for PreCut, indexed by it.
extern const char *const kPreCutName[kNPreCuts];

/// @brief A pair of strips whose sums form one classification plane.
struct GateSpec {
  Int_t sx; ///< Strip whose sum forms the x axis.
  Int_t sy; ///< Strip whose sum forms the y axis.
};

/**
 * @brief The fixed window every strip-sum scatter is built over.
 *
 * Building to a fixed range is what lets the display windows be retuned without
 * refilling.
 *
 * The axes get different ceilings because they sum different numbers of strips:
 * x runs over `X_LO..X_HI` (16 strips, so it reaches about 20), while y covers
 * only `POST_TRIGGER_SUM_STRIPS` after the trigger (6, so it never approaches
 * 40). Giving y its own lower ceiling doubles its resolution at the same bin
 * count — which matters, because y is the axis the reaction populations
 * separate along, and its display windows are much narrower than x's.
 */
namespace ScatterBuildRange {
const Double_t kXMin = 0.0;  ///< Lower x bound of the build window.
const Double_t kXMax = 40.0; ///< Upper x bound; x sums 16 strips.
const Double_t kYMin = 0.0;  ///< Lower y bound.
const Double_t kYMax =
    20.0; ///< Upper y bound; y sums the post-trigger strips only.
} // namespace ScatterBuildRange

/**
 * @brief The classification ellipses defining a pure-beam event.
 *
 * An entrance ellipse and an exit ellipse. An event is pure beam only if it
 * passes **both** — entering like beam and leaving like beam.
 */
struct BeamEllipses {
  BeamFit2D s0_s1;    ///< Entrance ellipse on strips 0 and 1.
  BeamFit2D s1_s2;    ///< Alternative entrance ellipse, per `PURE_BEAM_GATE`.
  BeamFit2D s16_s17;  ///< Exit ellipse on strips 16 and 17.
  BeamFit2D s15_s16;  ///< Alternative exit ellipse.
  Bool_t ok;          ///< Whether the fits succeeded.
  Bool_t use_s15_s16; ///< Which exit ellipse is in force.
};

/**
 * @brief One tagged event, kept in the reservoir.
 *
 * The reservoir holds every tagged event, so a change to the scatter plane or
 * its axes needs only a reprojection rather than another full pass over the
 * events files.
 */
struct TraceEvt {
  /// @name Calibrated ends, each already carrying its strip's alignment
  /// @{
  Float_t leftdE[16];  ///< Left end of strips 1-16, index `strip - 1`.
  Float_t rightdE[16]; ///< Right end of strips 1-16, index `strip - 1`.
  Float_t strip0dE;    ///< Strip 0, unsegmented.
  Float_t strip17dE;   ///< Strip 17, unsegmented.
  /// @}
  /// @name The same, raw ADC
  /// @{
  Float_t leftdE_adc[16];
  Float_t rightdE_adc[16];
  Float_t strip0_adc;
  Float_t strip17_adc;
  /// @}
  UInt_t reac_mask; ///< Bit per reaction strip this event was tagged at.
  Bool_t beam_flat; ///< Whether the trace looked flat, i.e. beam-like.
  Int_t both_mult;  ///< Split strips (1-16) with both ends above threshold.
  /// Timestamp of the grid hit that seeded the event, from the events tree.
  /// Unique per event and independent of how runs are split into files, so a
  /// cached event joins back to its source record. Zero when the source
  /// predates the seed-timestamp branch.
  ULong64_t seed_ts;

  /// @brief A strip's deposit: the sum of its ends, or the unsegmented value.
  Double_t Total(Int_t strip) const {
    if (strip <= 0)
      return strip0dE;
    if (strip >= 17)
      return strip17dE;
    return Double_t(leftdE[strip - 1]) + Double_t(rightdE[strip - 1]);
  }
  /// @brief A strip's raw ADC deposit, summed the same way.
  Double_t TotalAdc(Int_t strip) const {
    if (strip <= 0)
      return strip0_adc;
    if (strip >= 17)
      return strip17_adc;
    return Double_t(leftdE_adc[strip - 1]) + Double_t(rightdE_adc[strip - 1]);
  }
  /// @brief Every strip's deposit, for code that wants an array of 18.
  void Totals(Double_t *out) const {
    for (Int_t s = 0; s < 18; s++)
      out[s] = Total(s);
  }
  /// @brief Every strip's raw ADC deposit, for code that wants an array of 18.
  void TotalsAdc(Double_t *out) const {
    for (Int_t s = 0; s < 18; s++)
      out[s] = TotalAdc(s);
  }
};

/**
 * @brief One run's fitted beam gates.
 *
 * Both fill phases run one worker per run and merge afterwards **in run
 * order**, which is what keeps the threaded result bit-identical to the
 * sequential one.
 */
struct SingleRunFitResult {
  BeamEllipses pure_beam;              ///< Entrance and exit ellipses.
  std::vector<BeamFit2D> series_gates; ///< One gate per active GateSpec.
  Bool_t ok;                           ///< Whether the fits succeeded.
  /// @brief Construct unfitted.
  SingleRunFitResult() : ok(kFALSE) {}
};

/**
 * @brief One run's filled scatters, reservoir and normalisation counts.
 */
struct SingleRunFillResult {
  /// Private clones, one per reaction strip. **Owned by the caller after the
  /// merge.**
  std::vector<TH2F *> scatters;
  std::vector<TraceEvt> reservoir; ///< Tagged events from this run.
  Long64_t gated;                  ///< Events passing the beam gates.
  Long64_t seen;                   ///< Events examined.
  /// Events surviving every cut applied before reaction tagging: the beam
  /// gates, the pileup, noise and off-beam rejection, and the both-ends
  /// multiplicity. Every tagged event passed exactly this selection, which is
  /// what makes it the denominator in which those efficiencies cancel.
  Long64_t normed;
  /// Per-strip denominator: beam particles that reached that strip under
  /// exactly the conditions a reaction there would have had to satisfy — every
  /// strip fired, and strips `1..reac-1` beam-like.
  ///
  /// Those conditions are strip-dependent: strip 5 needs four upstream strips
  /// to look like beam where strip 3 needs two. A single flat denominator would
  /// leave that efficiency uncancelled and bias the excitation function.
  std::vector<Long64_t> normed_at;
  /// Events tagged at each reaction strip, indexed the same way.
  std::vector<Long64_t> tagged;
  /// Events tagged under each threshold variant, `[variant][ReacIndex]`, in
  /// the order of StripSumScatter::ThresholdVariants().
  std::vector<std::vector<Long64_t>> tagged_var;
  /// Sequential per-condition counts, `[ReacIndex(reac) * kNTagCuts + cut]`,
  /// for the tag-cut report. `kTagPass` entries equal `tagged`.
  std::vector<Long64_t> cut_counts;
  /// Sequential event-level counts, indexed by PreCut.
  std::vector<Long64_t> pre_counts;
  /// @brief Construct with zeroed counters.
  SingleRunFillResult() : gated(0), seen(0), normed(0) {}
};

/// @brief One simulated population to overlay on the data.
struct SimPop {
  TString file;  ///< Simulation ROOT file.
  TString label; ///< Legend label.
};

/**
 * @brief Builds the strip-sum scatters, tags reactions, and draws the results.
 *
 * The entry point behind the `strip-sum-scatter` binary. Fills per-strip
 * scatters from the events files, caches them alongside a reservoir of tagged
 * events, and draws the regions and trace overlays.
 *
 * @note The cache carries a fingerprint of every configuration knob that
 *       changes its contents — filters, gates, strip spans — so editing any of
 *       them rebuilds it automatically on the next run. There is no
 *       cache-clearing step to remember.
 */
class StripSumScatter {
public:
  /// @brief Construct with empty scatters and no cache loaded.
  StripSumScatter();
  /// @brief Frees the scatters and the reservoir.
  ~StripSumScatter();

  /**
   * @brief Build or load the scatters, tag reactions, and draw everything.
   *
   * The entry point called from `main_strip_sum_scatter.cpp`: Prepare(),
   * then the optional simulation overlays and the interactive region-trace
   * overlay for the candidate strip.
   */
  void Run();

  /**
   * @brief Measure the beam, then load the scatter cache or build it.
   *
   * Everything Run() does before any overlay: the first run's beam mean and
   * width per strip, the cache fingerprint, a load of a matching cache or a
   * fill and write of a fresh one, and the batch scatter figures. Leaves the
   * scatters and the trace reservoir in memory. compute-regions calls this in
   * the all-tagged mode so strip-sum-scatter itself never has to be run.
   *
   * @return `kFALSE` when there are no runs or the beam cannot be measured.
   */
  Bool_t Prepare();

  /**
   * @brief Region-trace overlays of every tagged event, one figure per
   *        reaction strip, against the beam band.
   *
   * The all-tagged mode's counterpart of the interactive overlay: no region
   * cut, and every event tagged at the strip is drawn with no cap (the
   * `TRACES_PER_CLASS` cap applies to the beam sample only), in calibrated
   * units, in raw ADC when `PLOT_ADC_TRACES`, and as mean traces when
   * `PLOT_REGION_MEAN_TRACES`. Needs the reservoir, so call after Prepare().
   *
   * @param subdir Plot subdirectory the figures go to.
   */
  void DrawAllTaggedTraces(const TString &subdir);

  /**
   * @brief Filename of the scatter cache for this configuration.
   * @note `compute-regions` reads the same file, so the name is defined here
   *       and nowhere else.
   */
  static TString CacheName();

  /**
   * @brief Where an event sits in the scatter plane for a given reaction strip.
   *
   * x is the sum over `X_LO..X_HI`; y is the sum over the post-trigger window.
   *
   * @param total Calibrated per-strip totals for the event.
   * @param reac  Reaction strip index.
   * @param[out] x Plane x coordinate.
   * @param[out] y Plane y coordinate.
   *
   * @note Every consumer of a region cut goes through this, which is what stops
   *       the fill, the overlay and the cross section from ever disagreeing
   *       about where an event sits.
   */
  static void PlaneXY(const Double_t *total, Int_t reac, Double_t &x,
                      Double_t &y);
  /// @brief First strip of the post-trigger window summed onto y.
  /// @param reac Reaction strip index.
  /// @return The strip index, inclusive, per the `POST_TRIGGER_SUM_STRIPS` /
  ///         `POST_WINDOW_LAST_STRIP` rule.
  static Int_t YLoOf(Int_t reac);
  /// @brief Last strip of the post-trigger window, inclusive.
  /// @param reac Reaction strip index.
  static Int_t YHiOf(Int_t reac);
  /**
   * @name Measured beam noise
   *
   * Measured once from the data and stamped into the cache. Set before any
   * tagging begins and read-only afterwards, which is what makes them safe to
   * share across the worker threads without synchronisation.
   * @{
   */

  /// @brief Width of the beam's deposit on a strip, from the gated beam
  ///        sample of the first run (MeasureBeamNoise).
  /// @param strip Strip index.
  static Double_t StripSigma(Int_t strip);
  /// @brief Mean of the beam's deposit on a strip over the same sample; the
  ///        level every cut measures from. 1.0 when none was installed.
  /// @param strip Strip index.
  static Double_t StripMean(Int_t strip);
  /// @brief Install the per-strip beam means and sigmas. Call before tagging
  ///        starts.
  /// @param mean  18 values, one per strip; null means 1.0 everywhere.
  /// @param sigma 18 values, one per strip.
  static void SetStripNoise(const Double_t *mean, const Double_t *sigma);
  /// @brief Minimum jump for a tag: `REAC_JUMP_NSIGMA * StripSigma(reac)`.
  /// @param reac Reaction strip index.
  static Double_t JumpMin(Int_t reac);
  /// @}

private:
  static Double_t s_stripSigma[18];
  static Double_t s_stripMean[18];
  std::map<Int_t, TH2F *> m_scatter;
  std::vector<TraceEvt> m_reservoir;
  // Normalization counts, merged over every run and persisted in the cache so
  // a cross section can be taken from it without a second pass over the data.
  Long64_t m_nSeen;
  Long64_t m_nNormed;
  std::vector<Long64_t> m_normedAt;
  std::vector<Long64_t> m_tagged;
  /// Tagged counts per threshold variant, `[variant][ReacIndex]`; the
  /// variant names alongside, persisted with the cache for the cross
  /// section's cut-variation systematic.
  std::vector<std::vector<Long64_t>> m_taggedVar;
  std::vector<TString> m_variantNames;
  // Per-condition tag counts from the last fill, for the cut report; not
  // cached, so a cache load leaves them empty and writes no report.
  std::vector<Long64_t> m_cutCounts;
  std::vector<Long64_t> m_preCounts;
  /// Print the per-condition counts and write them to
  /// plots/strip_sum_scatter/tag_cuts.txt.
  void WriteCutReport() const;
  Double_t m_yLo[64];
  Double_t m_yHi[64];

  const TString kSimCacheName = "StripSumScatter_simcache.root";

  static Int_t ReacIndex(Int_t reac);

  Bool_t TryLoadCache(const TString &cacheName, const TString &fingerprint);
  void WriteCache(const TString &cacheName, const TString &fingerprint);

  // Fresh, empty scatters for the configured plane and build range.
  void AllocateScatters();
  /// Fill the scatters from the reservoir alone. The reservoir keeps every
  /// tagged event, so a plane or axes change needs only this and not a
  /// 25-minute pass over the events files.
  void ReprojectFromReservoir();
  void FillScatters(const FileSet::GateGroups &groups);

  void PlotScatters();

  void InteractiveOverlay(Int_t reac);

  static void EnableEventBranches(TChain *chain);
  static Bool_t AllStripsFired(const EnergyView &ev);
  static Bool_t IsPureBeam(const EnergyView &ev, const BeamEllipses &be);
  /// The beam's mean and width on every strip: over the events of `chain`
  /// that pass the entrance and exit ellipses `be` (IsPureBeam), up to a cap,
  /// each strip's sigma-clipped mean and width. Every sigma-unit cut, level
  /// or step, resolves through these. False when too few events pass the
  /// gate to measure them.
  static Bool_t MeasureBeamNoise(TChain *chain, const BeamEllipses &be,
                                 Double_t *strip_mean, Double_t *strip_sigma);
  // Strips 1..reac-1 within BEAM_UPSTREAM_NSIGMA of the beam, or the
  // requirement is off. Shared by the tag and its per-strip denominator.
  static Bool_t BeamUpstreamOf(const EnergyView &ev, Int_t reac);
  static Bool_t BeamUpstreamOf(const EnergyView &ev, Int_t reac,
                               const TagThresholds &T);
  static Bool_t IsPileup(const EnergyView &ev);
  static Bool_t IsNoise(const EnergyView &ev);
  static Double_t SumRange(const Double_t *total, Int_t lo, Int_t hi);
  static std::vector<GateSpec> ActiveGates();

  static Bool_t PassesGate(const BeamFit2D &gate, const EnergyView &ev,
                           Int_t sx, Int_t sy);

  static BeamFit2D FindBeamGate(TChain *chain, Int_t sx, Int_t sy,
                                const std::vector<GateSpec> &prior_specs,
                                const std::vector<BeamFit2D> &prior_gates,
                                const TString &tag, const TString &subdir);

  static void DrawTraceSet(const std::vector<TGraph *> &traces, Int_t color);
  static TGraph *TraceFromTotal(const Double_t *total);
  /// Overlay of sampled traces per region. The beam is drawn as its mean
  /// with a +-1 sigma band rather than as individual traces: the sigma is
  /// the cache-measured per-strip beam spread (StripSigma, scaled by the
  /// sample mean so it applies in ADC too) when `beam_sigma_measured`, else
  /// the sample's own RMS (for smoothed traces, whose noise is narrower).
  /// `an_label` renames the red class in the legend (the all-tagged overlay
  /// draws every tagged event there); an empty class is left out of it.
  static void DrawRegionTraces(const TString &save_name, const TString &subdir,
                               const std::vector<TGraph *> &beam,
                               const std::vector<TGraph *> &aa,
                               const std::vector<TGraph *> &an, Double_t y_min,
                               Double_t y_max, const char *y_title,
                               Bool_t beam_sigma_measured = kTRUE,
                               const char *an_label = nullptr);

  static void DrawRegionMeanTraces(const TString &save_name,
                                   const TString &subdir,
                                   const std::vector<TGraph *> &beam,
                                   const std::vector<TGraph *> &aa,
                                   const std::vector<TGraph *> &an,
                                   Double_t y_min, Double_t y_max,
                                   const char *y_title);

  static void TraceYRange(const std::vector<TGraph *> &beam,
                          const std::vector<TGraph *> &aa,
                          const std::vector<TGraph *> &an, Double_t &y_min,
                          Double_t &y_max);
  /// One gate group (FileSet::GateGroups): a run's chunks on SOLARIS, one
  /// subfile on CoMPASS. `label` names it in logs and plot folders.
  static SingleRunFitResult
  FitRunGates(Int_t key, const TString &label, TChain *chain,
              const std::vector<GateSpec> &activeGates);
  static SingleRunFillResult
  FillRunScatters(Int_t key, const TString &label, TChain *chain,
                  const std::vector<GateSpec> &activeGates,
                  const std::vector<BeamFit2D> &runGates,
                  const BeamEllipses &runBeam);
  static TCutG *PromptCut(TCanvas *c, const char *name, const char *label);
  static void SaveRegionCuts(Int_t reac, TCutG *cut_an, TCutG *cut_aa);
  static TCutG *LoadRegionCut(const char *name, Int_t reac);

  static void SmoothTrace(const Double_t *in, Double_t *out, Int_t width);

  /// Savitzky-Golay smoothing: 3rd-degree polynomial, half-window of 2
  /// (5-point convolution). Uses standard SG coefficients [-3,12,17,12,-3]/35.
  /// At edges, the window shrinks and coefficients are renormalised.
  static void SavitzkyGolay(const Double_t *in, Double_t *out);

  /// CFD-style trigger finder: locate the first strip whose beam-subtracted
  /// signal (td[s]-1) exceeds both a fraction of the trace peak and a multiple
  /// of the beam sigma. Returns the strip index, or -1 if no trigger fires.
  static Int_t FindTrigger(const Double_t *td, const Double_t *base,
                           Double_t beam_sigma);

  // Build a TGraph from Savitzky-Golay-smoothed per-strip totals. Input is
  // the raw normed array; smoothing is applied internally before graph build.
  static TGraph *SmoothedTraceFromTotal(const Double_t *total);

  void ClusterVarHists(Int_t reac, TCutG *cut_aa, TCutG *cut_an,
                       const TString &subdir);

  static TString
  SimFingerprint(const std::vector<RemixSim::SimFileSpec> &specs);
  static TString BuildFingerprint(const FileSet::GateGroups &groups);
  static void YBounds(Double_t *y_lo, Double_t *y_hi);
  static TString PrettyLabel(const TString &tag);

public:
  /**
   * @brief Whether an event is tagged as a reaction at a given strip.
   *
   * @param ev   Decoded event.
   * @param reac Reaction strip index.
   * @return `kTRUE` if the event is tagged there.
   *
   * @note Public because TagEfficiency pushes bootstrapped traces through this
   *       same tag. An efficiency measured against a different selection than
   *       the one that produced the count would not apply to it.
   */
  static Bool_t PassesReaction(const EnergyView &ev, Int_t reac);
  /**
   * @brief The tail-shape part of the tag, applied inside PassesReaction.
   *
   * The conditions the published 87Rb per-strip macros put on the strips
   * downstream of the reaction: smoothness continued to the last strip, a
   * monotonically falling tail, no return to the beam and persistence of the
   * excess. Each is
   * off unless its `StripSumScatterConfig` value is set, so a dataset that
   * sets none of them tags exactly as before.
   *
   * @param ev    Decoded event.
   * @param reac  Reaction strip index.
   * @return `kTRUE` if the tail passes every enabled condition.
   */
  static Bool_t PassesTail(const EnergyView &ev, Int_t reac);
  /**
   * @brief Which condition of the tag an event fails first at a strip.
   * @param ev    Decoded event.
   * @param reac  Reaction strip index.
   * @return `kTagPass` if tagged, else the first failing TagCut. PassesReaction
   *         is exactly `RejectReason(ev, reac) == kTagPass`.
   */
  static TagCut RejectReason(const EnergyView &ev, Int_t reac);
  /// @brief RejectReason under an explicit threshold set.
  static TagCut RejectReason(const EnergyView &ev, Int_t reac,
                             const TagThresholds &T);
  /// @brief The tail-shape part of RejectReason: `kTagPass` or the first
  ///        failing tail condition.
  static TagCut TailReason(const EnergyView &ev, Int_t reac);
  static TagCut TailReason(const EnergyView &ev, Int_t reac,
                           const TagThresholds &T);
  /// @brief The thresholds the configuration sets.
  static TagThresholds NominalThresholds();
  /// @brief The cut-variation set: each active threshold shifted up and down
  ///        by its configured step (`CUT_VARIATION_*`), one at a time, named
  ///        `<threshold>+` and `<threshold>-`. Empty when `CUT_VARIATION` is
  ///        off.
  static std::vector<std::pair<TString, TagThresholds>> ThresholdVariants();

private:
  static Bool_t SimBeamGains(Double_t *gain);
  // Run n indexed tasks on a pool of workers, pulling from a shared queue.
  static void RunIndexedParallel(Int_t n, Int_t workers,
                                 const std::function<void(Int_t)> &task);
  static void SimTotal(const RemixSim::Event &e, const Double_t *gain,
                       Double_t *total);
  static TGraph *SimPopScatter(const TString &file, Int_t reac,
                               const Double_t *gain, Long64_t max_points);
  static std::vector<TGraph *>
  SimPopTraces(const TString &file, const Double_t *gain, Long64_t max_traces);
  void SimTraceOverlay();
  static Bool_t LoadSimCache(const TString &fp,
                             std::map<Int_t, std::vector<TGraph *>> &by_strip);
  static void
  WriteSimCache(const TString &fp,
                const std::map<Int_t, std::vector<TGraph *>> &by_strip);
  void SimOverlay();
};

#endif
