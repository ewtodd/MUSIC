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
  Float_t total[18];     ///< Calibrated per-strip totals.
  Float_t total_adc[18]; ///< Raw, un-normalised ADC sum per strip.
  /// @name Split-strip halves, in calibrated units
  /// Kept separately and independently of `IGNORE_SHORT_STRIPS`, which is a
  /// decode-time switch that zeroes one half of `EnergyView::left`/`right`.
  /// Reconstructed from the raw ADC and the per-channel gains, so one
  /// calibration can be rendered under either decode.
  /// @{
  Float_t long_au[18];  ///< Long end of each split strip.
  Float_t short_au[18]; ///< Short end.
  /// @}

  UInt_t reac_mask; ///< Bit per reaction strip this event was tagged at.
  Bool_t beam_flat; ///< Whether the trace looked flat, i.e. beam-like.
  Int_t both_mult;  ///< Split strips (1-16) with both ends above threshold.
  /// Timestamp of the grid hit that seeded the event, from the events tree.
  /// Unique per event and independent of how runs are split into files, so a
  /// cached event joins back to its source record. Zero when the source
  /// predates the seed-timestamp branch.
  ULong64_t seed_ts;
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
   * The entry point called from `main_strip_sum_scatter.cpp`.
   */
  void Run();

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
  ///         `POST_WINDOW_LAST_STRIP` / `POST_WINDOW_STRIPS` rule.
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

  /// @brief Sigma of the strip-to-strip difference `total[s] - total[s-1]`.
  /// @param strip Strip index.
  static Double_t JumpSigma(Int_t strip);
  /// @brief Sigma of a strip's own deposit.
  /// @param strip Strip index.
  static Double_t StripSigma(Int_t strip);
  /// @brief Install the jump sigmas. Call before tagging starts.
  /// @param sigma 18 values, one per strip.
  static void SetJumpSigma(const Double_t *sigma);
  /// @brief Install the per-strip sigmas. Call before tagging starts.
  /// @param sigma 18 values, one per strip.
  static void SetStripSigma(const Double_t *sigma);
  /// @brief Minimum jump for a tag: `REAC_JUMP_NSIGMA * JumpSigma(reac)`.
  /// @param reac Reaction strip index.
  static Double_t JumpMin(Int_t reac);
  /// @}

private:
  static Double_t s_jumpSigma[18];
  static Double_t s_stripSigma[18];
  std::map<Int_t, TH2F *> m_scatter;
  std::vector<TraceEvt> m_reservoir;
  // Normalization counts, merged over every run and persisted in the cache so
  // a cross section can be taken from it without a second pass over the data.
  Long64_t m_nSeen;
  Long64_t m_nNormed;
  std::vector<Long64_t> m_normedAt;
  std::vector<Long64_t> m_tagged;
  Double_t m_yLo[64];
  Double_t m_yHi[64];

  const TString kSimCacheName = "StripSumScatter_simcache.root";

  static Int_t ReacIndex(Int_t reac);

  Bool_t TryLoadCache(const TString &cacheName, const TString &fingerprint);
  void WriteCache(const TString &cacheName, const TString &fingerprint);

  // Fresh, empty scatters for the configured plane and build range.
  void AllocateScatters();
  // Fill the scatters from the reservoir alone. The reservoir keeps every
  // tagged event, so a plane or axes change needs only this and not a
  // 25-minute pass over the events files.
  void ReprojectFromReservoir();
  void FillScatters(const std::vector<Int_t> &runOrder,
                    std::map<Int_t, TChain *> &chains);

  void PlotScatters();

  // Diagnostic (PLOT_PARITY_REJECTED_GRID): grid #DeltaE of events that pass
  // the cheap pre-tag cuts and are then rejected by the parity cut. Saves both
  // an a.u. view (decoded grid) and an ADC view (raw trigger channel), each on
  // a log-y axis. Reads the Grid branch directly (a self-contained pass), so it
  // works regardless of whether the scatter cache has already been filled.
  void PlotParityRejectedGrid(const std::vector<Int_t> &run_order,
                              std::map<Int_t, TChain *> &chains);

  void InteractiveOverlay(Int_t reac);

  static void EnableEventBranches(TChain *chain);
  static Bool_t AllStripsFired(const EnergyView &ev);
  static Bool_t IsPureBeam(const EnergyView &ev, const BeamEllipses &be);
  // Sigma-clipped width of each strip-to-strip difference and of each strip's
  // deposit over a capped sample of `chain`, after the cheap pre-tag cuts.
  // False when too few events survive to measure them.
  static Bool_t MeasureBeamNoise(TChain *chain, Double_t *jump_sigma,
                                 Double_t *strip_sigma);
  // Strips 1..reac-1 within BEAM_UPSTREAM_NSIGMA of the beam, or the
  // requirement is off. Shared by the tag and its per-strip denominator.
  static Bool_t BeamUpstreamOf(const EnergyView &ev, Int_t reac);
  static Bool_t IsPileup(const EnergyView &ev);
  static Bool_t IsNoise(const EnergyView &ev);
  static Bool_t IsOffbeam(const EnergyView &ev);
  static Bool_t IsParityAsymmetric(const EnergyView &ev);
  static Double_t SumRange(const Double_t *total, Int_t lo, Int_t hi);
  static std::vector<GateSpec> ActiveGates();

  static Bool_t PassesGate(const BeamFit2D &gate, const EnergyView &ev,
                           Int_t sx, Int_t sy);

  static BeamFit2D FindBeamGate(TChain *chain, Int_t sx, Int_t sy,
                                const std::vector<GateSpec> &prior_specs,
                                const std::vector<BeamFit2D> &prior_gates,
                                const TString &tag, const TString &subdir);

  static void DrawTraceSet(const std::vector<TGraph *> &traces, Int_t color);
  void DrawAltDecodeRegionTraces(Int_t reac, TCutG *cutAn, TCutG *cutAa);
  static TGraph *TraceFromTotal(const Float_t *total);
  static void DrawRegionTraces(const TString &save_name, const TString &subdir,
                               const std::vector<TGraph *> &beam,
                               const std::vector<TGraph *> &aa,
                               const std::vector<TGraph *> &an, Double_t y_min,
                               Double_t y_max, const char *y_title);

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
  static SingleRunFitResult
  FitRunGates(Int_t run, TChain *chain,
              const std::vector<GateSpec> &activeGates);
  static SingleRunFillResult FillRunScatters(
      Int_t run, TChain *chain, const std::vector<GateSpec> &activeGates,
      const std::vector<BeamFit2D> &runGates, const BeamEllipses &runBeam);
  static TCutG *PromptCut(TCanvas *c, const char *name, const char *label);
  static void SaveRegionCuts(Int_t reac, TCutG *cut_an, TCutG *cut_aa);
  static TCutG *LoadRegionCut(const char *name, Int_t reac);

  static void SmoothTrace(const Double_t *in, Double_t *out, Int_t width);

  // Savitzky-Golay smoothing: 3rd-degree polynomial, half-window of 2
  // (5-point convolution). Uses standard SG coefficients [-3,12,17,12,-3]/35.
  // At edges, the window shrinks and coefficients are renormalised.
  static void SavitzkyGolay(const Double_t *in, Double_t *out);

  // CFD-style trigger finder: locate the first strip whose beam-subtracted
  // signal (td[s]-1) exceeds both a fraction of the trace peak and a multiple
  // of the beam sigma. Returns the strip index, or -1 if no trigger fires.
  static Int_t FindTrigger(const Double_t *td, const Double_t *base,
                           Double_t beam_sigma);

  // Build a TGraph from Savitzky-Golay-smoothed per-strip totals. Input is
  // the raw normed array; smoothing is applied internally before graph build.
  static TGraph *SmoothedTraceFromTotal(const Float_t *total);

  void ClusterVarHists(Int_t reac, TCutG *cut_aa, TCutG *cut_an,
                       const TString &subdir);

  static TString
  SimFingerprint(const std::vector<RemixSim::SimFileSpec> &specs);
  static TString BuildFingerprint(const std::vector<Int_t> &run_order,
                                  std::map<Int_t, TChain *> &chains);
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

private:
  static Bool_t SimBeamGains(Double_t *gain);
  static void SimTotal(const Float_t *left, const Float_t *right,
                       const Double_t *gain, Double_t *total);
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
