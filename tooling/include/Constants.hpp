#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include "DedupStrategy.hpp"
#include "PulseHistoryGroups.hpp"
#include "RunEpoch.hpp"
#include "SlotLayout.hpp"
#include <Rtypes.h>
#include <RtypesCore.h>
#include <TString.h>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

/**
 * @file Constants.hpp
 * @brief Dataset configuration.
 *
 * Configuration is split in two. The struct definitions and the tooling-wide
 * defaults live here and in `Constants.cpp`; each dataset then overrides only
 * the fields it cares about in `analysis/<dataset>/config/Constants.cpp`, which
 * is compiled into that dataset's binaries.
 *
 * There is additionally an option for epochs. A dataset spanning several
 * acquisition eras declares a RunEpoch per era, and the `Active*()` accessors
 * read the active epoch when one is set and the flat block otherwise.
 */

/**
 * @def MUSIC_HOT_PATH_LOGGING
 * @brief Whether per-hit and per-event progress logging is compiled in.
 *
 * Compile-time so the branch leaves the hot loop entirely when off. Define it
 * to `1` on the compiler command line to enable.
 *
 * @note Logging outside hot loops (per-file and per-run summaries) is
 *       unconditional and not covered by this option.
 */
#ifndef MUSIC_HOT_PATH_LOGGING
#define MUSIC_HOT_PATH_LOGGING 0
#endif

/**
 * @brief Options which govern the reaction search in the 2D histograms
 * (strip-sum-scatter).
 *
 * There are two kinds of controls; ones which define which events populating
 * the histograms, and ones which define how the histograms are constructed.
 * Each field says which it is. The cache fingerprints them all, so rebuilds
 * should happen automatically if they are needed.
 */
/// An axis of the pure-beam entrance ellipse (PURE_BEAM_GATE_S1_GRID): the
/// Frisch grid over its beam anchor (`GainGrid`), so it sits at 1.0 like a
/// strip.
const Int_t GATE_AXIS_GRID = -1;

struct StripSumScatterConfig {
  /// The entrance ellipse of the pure-beam sample (beam-noise reference and
  /// beam_flat reservoir): strips 0 and 1, strips 1 and 2, or strip 1 and
  /// the grid.
  enum PureBeamGate {
    PURE_BEAM_GATE_S0_S1,
    PURE_BEAM_GATE_S1_S2,
    PURE_BEAM_GATE_S1_GRID
  };
  PureBeamGate PURE_BEAM_GATE;
  /// One tag per event. When several reaction strips pass the tag for the
  /// same event, the reaction is attributed to one of them: the first (most
  /// upstream) strip that passes, where the trace first leaves the beam; or
  /// the strip with the largest jump from its predecessor, in sigma. The
  /// others count as "other strip" in the cut report.
  enum TagResolve { TAG_RESOLVE_FIRST_STRIP, TAG_RESOLVE_LARGEST_JUMP };
  TagResolve TAG_RESOLVE;
  /// Mahalanobis level of the pure-beam entrance and exit ellipses; separate
  /// from GATE_NSIGMA so the beam sample stays narrow whatever the gate does.
  Double_t PURE_BEAM_NSIGMA;

  /// Number of strips summed onto the scatter y-axis after the trigger strip:
  /// y spans reac+1 .. min(reac+POST_TRIGGER_SUM_STRIPS,
  /// POST_WINDOW_LAST_STRIP)
  Int_t POST_TRIGGER_SUM_STRIPS;
  Int_t POST_WINDOW_LAST_STRIP;

  /// Cap on worker threads for the scatter fill.
  Int_t MAX_STRIP_SUM_WORKERS;

  Int_t REACTION_STRIP_MIN;
  Int_t REACTION_STRIP_MAX;

  /// Event-level smoothness: every strip-to-strip step of the trace, from
  /// strips 1-2 to the last pair, must be under SMOOTHNESS_NSIGMA times the
  /// later strip's measured beam spread (StripSumScatter::StripSigma), the
  /// single largest rise excepted. A reaction is one rise (the jump at the
  /// reaction strip) followed by gradual change; a spike, a partial second
  /// particle or a glitch shows more than one abrupt step. So the biggest
  /// rise is left to the jump condition and every other step, falls
  /// included, is held to this limit, which keeps it independent of
  /// REAC_JUMP_NSIGMA. Applied once per event before any strip is asked
  /// about, so it enters the beam denominator too. A step is the difference
  /// of two strips, so its noise is about 1.4 spreads. 0 or less = off. Part
  /// of the cut variation.
  Double_t SMOOTHNESS_NSIGMA;

  /// From TAIL_FALL_FROM_STRIP to the last strip, allow no
  /// strip-to-strip rise above TAIL_RISE_NSIGMA x StripSigma of the strip.
  /// FROM_STRIP 0 or less = off.
  Int_t TAIL_FALL_FROM_STRIP;
  Double_t TAIL_RISE_NSIGMA;

  /// After the post-window peak, once the trace has come
  /// back down to within TAIL_RETURN_NSIGMA x StripSigma above the beam it must
  /// not rise above TAIL_RERISE_NSIGMA x StripSigma above the beam again.
  /// RERISE 0 or less = off.
  Double_t TAIL_RETURN_NSIGMA;
  Double_t TAIL_RERISE_NSIGMA;

  /// Every strip from reac+1 to reac+
  /// POST_ABOVE_STRIPS (capped at the last strip) must sit more than
  /// POST_ABOVE_NSIGMA x StripSigma(s) above the beam, total[s] > 1 + nsigma
  /// x StripSigma(s).
  Double_t POST_ABOVE_NSIGMA;
  Int_t POST_ABOVE_STRIPS;

  /// The trace may not read below the beam mean before this strip: strips
  /// reac+1 .. S all read at or above the beam. An absolute strip, not a
  /// count from the reaction: a residue born at reac inherits the beam's
  /// remaining range from there, so it stops, and the trace crosses the
  /// beam, at nearly the same strip whatever reac was (about
  /// c x R_beam + (1 - c) x reac with c ~ 0.85; on 87Rb around strip 14).
  /// One value therefore serves every reaction strip, and the window it
  /// checks grows toward the early strips by itself. A light product or a
  /// fluctuation that leaves the beam envelope again within a strip or two
  /// fails it. Vacuous for reac >= S; a trace that never reads below the
  /// beam passes (the end-strip condition decides it). 0 = off. Part of the
  /// cut variation, shifted by one strip.
  Int_t POST_CROSS_MIN_STRIP;

  /// A residue stops gradually; the beam nucleus after a large-angle elastic
  /// scatter keeps its excess to the second-to-last strip and falls in one
  /// step at the last. Of the total fall from the trace's peak (over reac to
  /// last-1) to the last strip, the share taken by the last step,
  /// (total[last-1] - total[last]) / (total[peak] - total[last]), must stay
  /// under this fraction. Simulated 37Cl residues give 0.2 to 0.3, the elastic
  /// class about 0.9, a noise fake past the end-strip test about 0.7. 0 or
  /// less = off.
  Double_t TAIL_CLIFF_MAX_FRACTION;

  /// Cut-variation systematic, after the published 87Rb analysis: every
  /// condition that selects the beam or identifies the reaction is varied
  /// independently by the detector's resolution, and the change in the cross
  /// section is the systematic. With CUT_VARIATION on, the fill also counts,
  /// per strip, the tagged events and the beam denominator with each active
  /// condition shifted up and down by its step, one at a time: the
  /// sigma-scaled ones by NSIGMA_STEP sigma of the measured beam spread (3 =
  /// the resolution at 3 sigma, the paper's "±10%"), the event-level beam
  /// selection (gate level, pileup, noise, smoothness) among them, and
  /// TAIL_CLIFF_MAX_FRACTION by CLIFF_STEP. The cross section takes, per
  /// condition, the larger change of tagged / denominator under the two
  /// shifts and adds them in quadrature, with the gas-pressure uncertainty
  /// (GAS_PRESSURE_TORR_ERR). The counts travel with the cache.
  Bool_t CUT_VARIATION;
  Double_t CUT_VARIATION_NSIGMA_STEP;
  Double_t CUT_VARIATION_CLIFF_STEP;

  /// Normalize the partial dE sum on y axis by the event's own mean deposit
  /// over strips 1 .. reac-1, so per-event beam
  /// energy and gain jitter cancel instead of making events harder to see.
  Bool_t Y_RATIO_TO_UPSTREAM;

  /// Minimum jump at the reaction strip for a tag, in sigma of the measured
  /// beam spread of that strip (StripSumScatter::StripSigma).
  Double_t REAC_JUMP_NSIGMA;

  /// The end strip (17, or 16 with REQUIRE_STRIP_16_BELOW_BEAM or
  /// IGNORE_STRIP_17) must read below the beam there by this many sigma of
  /// its measured spread: total[end] < mean[end] - n x StripSigma(end). A
  /// stopped or stopping residue is low there, the beam and a still-flying
  /// elastic are not. Always applied; 0 is "below the beam mean", a negative
  /// value allows that much above it.
  Double_t END_STRIP_NSIGMA;

  /// Event-level cuts before any reaction is asked about, in sigma of each
  /// strip's measured beam spread (StripSumScatter::StripSigma). An event is
  /// pileup when PILEUP_MIN_STRIPS or more of strips 1-16 read at or above
  /// 1 + PILEUP_NSIGMA x StripSigma(s), and noise when NOISE_MIN_STRIPS or
  /// more read at or below 1 - NOISE_NSIGMA x StripSigma(s). Applied to the
  /// beam denominator too, so they cancel in the cross section. Both are
  /// inactive while the noise itself is being measured, when the clipped
  /// widths stand in for them. Part of the tagging, so a change refills.
  Double_t PILEUP_NSIGMA;
  Int_t PILEUP_MIN_STRIPS;
  Double_t NOISE_NSIGMA;
  Int_t NOISE_MIN_STRIPS;

  /// Redraw the regions for this reaction strip even if saved ones exist, and
  /// overwrite only that strip's entry.
  Bool_t REGION_CUT_REDRAW;

  /// compute-regions: the (a,n) region is the reaction component's
  /// AN_REGION_NSIGMA Mahalanobis ellipse from the bivariate Gaussian mixture
  /// fitted to each strip's scatter; the (a,a') region is the beam-like
  /// component's AA_REGION_NSIGMA ellipse. Saved like drawn cuts, so nothing
  /// downstream knows the difference.
  Double_t AN_REGION_NSIGMA;
  Double_t AA_REGION_NSIGMA;

  /// How compute-regions defines the (a,n) region. MIXTURE: a 2D
  /// Gaussian mixture per strip, beam plus reaction; the (a,n) region is the
  /// reaction component's AN_REGION_NSIGMA ellipse and the (a,a') region
  /// the beam's AA_REGION_NSIGMA ellipse, and the cross section takes the
  /// count the fit attributes to the reaction. ALL_TAGGED: use selection
  /// from strip-sum-scatter if event selection is strict enough to eliminate
  /// scattering events.build window and the cross section takes the tagged
  /// count with no enclosed-fraction correction and no region systematic.
  enum AnRegionMode { AN_REGION_MIXTURE, AN_REGION_ALL_TAGGED };
  AnRegionMode AN_REGION_MODE;

  /// How close to require strips to be to beam prior to the reaction,
  /// in sigma of each strip's measured beam spread
  /// (StripSumScatter::StripSigma). Part of the tagging, so a change refills.
  Bool_t REQUIRE_BEAM_UPSTREAM_OF_REAC;
  Double_t BEAM_UPSTREAM_NSIGMA;

  /// Use the segmentation of strips to reject an event when more than MAX
  /// strips in 1..COUNT_TO had BOTH ends fire. Read off raw ADC, so it is
  /// independent of IGNORE_SHORT_STRIPS. MAX < 0 disables it. COUNT_TO < 16
  /// restricts the count to strips upstream of the reaction.
  Int_t BOTH_MULT_MAX;
  Int_t BOTH_MULT_COUNT_TO;

  Int_t SEED_HALF_BINS;

  Int_t TRACES_PER_CLASS;

  Int_t X_LO;
  Int_t X_HI;

  /// The event-level beam gate: one strip, within GATE_NSIGMA of its fitted
  /// beam peak per gate group. GATE_NSIGMA at or below zero is no gate.
  Int_t GATE_STRIP;
  Double_t GATE_NSIGMA;
  /// The peak to fit, in beam units; zero lets the fit find it.
  Double_t GATE_CENTER;
  Double_t GATE_MIN;
  Double_t GATE_MAX;
  Int_t GATE_BINS;

  /// Display-only windows for the strip-sum scatters (a.u.): the histograms are
  /// built over the fixed ScatterBuildRange, so these only zoom the drawn plot
  /// via SetRangeUser (no rebuild); XBINS/YBINS do require a rebuild.
  Double_t X_DISPLAY_MIN;
  Double_t X_DISPLAY_MAX;
  Int_t XBINS;
  Double_t Y_DISPLAY_MIN;
  Double_t Y_DISPLAY_MAX;
  Int_t YBINS;

  /// Per-reaction-strip y-axis display windows, overriding Y_DISPLAY_MIN/MAX
  /// for individual strips (display-only, same as the X_DISPLAY_* windows).
  std::map<Int_t, std::pair<Double_t, Double_t>> Y_DISPLAY_RANGE;

  Long64_t SAMPLE_MAX_POINTS;

  Bool_t RERUN_SIM;
  Int_t CANDIDATE_REAC_STRIP;

  /// Leave out the Savitzky-Golay-smoothed copies of the trace overlays (the
  /// `_sg` figures and cluster-variable histograms), in both the
  /// strip-sum-scatter region overlays and compute-regions' all-tagged ones.
  Bool_t SKIP_SAVGOL_PLOTS;

  /// Also draw the per-region mean traces with RMS bands (the
  /// region_mean_traces_* figures) next to the trace overlays. Off by
  /// default: the overlay already carries the beam mean and its measured
  /// sigma band.
  Bool_t PLOT_REGION_MEAN_TRACES;
  /// Also draw the trace overlays in raw ADC (the *_adc figures) next to the
  /// calibrated ones, in strip-sum-scatter and in compute-regions' all-tagged
  /// overlays alike. Off by default.
  Bool_t PLOT_ADC_TRACES;
  Bool_t REQUIRE_STRIP_16_BELOW_BEAM;

  void SetDefaults();
};

/// TALYS calculation with name for plot and parameters to actually
/// run TALYS (see CrossSectionConfig::TALYS_MODELS).
/// One measured reaction channel; see CrossSectionConfig::CHANNELS.
/**
 * @brief One reaction channel for which the cross section is extracted.
 *
 * A dataset declares these in `CrossSectionConfig::CHANNELS`, one per reaction
 * being counted.
 */
struct CrossSectionChannel {
  /// Names the region cut (`region_<name>`) and the
  /// output figure. "an", "ap".
  TString name;
  /// The reaction as it appears in the plot title after the dataset name,
  /// e.g. "(#alpha, n)". Empty: derived from talys_exits -- a single "n" gives
  /// "(#alpha, n)", all-neutron exits give "(#alpha, xn)", anything else
  /// joins the exits with commas.
  TString label;
  /// Which of the residual channels TALYS wrote make up this channel's curve,
  /// as exit channels by name: what leaves the compound nucleus (beam +
  /// alpha), with an optional multiplicity digit, e.g. {"n", "2n"} for
  /// (a,xn), {"n"} for (a,n) alone, {"p"} for (a,p), "g" for radiative
  /// capture. The residue is the compound minus what left, so "n" on 37Cl is
  /// 40K and "p" is 40Ar. Read at plot time, so changing it never needs a
  /// TALYS rerun. A name that does not parse stops the run.
  std::vector<TString> talys_exits;
  /// Published values to compare against, if any, one row per point: the
  /// effective centre-of-mass energy [MeV] with its upward and downward
  /// uncertainties [MeV] (the strip's extent, asymmetric about that energy;
  /// zero when the table gives none), then the cross section [mb] with its
  /// total uncertainty [mb]. Empty when there is nothing to compare to.
  std::vector<std::vector<Double_t>> reference_xs;
  TString reference_label;
  /// Index into TALYS_MODELS of the model for the second figure
  /// (`cross_section_<name>_fit`): that model's sum and per-exit curves,
  /// unscaled, plus the sum with one scale per exit fitted to this work's
  /// points, its 3-sigma band, and a deviation panel. -1: no second figure.
  Int_t fit_model = -1;
};

/// @brief One TALYS model variant to compare against.
struct TalysModel {
  TString label;                 ///< Legend label for this model's curve.
  std::vector<TString> keywords; ///< TALYS keywords selecting the variant.
};

/**
 * @brief The fill gas the beam reacts in.
 *
 * Only helium so far. The per-gas numbers the cross section needs come from
 * TargetGasA() and TargetGasAtomsPerMolecule() rather than being repeated at
 * each use.
 */
enum TargetGas { kHELIUM /**< Helium fill. */ };

/// @brief Mass number of the target nucleus in a fill gas.
/// @param gas Fill gas.
Int_t TargetGasA(TargetGas gas);

/// @brief Atoms of the reacting species per molecule of fill gas.
/// @param gas Fill gas.
Double_t TargetGasAtomsPerMolecule(TargetGas gas);

/**
 * @brief Experimental details needed to properly calculate the cross section.
 *
 * The gas in which the beam reacts, and which reactions are being counted.
 */
struct CrossSectionConfig {
  TargetGas TARGET_GAS;
  /// At the pressure the gas was actually at rather than the nominal one.
  Double_t GAS_PRESSURE_TORR;
  /// Its uncertainty [Torr]; a relative systematic on every point (the
  /// cross section scales as 1 / pressure), added in quadrature with the
  /// cut variation. 0 = none.
  Double_t GAS_PRESSURE_TORR_ERR;

  /// Beam mass number, for the lab-to-centre-of-mass conversion; its Z and
  /// element symbol name it to a reaction code.
  Int_t BEAM_A;
  Int_t BEAM_Z;
  TString BEAM_ELEMENT;

  /// Hauser-Feshbach predictions from TALYS. talys-xs runs the code for an
  /// alpha on this beam nucleus over the reported strips' energies, once per
  /// model, and writes every residual-production channel to
  /// root_files/talys/talys_xs.root; cross-section sums the (a,xn) ones. Each
  /// model is a legend label and the input lines passed to TALYS verbatim
  /// after the projectile, target and energies -- "alphaomp 6" picks the
  /// alpha optical potential by TALYS's own index, for instance -- so any of
  /// its options is reachable without this config knowing them. The first
  /// model is the one drawn in full and the one whose shape sets each
  /// strip's effective energy; the others are drawn dashed and give the
  /// spread of that energy across shapes. Empty: no overlay, midpoint
  /// energies.
  std::vector<TalysModel> TALYS_MODELS;
  /// Input lines written into every model's TALYS input, after the
  /// projectile, target and energies and before the model's own keywords:
  /// the settings all models share (level density, accuracy thresholds,
  /// width fluctuations, pre-equilibrium, excitation-energy binning). Lines
  /// starting with `#` are comments and go through as such. TALYS takes the
  /// last value it reads for a keyword, so a model keyword repeating one of
  /// these overrides it. The defaults (Constants.cpp) follow the group's
  /// 14O(a,p) input; a dataset may replace the whole list.
  std::vector<TString> TALYS_COMMON_KEYWORDS;

  /// Simulated unreacted beam, read for the energy at each strip. Relative to
  /// the dataset's sim_root_files directory.
  TString BEAM_SIM_FILE;

  /// Reaction strips for which to calculate a cross-section.
  Int_t XS_STRIP_MIN;
  Int_t XS_STRIP_MAX;

  /// Statistical error on a strip's count. Below this count the 68.27 percent
  /// Feldman-Cousins interval on the Poisson mean (asymmetric, never below
  /// zero, an upper limit at zero counts); at or above it, root-N. There is
  /// no exact crossover: the two differ by 1/sqrt(N) with the discreteness
  /// ripple on top, about 40 percent of the bar at 4 counts, 15 at 25, 10 at
  /// 50, 5 at 80. 50 is where the step between the two is under a tenth of
  /// the bar and the asymmetry under a twentieth.
  Int_t FELDMAN_COUSINS_MAX_COUNT;

  /// Report each strip at its effective centre-of-mass energy (the energy
  /// below which half the strip's yield is produced, with the first TALYS
  /// model's cross section taken linear between the strip's entrance and
  /// exit values; the spread across models is an energy systematic), or,
  /// when off, at the strip's midpoint with the strip's extent as the only
  /// energy error. Off also means the TALYS shape never enters the measured
  /// points.
  Bool_t EFFECTIVE_ENERGY;

  /// Stamp the cross-section figures "PRELIMINARY" (large, faint, across the
  /// plot area) so a figure that leaves the analysis before its systematics
  /// are settled says so on its face.
  Bool_t PRELIMINARY;

  /// Epochs whose runs feed the cross-section chain (strip-sum-scatter,
  /// compute-regions, cross-section), by name. Those binaries run with no
  /// active epoch, so this is what their run list is built from. Empty: every
  /// enabled epoch. Ignored by a dataset that uses the flat RUN_NUMBERS list.
  /// An epoch at another pressure is kept for event building and calibration
  /// but left out here, since the gas density below is a single number.
  std::vector<TString> EPOCHS;

  /// The reaction channels measured on this dataset. Everything up to the
  /// tag is shared -- a jump is a jump whatever the residue -- and everything
  /// after it is per channel: the region cut (`region_<name>`), the count in
  /// it, the TALYS curve, the label, the published
  /// table, the figure (`cross_section_<name>`).
  std::vector<CrossSectionChannel> CHANNELS;

  void SetDefaults();
};

class DatasetConfig {
public:
  /// Data source
  Bool_t USE_SOLARIS_DATA;
  TString SOL_BASE_DIR;
  TString SOL_SPLIT_DIR;
  /// Length of the time chunks the SOLARIS `.sol` files are split into by
  /// the preprocess step, in seconds; each chunk becomes one subfile.
  /// Non-positive (e.g. -1) turns splitting off: the preprocess step does
  /// nothing and the pipeline reads each whole file from #SOL_BASE_DIR as one
  /// subfile, ignoring anything in #SOL_SPLIT_DIR. With epochs declared,
  /// MakeEpoch() copies this into `RunEpoch::split_chunk_seconds`, which an
  /// epoch may override; the tools read Constants::ActiveSplitChunkSeconds().
  Double_t SOL_SPLIT_CHUNK_SECONDS;
  Int_t SOL_N_SPLIT_WORKERS;

  TString COMPASS_BASE_DIR;

  /// Flat run list, used when EPOCHS is empty. When EPOCHS is populated the
  /// pipeline walks the epochs instead and this is ignored; binaries that run
  /// with no active epoch then take the runs of the epochs named in
  /// CROSS_SECTION_CONFIG.EPOCHS (all enabled ones when that is empty).
  std::vector<Int_t> RUN_NUMBERS;
  std::vector<RunEpoch> EPOCHS;
  Int_t N_CHUNKS;

  TString SIM_BEAM_FILE;

  /// Hardware layout
  Int_t N_BOARDS;
  Int_t N_CHANNELS;
  UShort_t TIMING_REF_BOARD;
  std::vector<UShort_t> TIMING_REF_BOARD_CHANNELS;

  /// Timing
  Double_t TIMING_MIN_ENERGY;
  Double_t TIMING_MAX_ENERGY;
  Double_t TIMING_OVERLAP_MARGIN_S;
  Double_t TIMING_THRESH_DT_US;
  Double_t TIMING_MAX_ABS_SHIFT_S;
  Double_t TIMING_SHIFT_COARSE_STEP_US;
  Double_t TIMING_SHIFT_FINE_STEP_US;
  Double_t TIMING_SHIFT_FINE_HALF_WIDTH_US;
  Int_t TIMING_SHIFT_MIN_NPTS;
  Int_t TIMING_SHIFT_MAX_SCAN_CANDIDATES;
  Bool_t TIMING_DO_BOARD_SYNC;
  Bool_t TIMING_DO_SORT;

  Bool_t REJECT_FLAGGED_EVENTS;

  /// Pole-zero pulse-history correction on the raw hits, before event
  /// building (see PulseHistory.hpp). The kernel is measured per subfile on
  /// its own beam-like events; MIN_EVENTS is the smallest sample that is
  /// trusted, BEAM_LO/HI the window (x the channel's beam peak) every long
  /// end must sit in for an event to count as beam, and APPLY_MAX_US how far
  /// back the correction looks. CORRECTION is the global switch; the groups
  /// are per epoch (see PULSE_HISTORY_GROUPS).
  Bool_t PULSE_HISTORY_CORRECTION;
  Long64_t PULSE_HISTORY_MIN_EVENTS;
  Double_t PULSE_HISTORY_BEAM_LO;
  Double_t PULSE_HISTORY_BEAM_HI;

  /// Looser window for the long ends of the chain whose kernel is being
  /// fit; the one above is used on strips that are not actively being fit
  /// in order to find actual beam-like events.
  Double_t PULSE_HISTORY_OWN_LO;
  Double_t PULSE_HISTORY_OWN_HI;

  // How far back to check and apply a correction. The trapezoid baseline is an
  // average over a
  /// few hundred microseconds, so pulses beyond 100 us still move it.
  Double_t PULSE_HISTORY_APPLY_MAX_US;

  /// Kernel bands in the previous pulse's amplitude, in units of the channel's
  /// beam peak: 1 is one kernel linear in the amplitude; N > 1 fits one kernel
  /// per band [0, 0.5), [0.5, 1.5), ..., [N-1.5, inf), i.e. single beam, twice,
  /// three times the beam pulse, so a nonlinear undershoot can be followed.
  Int_t PULSE_HISTORY_AMP_BINS;

  /// Which channel groups are corrected and with which kernel shape; all off
  /// by default, so a dataset that turns PULSE_HISTORY_CORRECTION on must
  /// also pick its groups, e.g.
  /// `gInstance.PULSE_HISTORY_GROUPS.long_left.enabled = kTRUE;` and
  /// `gInstance.PULSE_HISTORY_GROUPS.long_left.kernel = kPulseHistoryForm;`.
  /// With epochs declared, MakeEpoch() copies this set into
  /// `RunEpoch::pulse_history`, which an epoch may override group by group
  /// (a known tau for one era, a group off for another); the correction
  /// reads Constants::ActivePulseHistoryGroups().
  PulseHistoryGroups PULSE_HISTORY_GROUPS;

  /// Trapezoid filter settings of the DPP-PHA, in microseconds, which fix the
  /// two times in the pole-zero kernel form: the energy is read at
  /// `t_m = rise + peaking` after a pulse, and the baseline stays frozen for
  /// `t_f = 2 rise + flat` after a trigger. The form holds for gaps of at
  /// least `(2 rise + flat) - t_m`, when the previous trapezoid has passed
  /// under the read point. Defaults are board 66222's (2 / 3 / 1.5 us).
  Double_t PULSE_HISTORY_TRAP_RISE_US;
  Double_t PULSE_HISTORY_TRAP_FLAT_US;
  Double_t PULSE_HISTORY_PEAKING_US;

  /// Decode the long end only, dropping every short end. With epochs
  /// declared, MakeEpoch() copies this into `RunEpoch::ignore_short_strips`,
  /// which an epoch may override; the decode reads
  /// Constants::ActiveIgnoreShortStrips().
  Bool_t IGNORE_SHORT_STRIPS;
  /// Drop strip 0 / strip 17 from the analysis altogether: not required for
  /// a complete event, not in the all-strips cut, left off the trace figures.
  Bool_t IGNORE_STRIP_0;
  Bool_t IGNORE_STRIP_17;
  /// Whether a complete event must carry a strip 0 hit (EventBuilder's
  /// CheckEventComplete at build time, and the analysis' all-strips cut).
  /// Off, an event without one is kept; strip 0 is still stored, drawn and
  /// gated on where it fired. On 37Cl the strip 0 deposit sits at the DAQ
  /// threshold and fires in a third of the beam events, so requiring it
  /// costs two thirds of the statistics for a strip no tag condition reads.
  /// Moot when IGNORE_STRIP_0 is set.
  Bool_t REQUIRE_STRIP_0;

  Bool_t HAS_CATHODE;
  Bool_t HAS_GRID;
  Bool_t HAS_STRIP0;
  Bool_t HAS_STRIP17;

  Bool_t SKIP_EXISTING;

  /// Per-subfile output: the events summaries, timing, beam calibration and
  /// pulse-history figure folders of the pipeline, the per-group beam-gate
  /// folders of strip-sum-scatter, and the full per-file report in the log.
  /// Every subfile when SAVE_FULL_PLOTS; otherwise only the first
  /// PLOT_SAMPLE_FILES subfiles (chunks, or gate groups) of each epoch, in
  /// file order, enough to check every stage without a plots directory the
  /// size of the data; the other files write one summary line per stage.
  /// Off by default. Dataset-level figures are always drawn. Read through
  /// Constants::FileInSample() and Constants::SavePlots().
  Bool_t SAVE_FULL_PLOTS;
  Int_t PLOT_SAMPLE_FILES;

  /// Number of sample traces to save during event build and normed summary
  /// passes (0 = disabled). Saved as overlays in events_summary and
  /// events_summary_normed alongside the histograms.
  Int_t SAVE_SAMPLE_TRACES;

  Int_t MAX_FUSED_WORKERS;

  TString REFERENCE_CHANNEL;
  Double_t EVENT_TIME_WINDOW_US;
  /// Seed holdoff: a reference hit arriving within SEED_HOLDOFF_US after the
  /// current seed, when that seed is under SEED_HOLDOFF_MAX_RATIO times the
  /// new hit, is the pulse behind a pre-trigger -- the grid's slow rise trips
  /// the trigger early at a small energy, then the real pulse trips it again
  /// a few microseconds later -- and re-seeds the open event (the larger of
  /// the two energies is kept, SeedTs becomes the later hit's, the window is
  /// extended from it, the anodes already collected stay) instead of opening
  /// a second event that
  /// splits the anodes between the two and leaves both incomplete. On 37Cl
  /// run 97 the pre-triggers sit 3.5-8 us early and about 6 percent of grid
  /// seeds come out with no anodes at all. A pre-trigger is recognised by
  /// the open event holding no anode hits when the new reference hit
  /// arrives -- the particle's strips fire after the real grid trigger --
  /// with the size ratio only as a guard: on its own the ratio cannot tell a
  /// pre-trigger from a first particle read low by its follower, which lands
  /// at the same ratios (run 97: no gap between the pre-trigger pile at
  /// 0.0-0.2 and the follower-deficit floor from 0.3). 0 = off; the build
  /// summary counts the candidates either way, split by anodes pending, so
  /// the setting can be judged before it is on.
  Double_t SEED_HOLDOFF_US;
  Double_t SEED_HOLDOFF_MAX_RATIO;
  DedupStrategy DEDUP_STRATEGY;

  Bool_t USE_GPU_ACCELERATION;
  Int_t MAX_GPU_CONCURRENT_SORTS;

  /// n-sigma of the per-strip beam gate in the beam calibration: strip s is
  /// gated by the ellipse on the (strip s-1, strip s) raw totals, which is what
  /// defines that strip's beam sample. One level, in the fitted sigma: the two
  /// axes are the same kind of quantity and the ridge between them is tilted,
  /// so the correlation-aware contour is the gate and separate x and y levels
  /// would mean nothing.
  Double_t BEAM_GATE_NSIGMA;

  Double_t STRIP_DE_MIN_NORMED;
  Double_t STRIP_DE_MAX_NORMED;

  StripSumScatterConfig STRIP_SUM_SCATTER_CONFIG;
  CrossSectionConfig CROSS_SECTION_CONFIG;

  Double_t STRIP_E_MIN_ADC;
  Double_t STRIP_E_MAX_ADC;

  Double_t GRID_MAX_ADC;

  Double_t REFERENCE_CHANNEL_MIN_ADC;
  Double_t REFERENCE_CHANNEL_MAX_ADC;

  Double_t STRIP0_MAX_ADC;
  Double_t STRIP17_MAX_ADC;
  Double_t CATHODE_MAX_ADC;

  Double_t LEFT_EVEN_MAX_ADC;
  Double_t LEFT_ODD_MAX_ADC;
  Double_t RIGHT_EVEN_MAX_ADC;
  Double_t RIGHT_ODD_MAX_ADC;

  std::map<std::pair<Int_t, Int_t>, TString> channelMap;
  std::map<std::pair<Int_t, Int_t>, TString> channelMap64;

  DatasetConfig();
};

namespace Constants {
/// @brief The active dataset's configuration, flat block.
/// Prefer the `Active*()` accessors, which respect the epoch layer.
extern const DatasetConfig &cfg;

/**
 * @brief Set the epoch the `Active*()` accessors read from.
 *
 * Pipeline sets this around each epoch's work. Every `Active*()` accessor reads
 * it when set and falls back to the flat configuration block otherwise, so a
 * dataset declaring no epochs behaves exactly as it did before epochs existed.
 *
 * @param epoch Epoch to activate, or null to fall back to the flat block.
 */
void SetActiveEpoch(const RunEpoch *epoch);

/// @brief The currently active epoch, or null when none is set.
const RunEpoch *GetActiveEpoch();
/**
 * @brief The epoch owning a run number.
 * @param run Run number.
 * @return The epoch, or null when no untagged epoch declares it.
 * @note A tagged epoch reuses another era's run numbers and is addressed by
 *       tag, never by number alone — see RunEpoch::file_tag.
 */
const RunEpoch *EpochForRun(Int_t run);
/**
 * @brief Build an epoch prefilled from the flat configuration block.
 *
 * Carries the hardware layout, timing, reference channel, dedup strategy, ADC
 * caps and channel map, so a dataset states its detector settings once and its
 * epochs as one line each. Override fields on the result for an epoch that
 * genuinely differs.
 *
 * @param name Epoch name.
 * @param runs Run numbers belonging to it.
 * @return The prefilled epoch.
 *
 * @warning Call this **after** the flat block is set; it copies whatever is
 *          there at the time.
 */
RunEpoch MakeEpoch(const TString &name, const std::vector<Int_t> &runs);
/// @brief Consecutive run numbers, for MakeEpoch().
/// @param first First run, inclusive.
/// @param last  Last run, inclusive.
std::vector<Int_t> RunRange(Int_t first, Int_t last);
/// @brief Output-name prefix of the active epoch.
/// @return The tag, or an empty string when none is set — leaving existing
///         filenames unchanged.
const TString &ActiveFileTag();

/// @brief Boards in the active epoch's setup.
Int_t ActiveNBoards();
/// @brief Channels per board.
Int_t ActiveNChannels();
/// @brief Board the others are timing-aligned against.
UShort_t ActiveTimingRefBoard();
/// @brief Reference channel per board, for the alignment.
const std::vector<UShort_t> &ActiveTimingRefBoardChannels();
/// @brief Whether to run the multi-board timing alignment.
Bool_t ActiveDoBoardSync();
/// @brief Whether to time-sort hits before event building.
Bool_t ActiveDoSort();
/// @brief Whether this era instrumented the cathode.
Bool_t ActiveHasCathode();
/// @brief Whether this era's data is SOLARIS rather than CoMPASS.
Bool_t ActiveUseSolarisData();
/// @brief Coincidence window for event building, in microseconds.
Double_t ActiveEventTimeWindowUs();
/// @brief Seed holdoff in microseconds; 0 = off. See SEED_HOLDOFF_US.
Double_t ActiveSeedHoldoffUs();
/// @brief The size, as a fraction of the new hit, under which the current
///        seed counts as a pre-trigger. See SEED_HOLDOFF_MAX_RATIO.
Double_t ActiveSeedHoldoffMaxRatio();
/// @brief Channel whose hits seed events.
const TString &ActiveReferenceChannel();
/// @brief Lower energy gate on the seed channel, in ADC.
Double_t ActiveReferenceChannelMinAdc();
/// @brief Upper energy gate on the seed channel, in ADC.
Double_t ActiveReferenceChannelMaxAdc();
/// @brief How repeated hits on one slot are resolved.
DedupStrategy ActiveDedupStrategy();
/// @brief Lower bound of the per-strip energy range, in ADC.
Double_t ActiveStripEMinAdc();
/// @brief Upper bound of the per-strip energy range, in ADC.
Double_t ActiveStripEMaxAdc();
/// @brief Cathode full scale, in ADC.
Double_t ActiveCathodeMaxAdc();
/// @brief Grid full scale, in ADC.
Double_t ActiveGridMaxAdc();
/// @brief Strip 0 full scale, in ADC.
Double_t ActiveStrip0MaxAdc();
/// @brief Strip 17 full scale, in ADC.
Double_t ActiveStrip17MaxAdc();
/// @brief Left-side ceiling for even strips, in ADC.
Double_t ActiveLeftEvenMaxAdc();
/// @brief Left-side ceiling for odd strips, in ADC.
Double_t ActiveLeftOddMaxAdc();
/// @brief Right-side ceiling for even strips, in ADC.
Double_t ActiveRightEvenMaxAdc();
/// @brief Right-side ceiling for odd strips, in ADC.
Double_t ActiveRightOddMaxAdc();

/// @brief Runs of the active epoch, or the flat `RUN_NUMBERS` when none is set.
const std::vector<Int_t> &ActiveRunNumbers();

/// @brief Cap on subfiles processed per run for the active epoch.
/// @return The cap, or `-1` for all — which is also the answer with no epoch
/// set.
Int_t ActiveMaxFiles();

/// @brief Chunk length the active epoch's `.sol` files are split into, in
///        seconds; non-positive means whole files. The flat
///        `SOL_SPLIT_CHUNK_SECONDS` when no epoch is set.
Double_t ActiveSplitChunkSeconds();
/// @brief The pulse-history groups in force: the active epoch's set, or the
///        flat `PULSE_HISTORY_GROUPS` when no epoch is set.
const PulseHistoryGroups &ActivePulseHistoryGroups();
/// @brief Whether the decode keeps the long end only: the active epoch's
///        setting; with none active, the setting of the one epoch the
///        cross-section chain works on (AnalysisEpoch()), else the flat
///        `IGNORE_SHORT_STRIPS`.
Bool_t ActiveIgnoreShortStrips();
/// @brief Make AnalysisEpoch() the active epoch, so a binary of the
///        cross-section chain resolves the file tag, the channel map, the
///        ADC scales and the decode flags of the one epoch it works on. A
///        no-op when the chain spans several epochs or none are declared,
///        in which case the flat block serves as before. Called at the top
///        of strip-sum-scatter, compute-regions and cross-section.
void ActivateAnalysisEpoch();

/// @brief Whether the file this thread is processing is in the per-file
///        sample: every file with `SAVE_FULL_PLOTS`, otherwise the ones the
///        tool marked with SetPlotsThisFile(), the first `PLOT_SAMPLE_FILES`
///        of each epoch. A sample file draws its figures and writes its full
///        report to the log; the rest write one summary line per stage.
Bool_t FileInSample();
/// @brief FileInSample(), by the name the plot call sites use.
Bool_t SavePlots();
/// @brief `std::cout` for a file in the sample, a stream that discards
///        otherwise: the per-file detail lines (per strip, per channel) go
///        through this, so a file outside the sample keeps only its summary
///        lines. DetailErr() is the same for `std::cerr`.
std::ostream &Detail();
std::ostream &DetailErr();
/// @brief Whether the k-th file (0-based, in file order) of an epoch's list
///        is in the plot sample.
Bool_t InPlotSample(Int_t k);
/// @brief Mark the file this thread is about to process as in (or out of)
///        the plot sample. Per thread: a tool's worker calls it before each
///        file; the main thread starts out of the sample.
void SetPlotsThisFile(Bool_t on);

/**
 * @brief The channel map in force.
 * @return The active epoch's map when one is set, else the 64-channel map if
 *         populated, else the base map.
 * @warning Aborts the process when none is configured. A run with no channel
 *          map cannot be interpreted at all, so continuing would produce
 *          confidently wrong output.
 */
const std::map<std::pair<Int_t, Int_t>, TString> &ActiveChannelMap();
} // namespace Constants

#endif
