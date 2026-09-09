#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include "DedupStrategy.hpp"
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
 * @brief The dataset configuration, and how it is layered.
 *
 * Configuration is split in two. The struct definitions and the tooling-wide
 * defaults live here and in `Constants.cpp`; each dataset then overrides only
 * the fields it cares about in `analysis/<dataset>/config/Constants.cpp`, which
 * is compiled into that dataset's binaries.
 *
 * On top of that sits the epoch layer. A dataset spanning several acquisition
 * eras declares a RunEpoch per era, and the `Active*()` accessors read the
 * active epoch when one is set and the flat block otherwise. See RunEpoch for
 * why a half-populated epoch is an error rather than a fallback.
 */

/**
 * @def MUSIC_HOT_PATH_LOGGING
 * @brief Whether per-hit and per-event progress logging is compiled in.
 *
 * Compile-time so the branch leaves the hot loop entirely when off. Define it
 * to `0` on the compiler command line to disable.
 *
 * @note Logging outside hot loops — per-file and per-run summaries — is
 *       unconditional and not covered by this.
 */
#ifndef MUSIC_HOT_PATH_LOGGING
#define MUSIC_HOT_PATH_LOGGING 1
#endif

/**
 * @brief Everything governing the reaction search in the strip-sum scatters.
 *
 * Two kinds of knob live here, and the distinction decides what a change costs:
 * those that are part of the **built** quantity, where a change reprojects the
 * cache, and those that are part of the **tagging**, where a change refills it.
 * Each field says which it is. The cache fingerprints them all, so either
 * happens automatically.
 */
struct StripSumScatterConfig {
  enum PureBeamGate { PURE_BEAM_GATE_S0_S1, PURE_BEAM_GATE_S1_S2 };
  PureBeamGate PURE_BEAM_GATE;

  /// Strips summed onto the scatter y-axis after the trigger strip: y spans
  /// reac+1 .. min(reac+POST_TRIGGER_SUM_STRIPS, POST_WINDOW_LAST_STRIP), so
  /// the window shrinks at deep strips instead of running into the region
  /// where the recoil slows and its excess turns into a deficit. ApJ 983:142
  /// did this by hand (five strips at 3-8, two at 9-12, one at 13);
  /// POST_WINDOW_STRIPS sets a strip's window length outright, to reproduce
  /// such a table. Part of the built quantity, so a change re-projects the
  /// cache rather than refilling.
  Int_t POST_TRIGGER_SUM_STRIPS;
  Int_t POST_WINDOW_LAST_STRIP;
  std::map<Int_t, Int_t> POST_WINDOW_STRIPS;
  /// Cap on worker threads for the scatter fill.
  Int_t MAX_STRIP_SUM_WORKERS;

  Int_t REACTION_STRIP_MIN;
  Int_t REACTION_STRIP_MAX;

  Int_t REQUIRE_SMOOTHNESS_END_STRIP;
  Double_t REQUIRE_SMOOTHNESS_MAX_STEP;

  /// Minimum jump at the reaction strip for a tag, in sigma of the measured
  /// strip-to-strip beam noise (StripSumScatter::JumpSigma). Part of the
  /// tagging, so a change refills.
  Double_t REAC_JUMP_NSIGMA;
  Double_t REAC_JUMP_MAX;
  Double_t END_STRIP_MAX;

  Double_t PILEUP_THRESHOLD;
  Double_t NOISE_THRESHOLD;

  Bool_t REJECT_NOISE;
  Double_t NOISE_THRESH_PY;
  Int_t NOISE_MIN_STRIPS;

  Bool_t REJECT_PILEUP;
  Double_t PILEUP_THRESH_PY;
  Int_t PILEUP_MIN_STRIPS;

  /// Redraw the regions for this reaction strip even if saved ones exist, and
  /// overwrite only that strip's entry. Without it the only way to redraw was
  /// to delete the whole cut file, which discarded every other strip's work.
  Bool_t REGION_CUT_REDRAW;

  /// compute-regions: the (a,n) region is the reaction component's
  /// AN_REGION_NSIGMA Mahalanobis ellipse from the bivariate Gaussian mixture
  /// fitted to each strip's scatter; the (a,a') region is the beam-like
  /// component's AA_REGION_NSIGMA ellipse. Saved like drawn cuts, so nothing
  /// downstream knows the difference.
  Double_t AN_REGION_NSIGMA;
  Double_t AA_REGION_NSIGMA;
  /// How compute-regions defines the (a,n) region. MIXTURE: the ellipse of
  /// the mixture's reaction component, for a compact island beyond the beam
  /// (87Rb). RIDGE_BAND: everything between AN_RIDGE_NSIGMA_LO and _HI
  /// conditional sigma above the beam ridge, free in x across the window, for
  /// a reaction cloud that spreads along x and sits below the beam in total
  /// energy (37Cl, where the neutron carries energy out). The (a,a') region
  /// is the beam ellipse in both.
  enum AnRegionMode { AN_REGION_MIXTURE, AN_REGION_RIDGE_BAND };
  AnRegionMode AN_REGION_MODE;
  Double_t AN_RIDGE_NSIGMA_LO;
  Double_t AN_RIDGE_NSIGMA_HI;

  /// Tolerance in sigma of each strip's measured beam spread
  /// (StripSumScatter::StripSigma). Part of the tagging, so a change refills.
  Bool_t REQUIRE_BEAM_UPSTREAM_OF_REAC;
  Double_t BEAM_UPSTREAM_NSIGMA;

  /// Both-ends multiplicity cut: reject an event when more than MAX strips in
  /// 1..COUNT_TO had BOTH ends fire. Read off raw ADC, so it is independent of
  /// IGNORE_SHORT_STRIPS. Charge sharing means a displaced track lights the
  /// short end on every strip of one parity at once, and tagging on a summed
  /// trace then enriches for those -- this removes them. MAX < 0 disables it.
  /// COUNT_TO < 16 restricts the count to strips upstream of the reaction,
  /// where the sim has no sharing at all, without penalising products
  /// downstream.
  Int_t BOTH_MULT_MAX;
  Int_t BOTH_MULT_COUNT_TO;

  Bool_t REJECT_OFFBEAM;
  Double_t OFFBEAM_DIST;
  Int_t OFFBEAM_MIN_STRIPS;
  /// Reject events whose even strips and odd strips disagree by more than this
  /// fraction, |mean(even 2..16) / mean(odd 1..16) - 1|, before the tag. A beam
  /// particle or a residue reads the same on both parities to within the
  /// noise; an event that arrives on the pole-zero undershoot of the previous
  /// pulse reads low on every channel of the group with the wrong pole-zero
  /// (the R channels in the SOLARIS 37Cl runs), a sawtooth at 0.65 on even
  /// strips and 1.0 on odd. Applied to the beam denominator too, so it
  /// cancels. 0 or less: off.
  Double_t PARITY_ASYM_MAX;

  /// Diagnostic only: when set, run an extra pass over the events (gated by
  /// PARITY_ASYM_MAX > 0) that fills a histogram of the Grid `#DeltaE` of
  /// events that pass the cheap pre-tag cuts (all strips fired, pileup, noise,
  /// offbeam) and are then rejected by the parity cut, and save BOTH a decoded
  /// a.u. view (grid_adc / 16384, [0,1]) and a raw ADC view ([0,
  /// GRID_MAX_ADC]), each with a log-y axis. Purely visual; does not change
  /// what is tagged and is NOT part of the cache fingerprint (the reservoir
  /// keeps only tagged + beam events, so it cannot be rebuilt from cache
  /// anyway). Requires the Grid branch to be enabled. Off by default.
  Bool_t PLOT_PARITY_REJECTED_GRID;

  Double_t TRIGGER_NSIGMA;
  Double_t TRIGGER_CFD_FRAC;
  Int_t PLATEAU_POST;
  Int_t CLUSTER_SMOOTH_WINDOW;
  Int_t SEED_HALF_BINS;

  Int_t SAVGOL_HALF;

  Int_t TRACES_PER_CLASS;

  Int_t X_LO;
  Int_t X_HI;

  Int_t GATE_STRIP_X;
  Int_t GATE_STRIP_Y;
  Double_t GATE_NSIGMA_X;
  Double_t GATE_NSIGMA_Y;
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

  Bool_t REQUIRE_SMOOTHNESS;
  Bool_t REQUIRE_GATE_S3_S4;
  Bool_t REQUIRE_GATE_S5_S6;
  Bool_t SKIP_SAVGOL_PLOTS;
  Bool_t REQUIRE_STRIP_16_BELOW_BEAM;

  /// Also render the selected region traces under the OTHER decode
  /// (long-side-only <-> L+R sum) from the same calibration and the
  /// same selected events, into a separate plot plus a text dump.
  Bool_t ALT_DECODE_REGION_TRACES;

  void SetDefaults();
};

// One TALYS calculation: what the plot calls it, and the input lines that
// select it (see CrossSectionConfig::TALYS_MODELS).
// One measured reaction channel; see CrossSectionConfig::CHANNELS.
/**
 * @brief One reaction channel the cross section is extracted for.
 *
 * A dataset declares these in `CrossSectionConfig::CHANNELS`, one per reaction
 * being counted.
 */
struct CrossSectionChannel {
  /// Names the region cut (region_<name>), the tag-efficiency records and the
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
 * @brief What the cross section needs about the *experiment* rather than the
 * analysis.
 *
 * The gas the beam reacts in, and which reactions are being counted. Nothing
 * here is specific to any one dataset — a second experiment supplies its own
 * values and reuses the tool unchanged.
 */
struct CrossSectionConfig {
  TargetGas TARGET_GAS;
  /// At the pressure the gas was actually at rather than the nominal one.
  Double_t GAS_PRESSURE_TORR;

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

  /// Simulated unreacted beam, read for the energy at each strip. Relative to
  /// the dataset's sim_root_files directory.
  TString BEAM_SIM_FILE;

  /// Reaction strips to report a cross section for.
  Int_t XS_STRIP_MIN;
  Int_t XS_STRIP_MAX;

  /// Report each strip at its effective centre-of-mass energy (Szegedi et al.
  /// 2021: the energy at which the first TALYS model's cross section equals
  /// its average over the strip, with the spread across models as an energy
  /// systematic), or, when off, at the strip's midpoint with the strip's
  /// extent as the only energy error. Off also means the TALYS shape never
  /// enters the measured points.
  Bool_t EFFECTIVE_ENERGY;

  /// Epochs whose runs feed the cross-section chain (strip-sum-scatter,
  /// compute-regions, cross-section), by name. Those binaries run with no
  /// active epoch, so this is what their run list is built from. Empty: every
  /// enabled epoch. Ignored by a dataset that uses the flat RUN_NUMBERS list.
  /// An epoch at another pressure is kept for event building and calibration
  /// but left out here, since the gas density below is a single number.
  std::vector<TString> EPOCHS;

  /// The reaction channels measured on this dataset. Everything up to the
  /// tag is shared -- a jump is a jump whatever the residue -- and everything
  /// after it is per channel: the region cut (region_<name>), the count in it,
  /// the tag-efficiency record, the TALYS curve, the label, the published
  /// table, the figure (cross_section_<name>).
  std::vector<CrossSectionChannel> CHANNELS;

  void SetDefaults();
};

class DatasetConfig {
public:
  /// Data source
  Bool_t USE_SOLARIS_DATA;
  TString SOL_BASE_DIR;
  TString SOL_SPLIT_DIR;
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
  /// back the correction looks (the undershoot is gone by 40 us; beyond ~100
  /// us the fitted bins are degenerate with the intercept).
  Bool_t PULSE_HISTORY_CORRECTION;
  Long64_t PULSE_HISTORY_MIN_EVENTS;
  Double_t PULSE_HISTORY_BEAM_LO;
  Double_t PULSE_HISTORY_BEAM_HI;
  /// Looser window for the long ends of the chain whose kernel is being
  /// fitted; the tight one above applies to the other chain. Tight on both
  /// would cut off the large undershoots the kernel exists to describe.
  Double_t PULSE_HISTORY_OWN_LO;
  Double_t PULSE_HISTORY_OWN_HI;
  Double_t PULSE_HISTORY_APPLY_MAX_US;
  /// Kernel bands in the previous pulse's amplitude, in units of the channel's
  /// beam peak: 1 is one kernel linear in the amplitude; N > 1 fits one kernel
  /// per band [0, 0.5), [0.5, 1.5), ..., [N-1.5, inf), i.e. single beam, twice,
  /// three times the beam pulse, so a nonlinear undershoot can be followed.
  Int_t PULSE_HISTORY_AMP_BINS;

  Bool_t IGNORE_SHORT_STRIPS;
  Bool_t IGNORE_STRIP_0;
  Bool_t IGNORE_STRIP_17;

  Bool_t HAS_CATHODE;
  Bool_t HAS_GRID;
  Bool_t HAS_STRIP0;
  Bool_t HAS_STRIP17;

  Bool_t SKIP_EXISTING;
  Bool_t SAVE_PLOTS;

  Bool_t SKIP_CALIBRATION;

  /// Number of sample traces to save during event build and normed summary
  /// passes (0 = disabled). Saved as overlays in events_summary and
  /// events_summary_normed alongside the histograms.
  Int_t SAVE_SAMPLE_TRACES;

  Int_t MAX_FUSED_WORKERS;

  TString REFERENCE_CHANNEL;
  Double_t EVENT_TIME_WINDOW_US;
  DedupStrategy DEDUP_STRATEGY;

  Bool_t USE_GPU_ACCELERATION;
  Int_t MAX_GPU_CONCURRENT_SORTS;

  Double_t STRIP_DE_OVERVIEW_MIN_NORMED;
  Double_t STRIP_DE_OVERVIEW_MAX_NORMED;
  /// n-sigma of the per-strip beam gate in the beam calibration: strip s is
  /// gated by the ellipse on the (strip s-1, strip s) raw totals, which is what
  /// defines that strip's beam sample.
  Double_t BEAM_GATE_NSIGMA_X;
  Double_t BEAM_GATE_NSIGMA_Y;

  Double_t STRIP_DE_MIN_NORMED;
  Double_t STRIP_DE_MAX_NORMED;
  Double_t CATHODE_E_MAX_NORMED;
  Double_t TOTAL_E_MIN_NORMED;
  Double_t TOTAL_E_MAX_NORMED;

  StripSumScatterConfig STRIP_SUM_SCATTER_CONFIG;
  CrossSectionConfig CROSS_SECTION_CONFIG;

  Double_t STRIP_E_MIN_ADC;
  Double_t STRIP_E_MAX_ADC;
  Double_t TOTAL_E_MIN_ADC;
  Double_t TOTAL_E_MAX_ADC;

  Double_t GRID_MIN_ADC;
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
  std::map<std::pair<Int_t, Int_t>, Long64_t> ttfOffsetPs;

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
