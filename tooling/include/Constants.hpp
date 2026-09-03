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

// Progress logging on per-hit / per-event loops. Compile-time so the branch
// leaves the hot loop entirely when it is off; define MUSIC_HOT_PATH_LOGGING=0
// on the compiler command line to disable. Logging outside hot loops (per-file
// and per-run summaries) is unconditional and not covered by this.
#ifndef MUSIC_HOT_PATH_LOGGING
#define MUSIC_HOT_PATH_LOGGING 1
#endif

struct StripSumScatterConfig {
  enum PureBeamGate { PURE_BEAM_GATE_S0_S1, PURE_BEAM_GATE_S1_S2 };
  PureBeamGate PURE_BEAM_GATE;

  // Strips summed onto the scatter y-axis after the trigger strip: y spans
  // reac+1 .. min(reac+POST_TRIGGER_SUM_STRIPS, POST_WINDOW_LAST_STRIP), so
  // the window shrinks at deep strips instead of running into the region
  // where the recoil slows and its excess turns into a deficit. ApJ 983:142
  // did this by hand (five strips at 3-8, two at 9-12, one at 13);
  // POST_WINDOW_STRIPS sets a strip's window length outright, to reproduce
  // such a table. Part of the built quantity, so a change re-projects the
  // cache rather than refilling.
  Int_t POST_TRIGGER_SUM_STRIPS;
  Int_t POST_WINDOW_LAST_STRIP;
  std::map<Int_t, Int_t> POST_WINDOW_STRIPS;
  // Cap on worker threads for the scatter fill.
  Int_t MAX_STRIP_SUM_WORKERS;

  Int_t REACTION_STRIP_MIN;
  Int_t REACTION_STRIP_MAX;

  Int_t REQUIRE_SMOOTHNESS_END_STRIP;
  Double_t REQUIRE_SMOOTHNESS_MAX_STEP;

  // Minimum jump at the reaction strip for a tag, in sigma of the measured
  // strip-to-strip beam noise (StripSumScatter::JumpSigma). Part of the
  // tagging, so a change refills.
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

  // Redraw the regions for this reaction strip even if saved ones exist, and
  // overwrite only that strip's entry. Without it the only way to redraw was
  // to delete the whole cut file, which discarded every other strip's work.
  Bool_t REGION_CUT_REDRAW;

  // compute-regions: the (a,n) region is the reaction component's
  // AN_REGION_NSIGMA Mahalanobis ellipse from the bivariate Gaussian mixture
  // fitted to each strip's scatter; the (a,a') region is the beam-like
  // component's AA_REGION_NSIGMA ellipse. Saved like drawn cuts, so nothing
  // downstream knows the difference.
  Double_t AN_REGION_NSIGMA;
  Double_t AA_REGION_NSIGMA;

  // Tolerance in sigma of each strip's measured beam spread
  // (StripSumScatter::StripSigma). Part of the tagging, so a change refills.
  Bool_t REQUIRE_BEAM_UPSTREAM_OF_REAC;
  Double_t BEAM_UPSTREAM_NSIGMA;

  // Both-ends multiplicity cut: reject an event when more than MAX strips in
  // 1..COUNT_TO had BOTH ends fire. Read off raw ADC, so it is independent of
  // IGNORE_SHORT_STRIPS. Charge sharing means a displaced track lights the
  // short end on every strip of one parity at once, and tagging on a summed
  // trace then enriches for those -- this removes them. MAX < 0 disables it.
  // COUNT_TO < 16 restricts the count to strips upstream of the reaction, where
  // the sim has no sharing at all, without penalising products downstream.
  Int_t BOTH_MULT_MAX;
  Int_t BOTH_MULT_COUNT_TO;

  Bool_t REJECT_OFFBEAM;
  Double_t OFFBEAM_DIST;
  Int_t OFFBEAM_MIN_STRIPS;

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

  // Display-only windows for the strip-sum scatters (a.u.): the histograms are
  // built over the fixed ScatterBuildRange, so these only zoom the drawn plot
  // via SetRangeUser (no rebuild); XBINS/YBINS do require a rebuild.
  Double_t X_DISPLAY_MIN;
  Double_t X_DISPLAY_MAX;
  Int_t XBINS;
  Double_t Y_DISPLAY_MIN;
  Double_t Y_DISPLAY_MAX;
  Int_t YBINS;

  // Per-reaction-strip y-axis display windows, overriding Y_DISPLAY_MIN/MAX for
  // individual strips (display-only, same as the X_DISPLAY_* windows).
  std::map<Int_t, std::pair<Double_t, Double_t>> Y_DISPLAY_RANGE;

  Long64_t SAMPLE_MAX_POINTS;

  Bool_t RERUN_SIM;
  Int_t CANDIDATE_REAC_STRIP;

  Bool_t REQUIRE_SMOOTHNESS;
  Bool_t REQUIRE_GATE_S3_S4;
  Bool_t REQUIRE_GATE_S5_S6;
  Bool_t SKIP_SAVGOL_PLOTS;
  Bool_t REQUIRE_STRIP_16_BELOW_BEAM;

  // Also render the selected region traces under the OTHER decode
  // (long-side-only <-> L+R sum) from the same calibration and the
  // same selected events, into a separate plot plus a text dump.
  Bool_t ALT_DECODE_REGION_TRACES;

  void SetDefaults();
};

// One TALYS calculation: what the plot calls it, and the input lines that
// select it (see CrossSectionConfig::TALYS_MODELS).
struct TalysModel {
  TString label;
  std::vector<TString> keywords;
};

// The fill gas the beam reacts in. Only helium so far; the per-gas numbers
// the cross section needs come from TargetGasA and TargetGasAtomsPerMolecule.
enum TargetGas { kHELIUM };
// Mass number of the target nucleus, and atoms of the reacting species per
// gas molecule.
Int_t TargetGasA(TargetGas gas);
Double_t TargetGasAtomsPerMolecule(TargetGas gas);

// Everything the cross section needs that is a property of the experiment
// rather than of the analysis: the gas the beam reacts in and which reaction
// is being counted. Nothing here is 87Rb-specific -- a second dataset supplies
// its own values and reuses the tool unchanged.
struct CrossSectionConfig {
  TargetGas TARGET_GAS;
  // At the pressure the gas was actually at rather than the nominal one.
  Double_t GAS_PRESSURE_TORR;

  // Beam mass number, for the lab-to-centre-of-mass conversion; its Z and
  // element symbol name it to a reaction code.
  Int_t BEAM_A;
  Int_t BEAM_Z;
  TString BEAM_ELEMENT;

  // Hauser-Feshbach predictions from TALYS. talys-xs runs the code for an
  // alpha on this beam nucleus over the reported strips' energies, once per
  // model, and writes every residual-production channel to
  // root_files/talys/talys_xs.root; cross-section sums the (a,xn) ones. Each
  // model is a legend label and the input lines passed to TALYS verbatim
  // after the projectile, target and energies -- "alphaomp 6" picks the
  // alpha optical potential by TALYS's own index, for instance -- so any of
  // its options is reachable without this config knowing them. The first
  // model is the one drawn in full and the one whose shape sets each
  // strip's effective energy; the others are drawn dashed and give the
  // spread of that energy across shapes. Empty: no overlay, midpoint
  // energies.
  std::vector<TalysModel> TALYS_MODELS;

  // Simulated unreacted beam, read for the energy at each strip. Relative to
  // the dataset's sim_root_files directory.
  TString BEAM_SIM_FILE;

  // Reaction strips to report a cross section for.
  Int_t XS_STRIP_MIN;
  Int_t XS_STRIP_MAX;

  // Published values to compare against, if any, one row per point: the
  // effective centre-of-mass energy [MeV] with its upward and downward
  // uncertainties [MeV] (the strip's extent, asymmetric about that energy;
  // zero when the table gives none), then the cross section [mb] with its
  // total uncertainty [mb]. Left empty when there is nothing to compare to.
  std::vector<std::vector<Double_t>> REFERENCE_XS;
  TString REFERENCE_LABEL;

  void SetDefaults();
};

class DatasetConfig {
public:
  // Data source
  Bool_t USE_SOLARIS_DATA;
  TString SOL_BASE_DIR;
  TString SOL_SPLIT_DIR;
  Double_t SOL_SPLIT_CHUNK_SECONDS;
  Int_t SOL_N_SPLIT_WORKERS;

  TString COMPASS_BASE_DIR;
  // Flat run list, used when EPOCHS is empty. When EPOCHS is populated the
  // pipeline walks the epochs instead and this is ignored.
  std::vector<Int_t> RUN_NUMBERS;
  std::vector<RunEpoch> EPOCHS;
  Int_t N_CHUNKS;

  TString SIM_BEAM_FILE;

  // Hardware layout
  Int_t N_BOARDS;
  Int_t N_CHANNELS;
  UShort_t TIMING_REF_BOARD;
  std::vector<UShort_t> TIMING_REF_BOARD_CHANNELS;

  // Timing
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

  // Number of sample traces to save during event build and normed summary
  // passes (0 = disabled). Saved as overlays in events_summary and
  // events_summary_normed alongside the histograms.
  Int_t SAVE_SAMPLE_TRACES;

  Int_t MAX_FUSED_WORKERS;

  TString REFERENCE_CHANNEL;
  Double_t EVENT_TIME_WINDOW_US;
  DedupStrategy DEDUP_STRATEGY;

  Bool_t USE_GPU_ACCELERATION;
  Int_t MAX_GPU_CONCURRENT_SORTS;

  Double_t STRIP_DE_OVERVIEW_MIN_NORMED;
  Double_t STRIP_DE_OVERVIEW_MAX_NORMED;
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
extern const DatasetConfig &cfg;

// Epoch selection. Pipeline sets the active epoch around each epoch's work;
// every Active*() below reads it when set and falls back to the flat cfg block
// otherwise, so a dataset that declares no epochs behaves exactly as before.
void SetActiveEpoch(const RunEpoch *epoch);
const RunEpoch *GetActiveEpoch();
// Epoch owning a run number, or null when no UNTAGGED epoch declares it. A
// tagged epoch reuses another era's run numbers and is addressed by tag, never
// by number alone.
const RunEpoch *EpochForRun(Int_t run);
// Output-name prefix of the active epoch ("" when none, so existing filenames
// are unchanged).
const TString &ActiveFileTag();

Int_t ActiveNBoards();
Int_t ActiveNChannels();
UShort_t ActiveTimingRefBoard();
const std::vector<UShort_t> &ActiveTimingRefBoardChannels();
Bool_t ActiveDoBoardSync();
Bool_t ActiveDoSort();
Bool_t ActiveHasCathode();
Bool_t ActiveUseSolarisData();
Double_t ActiveEventTimeWindowUs();
const TString &ActiveReferenceChannel();
Double_t ActiveReferenceChannelMinAdc();
Double_t ActiveReferenceChannelMaxAdc();
DedupStrategy ActiveDedupStrategy();
Double_t ActiveStripEMinAdc();
Double_t ActiveStripEMaxAdc();
Double_t ActiveCathodeMaxAdc();
Double_t ActiveGridMaxAdc();
Double_t ActiveStrip0MaxAdc();
Double_t ActiveStrip17MaxAdc();
Double_t ActiveLeftEvenMaxAdc();
Double_t ActiveLeftOddMaxAdc();
Double_t ActiveRightEvenMaxAdc();
Double_t ActiveRightOddMaxAdc();

// Runs of the active epoch, else the flat RUN_NUMBERS. ActiveMaxFiles caps
// subfiles per run for the active epoch (-1 = all, and -1 without an epoch).
const std::vector<Int_t> &ActiveRunNumbers();
Int_t ActiveMaxFiles();

// Returns the active epoch's channel map when an epoch is set, else
// channelMap64 if populated, else channelMap. Aborts if none is configured.
const std::map<std::pair<Int_t, Int_t>, TString> &ActiveChannelMap();
} // namespace Constants

#endif
