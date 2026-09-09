#include "Constants.hpp"
#include <DedupStrategy.hpp>
#include <RtypesCore.h>

namespace Constants {

static DatasetConfig gInstance;

const DatasetConfig &cfg = gInstance;

void InitDatasetConfig() {
  gInstance.USE_SOLARIS_DATA = kTRUE;
  gInstance.SOL_BASE_DIR = "/labdata/MUSIC/New-37Cl/data_raw/";
  gInstance.SOL_SPLIT_DIR = "/labdata/MUSIC/New-37Cl/data_split/";
  gInstance.SOL_N_SPLIT_WORKERS = 32;
  gInstance.SOL_SPLIT_CHUNK_SECONDS = 300;
  gInstance.COMPASS_BASE_DIR = "/labdata/MUSIC/37Cl/";

  gInstance.N_CHUNKS = -1;
  gInstance.SIM_BEAM_FILE = "traces_37Cl_beam.root";

  gInstance.N_BOARDS = 1;
  gInstance.N_CHANNELS = 64;

  gInstance.SAVE_PLOTS = kFALSE;
  gInstance.SKIP_EXISTING = kTRUE;
  gInstance.SAVE_SAMPLE_TRACES = 10;

  gInstance.HAS_CATHODE = kFALSE;
  gInstance.REFERENCE_CHANNEL = "Grid";
  gInstance.REFERENCE_CHANNEL_MIN_ADC = 1000;
  gInstance.REFERENCE_CHANNEL_MAX_ADC = 3000;
  gInstance.EVENT_TIME_WINDOW_US = 5.0;
  gInstance.DEDUP_STRATEGY = kLARGEST_ENERGY;

  // The long-end channels undershoot after every previous pulse (pole-zero
  // mismatch, ~12 us recovery); at 46 kHz that is most of the even-strip
  // width. Measured and removed per subfile before event building.
  gInstance.PULSE_HISTORY_CORRECTION = kTRUE;

  gInstance.BEAM_GATE_NSIGMA_X = 3;
  gInstance.BEAM_GATE_NSIGMA_Y = 3;
  gInstance.IGNORE_STRIP_0 = kTRUE;
  gInstance.IGNORE_STRIP_17 = kTRUE;
  gInstance.MAX_FUSED_WORKERS = 16;

  gInstance.STRIP_SUM_SCATTER_CONFIG.MAX_STRIP_SUM_WORKERS = 8;

  gInstance.STRIP_SUM_SCATTER_CONFIG.RERUN_SIM = kFALSE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.POST_TRIGGER_SUM_STRIPS = 6;
  gInstance.STRIP_SUM_SCATTER_CONFIG.POST_WINDOW_LAST_STRIP = 14;
  gInstance.STRIP_SUM_SCATTER_CONFIG.CANDIDATE_REAC_STRIP = 5;

  gInstance.STRIP_SUM_SCATTER_CONFIG.X_LO = 1;
  gInstance.STRIP_SUM_SCATTER_CONFIG.X_HI = 16;

  gInstance.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN = 2;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX = 8;
  gInstance.STRIP_SUM_SCATTER_CONFIG.X_DISPLAY_MIN = 13;
  gInstance.STRIP_SUM_SCATTER_CONFIG.X_DISPLAY_MAX = 19;
  gInstance.STRIP_SUM_SCATTER_CONFIG.Y_DISPLAY_MIN = 4;
  gInstance.STRIP_SUM_SCATTER_CONFIG.Y_DISPLAY_MAX = 9;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REGION_CUT_REDRAW = kTRUE;
  // The (a,n) cloud sits above the ridge and to the LEFT of the beam in
  // total energy (the neutron leaves with some), spread over ~2 a.u. in x:
  // not an island a Gaussian component can hold. Region = the band above
  // the ridge; the beam tail is heavier than Gaussian out to ~5 sigma.
  gInstance.STRIP_SUM_SCATTER_CONFIG.AN_REGION_MODE =
      StripSumScatterConfig::AN_REGION_RIDGE_BAND;
  gInstance.STRIP_SUM_SCATTER_CONFIG.AN_RIDGE_NSIGMA_LO = 5.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.AN_RIDGE_NSIGMA_HI = 15.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.SKIP_SAVGOL_PLOTS = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REJECT_NOISE = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REJECT_PILEUP = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.NOISE_THRESHOLD = 0.5;
  gInstance.STRIP_SUM_SCATTER_CONFIG.NOISE_MIN_STRIPS = 2;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REQUIRE_BEAM_UPSTREAM_OF_REAC = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REQUIRE_STRIP_16_BELOW_BEAM = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.END_STRIP_MAX = 1.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.PILEUP_THRESHOLD = 2;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REQUIRE_SMOOTHNESS = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.PARITY_ASYM_MAX = 0.15;
  gInstance.STRIP_SUM_SCATTER_CONFIG.PLOT_PARITY_REJECTED_GRID = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.GATE_STRIP_X = 0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.GATE_STRIP_Y = 1;
  gInstance.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_X = 5.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_Y = 5.0;

  gInstance.CROSS_SECTION_CONFIG.TARGET_GAS = kHELIUM;
  gInstance.CROSS_SECTION_CONFIG.GAS_PRESSURE_TORR = 400.0;
  gInstance.CROSS_SECTION_CONFIG.BEAM_A = 37;
  gInstance.CROSS_SECTION_CONFIG.BEAM_Z = 17;
  gInstance.CROSS_SECTION_CONFIG.BEAM_ELEMENT = "Cl";
  gInstance.CROSS_SECTION_CONFIG.TALYS_MODELS = {
      {"TALYS HF, Avrigeanu", {"alphaomp 6"}},
      {"TALYS HF, McFadden-Satchler", {"alphaomp 2"}}};
  gInstance.CROSS_SECTION_CONFIG.BEAM_SIM_FILE = "traces_37Cl_beam.root";
  gInstance.CROSS_SECTION_CONFIG.XS_STRIP_MIN = 2;
  gInstance.CROSS_SECTION_CONFIG.XS_STRIP_MAX = 8;
  gInstance.CROSS_SECTION_CONFIG.EFFECTIVE_ENERGY = kFALSE;
  // Two channels on the same tagged sample: the 40K and 40Ar residues, each
  // its own region cut, efficiency, TALYS curve and figure. No published
  // table yet.
  gInstance.CROSS_SECTION_CONFIG.CHANNELS = {
      {"an", "(#alpha, n)", {"n"}, {}, ""},
      {"ap", "(#alpha, p)", {"p"}, {}, ""}};
  // 400 Torr only; the 450 Torr epoch stays out of the scatter and the
  // cross section.
  gInstance.CROSS_SECTION_CONFIG.EPOCHS = {"late"};

  gInstance.STRIP0_MAX_ADC = 1000;
  gInstance.STRIP17_MAX_ADC = 10000;
  gInstance.GRID_MAX_ADC = 10000;

  gInstance.RIGHT_ODD_MAX_ADC = 2500;
  gInstance.RIGHT_EVEN_MAX_ADC = 8000;
  gInstance.LEFT_ODD_MAX_ADC = 8000;
  gInstance.LEFT_EVEN_MAX_ADC = 2500;

  gInstance.STRIP_DE_MIN_NORMED = 0;
  gInstance.STRIP_DE_MAX_NORMED = 4;

  gInstance.STRIP_E_MAX_ADC = 4096.0;
  gInstance.TOTAL_E_MAX_ADC = 60000.0;

  gInstance.channelMap64 = {
      {{0, 0}, "R1"},    {{0, 1}, "L1"},       {{0, 4}, "R2"},
      {{0, 5}, "L2"},    {{0, 8}, "R3"},       {{0, 9}, "L3"},
      {{0, 12}, "R4"},   {{0, 13}, "L4"},      {{0, 16}, "R5"},
      {{0, 17}, "L5"},   {{0, 20}, "R6"},      {{0, 21}, "L6"},
      {{0, 24}, "R7"},   {{0, 25}, "L7"},      {{0, 28}, "R8"},
      {{0, 29}, "L8"},   {{0, 32}, "R9"},      {{0, 33}, "L9"},
      {{0, 36}, "R10"},  {{0, 37}, "L10"},     {{0, 40}, "R11"},
      {{0, 41}, "L11"},  {{0, 44}, "R12"},     {{0, 45}, "L12"},
      {{0, 48}, "R13"},  {{0, 49}, "L13"},     {{0, 52}, "R14"},
      {{0, 53}, "L14"},  {{0, 56}, "R15"},     {{0, 57}, "L15"},
      {{0, 58}, "Grid"}, {{0, 59}, "Cathode"}, {{0, 60}, "R16"},
      {{0, 61}, "L16"},  {{0, 62}, "Strip17"}, {{0, 63}, "Strip0"}};

  // Epochs carry the flat block above; only the name and the run list differ.
  // Run 82 is the same detector at 450 Torr: built and calibrated like the
  // rest, kept apart for the cross section (CROSS_SECTION_CONFIG.EPOCHS).
  RunEpoch late = MakeEpoch("late", RunRange(97, 137));
  // RunEpoch Pressure450Torr = MakeEpoch("Pressure450Torr", {82});

  // gInstance.EPOCHS.push_back(Pressure450Torr);
  gInstance.EPOCHS.push_back(late);
}

// Static initializer runs before main
struct InitGuard {
  InitGuard() { InitDatasetConfig(); }
};
static InitGuard gInit;

} // namespace Constants
