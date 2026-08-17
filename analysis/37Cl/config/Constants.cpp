#include "Constants.hpp"
#include <RtypesCore.h>

namespace Constants {

static DatasetConfig gInstance;

const DatasetConfig &cfg = gInstance;

void InitDatasetConfig() {
  gInstance.USE_SOLARIS_DATA = kTRUE;
  gInstance.SOL_BASE_DIR = "/labdata/MUSIC/New-37Cl/data_raw/";
  gInstance.SOL_SPLIT_DIR = "/labdata/MUSIC/New-37Cl/data_split/";
  gInstance.SOL_N_SPLIT_WORKERS = 32;
  gInstance.SOL_SPLIT_CHUNK_SECONDS = 60;
  gInstance.COMPASS_BASE_DIR = "/labdata/MUSIC/37Cl/";
  std::vector<Int_t> RUNS;
  for (Int_t i = 30; i < 53; i++) {
    RUNS.push_back(i);
  };
  for (Int_t i = 84; i < 138; i++) {
    RUNS.push_back(i);
  };
  gInstance.RUN_NUMBERS = RUNS;
  gInstance.N_CHUNKS = -1;
  gInstance.SIM_BEAM_FILE = "traces_37Cl_beam.root";

  gInstance.N_BOARDS = 1;
  gInstance.N_CHANNELS = 64;
  gInstance.TIMING_REF_BOARD = 0;
  gInstance.TIMING_REF_BOARD_CHANNELS = {8};

  gInstance.SAVE_PLOTS = kFALSE;
  gInstance.SAVE_SAMPLE_TRACES = 10;

  gInstance.HAS_CATHODE = kFALSE;
  gInstance.REFERENCE_CHANNEL = "Grid";
  gInstance.REFERENCE_CHANNEL_MIN_ADC = 500;
  gInstance.REFERENCE_CHANNEL_MAX_ADC = 4500;
  gInstance.EVENT_TIME_WINDOW_US = 8.0;
  gInstance.DEDUP_STRATEGY = kLARGEST_ENERGY;

  gInstance.TIMING_DO_BOARD_SYNC = kFALSE;
  gInstance.TIMING_DO_SORT = kFALSE;

  // Automatic pileup removal: any event where an anode channel fired twice
  // inside the 8 us window is two ions, not one — drop it (the per-channel
  // dedup would otherwise keep the larger hit and inflate one strip).
  gInstance.REJECT_MULTI_HIT_EVENTS = kTRUE;

  gInstance.IGNORE_STRIP_0 = kTRUE;
  gInstance.IGNORE_STRIP_17 = kTRUE;
  gInstance.MAX_FUSED_WORKERS = 32;

  gInstance.STRIP_SUM_SCATTER_CONFIG.MAX_STRIP_SUM_WORKERS = 8;
  gInstance.STRIP_SUM_SCATTER_CONFIG.PURE_BEAM_GATE =
      StripSumScatterConfig::PURE_BEAM_GATE_S1_S2;
  gInstance.STRIP_SUM_SCATTER_CONFIG.RERUN_SIM = kFALSE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.CANDIDATE_REAC_STRIP = 4;
  gInstance.STRIP_SUM_SCATTER_CONFIG.XMIN = 7;
  gInstance.STRIP_SUM_SCATTER_CONFIG.XMAX = 16;
  gInstance.STRIP_SUM_SCATTER_CONFIG.X_HI = 12;
  gInstance.STRIP_SUM_SCATTER_CONFIG.Y_RANGE = {
      {3, {0, 17.5}},  {4, {0, 17.5}}, {5, {0, 17.5}}, {6, {0, 17.5}},
      {7, {0, 17.5}},  {8, {0, 17.5}}, {9, {0, 17.5}}, {10, {0, 17.5}},
      {11, {0, 17.5}}, {12, {0, 10}},  {13, {0, 10}},  {14, {0, 10}},
      {15, {0, 10}}};
  gInstance.STRIP_SUM_SCATTER_CONFIG.SKIP_SAVGOL_PLOTS = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.SKIP_CLUSTER_HISTS = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REJECT_NOISE = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REJECT_PILEUP = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.NOISE_THRESHOLD = 0.85;
  gInstance.STRIP_SUM_SCATTER_CONFIG.NOISE_MIN_STRIPS = 5;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REQUIRE_STRIP_16_BELOW_BEAM = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.END_STRIP_MAX = 1.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.PILEUP_THRESHOLD = 1.4;
  gInstance.STRIP_SUM_SCATTER_CONFIG.PILEUP_MIN_STRIPS = 2;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REJECT_HIGH_STRIP = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.HIGH_STRIP_THRESHOLD = 2.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REQUIRE_SMOOTHNESS = kFALSE;

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
}

// Static initializer runs before main
struct InitGuard {
  InitGuard() { InitDatasetConfig(); }
};
static InitGuard gInit;

} // namespace Constants
