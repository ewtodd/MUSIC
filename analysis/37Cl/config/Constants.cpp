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

  gInstance.IGNORE_STRIP_0 = kTRUE;
  gInstance.IGNORE_STRIP_17 = kTRUE;
  gInstance.MAX_FUSED_WORKERS = 32;

  gInstance.STRIP_SUM_SCATTER_CONFIG.PURE_BEAM_GATE =
      StripSumScatterConfig::PURE_BEAM_GATE_S1_S2;
  gInstance.STRIP_SUM_SCATTER_CONFIG.RERUN_SIM = kFALSE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.CANDIDATE_REAC_STRIP = 5;
  gInstance.STRIP_SUM_SCATTER_CONFIG.XMIN = 11;
  gInstance.STRIP_SUM_SCATTER_CONFIG.XMAX = 19;
  gInstance.STRIP_SUM_SCATTER_CONFIG.Y_RANGE = {
      {3, {0, 17.5}},  {4, {0, 17.5}}, {5, {0, 17.5}}, {6, {0, 17.5}},
      {7, {0, 17.5}},  {8, {0, 17.5}}, {9, {0, 17.5}}, {10, {0, 17.5}},
      {11, {0, 17.5}}, {12, {0, 10}},  {13, {0, 10}},  {14, {0, 10}},
      {15, {0, 10}}};
  gInstance.STRIP_SUM_SCATTER_CONFIG.SKIP_SAVGOL_PLOTS = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REJECT_NOISE = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REJECT_PILEUP = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.NOISE_THRESHOLD = 0.5;
  gInstance.STRIP_SUM_SCATTER_CONFIG.NOISE_MIN_STRIPS = 4;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REQUIRE_STRIP_16_BELOW_BEAM = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.END_STRIP_MAX = 1.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.PILEUP_THRESHOLD = 2;
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

  // ---------------------------------------------------------------------
  // Epochs. The flat block above stays as the fallback for tools that run
  // without an epoch set; the values below are what the pipeline actually
  // uses. `late` restates the flat block exactly, so enabling epochs does not
  // change the late calibration that reproduces the reference traces.
  //
  // The SOLARIS eras share one 1x64 map; CoMPASS run 37 used four 16-channel
  // boards with a different map (0bfa4e6:ProductionMode_37Cl/macros).
  // ---------------------------------------------------------------------
  std::map<std::pair<Int_t, Int_t>, TString> solaris_map =
      gInstance.channelMap64;

  RunEpoch late;
  late.name = "late";
  late.source = kSolaris;
  late.enabled = kTRUE;
  // Physics runs only: 84-96 are tuning runs. The flat run list this replaced
  // started at 84 and so included them.
  for (Int_t i = 97; i < 138; i++)
    late.runs.push_back(i);
  late.max_files = -1;
  late.n_boards = 1;
  late.n_channels = 64;
  late.channel_map = solaris_map;
  late.timing_ref_board = 0;
  late.timing_ref_board_channels = {8};
  late.do_board_sync = kFALSE;
  late.do_sort = kFALSE;
  late.event_time_window_us = 8.0;
  late.reference_channel = "Grid";
  late.reference_channel_min_adc = 500;
  late.reference_channel_max_adc = 4500;
  late.dedup_strategy = kLARGEST_ENERGY;
  late.has_cathode = kFALSE;
  late.strip_e_min_adc = 0.0;
  late.strip_e_max_adc = 4096.0;
  late.cathode_max_adc = 16384.0;
  late.grid_max_adc = 10000.0;
  late.strip0_max_adc = 1000.0;
  late.strip17_max_adc = 10000.0;
  late.left_even_max_adc = 2500.0;
  late.left_odd_max_adc = 8000.0;
  late.right_even_max_adc = 8000.0;
  late.right_odd_max_adc = 2500.0;

  // Same detector and DAQ as `late`, earlier in the run. Split the .sol files
  // for these runs with preprocess-sol before enabling.
  RunEpoch early = late;
  early.name = "early";
  early.enabled = kFALSE;
  early.runs.clear();
  for (Int_t i = 30; i <= 53; i++)
    early.runs.push_back(i);

  // CoMPASS era. Different digitiser: 4x16 channels, full 16384 ADC scale, a
  // cathode, board sync and sorting on, and no Grid amplitude window.
  RunEpoch compass;
  compass.name = "compass";
  // A different experiment years earlier that happens to reuse run numbers the
  // SOLARIS eras also use, so its output has to be tagged apart.
  compass.file_tag = "compass";
  compass.source = kCoMPASS;
  compass.enabled = kFALSE;
  compass.runs.push_back(37);
  compass.max_files = -1;
  compass.n_boards = 4;
  compass.n_channels = 16;
  compass.timing_ref_board = 1;
  compass.timing_ref_board_channels = {12, 0, 0, 0};
  compass.do_board_sync = kTRUE;
  compass.do_sort = kTRUE;
  compass.event_time_window_us = 8.0;
  compass.reference_channel = "Grid";
  compass.reference_channel_min_adc = 0;
  compass.reference_channel_max_adc = 16384;
  compass.dedup_strategy = kCLOSEST_TO_FIRST_GRID;
  compass.has_cathode = kTRUE;
  compass.strip_e_min_adc = 0.0;
  compass.strip_e_max_adc = 5500.0;
  compass.cathode_max_adc = 16384.0;
  compass.grid_max_adc = 16384.0;
  compass.strip0_max_adc = 16384.0;
  compass.strip17_max_adc = 16384.0;
  compass.left_even_max_adc = 16384.0;
  compass.left_odd_max_adc = 16384.0;
  compass.right_even_max_adc = 16384.0;
  compass.right_odd_max_adc = 16384.0;
  compass.channel_map = {
      {{0, 0}, "Cathode"}, {{0, 1}, ""},         {{0, 2}, "L2"},
      {{0, 3}, ""},        {{0, 4}, "Strip0"},   {{0, 5}, "100HzPulserBoard0"},
      {{0, 6}, "L6"},      {{0, 7}, ""},         {{0, 8}, "L1"},
      {{0, 9}, ""},        {{0, 10}, "L10"},     {{0, 11}, ""},
      {{0, 12}, "R2"},     {{0, 13}, "L14"},     {{0, 14}, ""},
      {{0, 15}, "Grid"},

      {{1, 0}, "L3"},      {{1, 1}, ""},         {{1, 2}, "R1"},
      {{1, 3}, ""},        {{1, 4}, "R6"},       {{1, 5}, "100HzPulserBoard1"},
      {{1, 6}, "R5"},      {{1, 7}, ""},         {{1, 8}, "L9"},
      {{1, 9}, ""},        {{1, 10}, "R9"},      {{1, 11}, ""},
      {{1, 12}, "R12"},    {{1, 13}, "R13"},     {{1, 14}, ""},
      {{1, 15}, "L15"},

      {{2, 0}, "R4"},      {{2, 1}, ""},         {{2, 2}, "L4"},
      {{2, 3}, ""},        {{2, 4}, "L7"},       {{2, 5}, "100HzPulserBoard2"},
      {{2, 6}, "L8"},      {{2, 7}, ""},         {{2, 8}, "R10"},
      {{2, 9}, ""},        {{2, 10}, "L12"},     {{2, 11}, ""},
      {{2, 12}, "L13"},    {{2, 13}, "L16"},     {{2, 14}, ""},
      {{2, 15}, "L11"},

      {{3, 0}, "L5"},      {{3, 1}, ""},         {{3, 2}, "R3"},
      {{3, 3}, ""},        {{3, 4}, "R8"},       {{3, 5}, "100HzPulserBoard3"},
      {{3, 6}, "R7"},      {{3, 7}, "SidETime"}, {{3, 8}, "R14"},
      {{3, 9}, ""},        {{3, 10}, "R11"},     {{3, 11}, "SidE"},
      {{3, 12}, "R16"},    {{3, 13}, "R15"},     {{3, 14}, ""},
      {{3, 15}, "Strip17"}};

  gInstance.EPOCHS.push_back(compass);
  gInstance.EPOCHS.push_back(early);
  gInstance.EPOCHS.push_back(late);
}

// Static initializer runs before main
struct InitGuard {
  InitGuard() { InitDatasetConfig(); }
};
static InitGuard gInit;

} // namespace Constants
