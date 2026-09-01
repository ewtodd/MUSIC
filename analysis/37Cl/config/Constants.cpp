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

  gInstance.N_CHUNKS = 1;
  gInstance.SIM_BEAM_FILE = "traces_37Cl_beam.root";

  gInstance.N_BOARDS = 1;
  gInstance.N_CHANNELS = 64;
  gInstance.TIMING_REF_BOARD = 0;
  gInstance.TIMING_REF_BOARD_CHANNELS = {8};

  gInstance.SAVE_PLOTS = kTRUE;
  gInstance.SKIP_EXISTING = kTRUE;
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
  // 6 keeps the y-sum span main has always used; the tooling default is 3.
  gInstance.STRIP_SUM_SCATTER_CONFIG.POST_TRIGGER_SUM_STRIPS = 6;
  gInstance.STRIP_SUM_SCATTER_CONFIG.CANDIDATE_REAC_STRIP = 5;
  gInstance.STRIP_SUM_SCATTER_CONFIG.X_DISPLAY_MIN = 11;
  gInstance.STRIP_SUM_SCATTER_CONFIG.X_DISPLAY_MAX = 21;
  gInstance.STRIP_SUM_SCATTER_CONFIG.Y_DISPLAY_RANGE = {
      {2, {3.3, 8.7}},  {3, {3.3, 8.7}},  {4, {3.3, 8.7}},  {5, {3.3, 8.7}},
      {6, {3.3, 8.7}},  {7, {3.3, 8.7}},  {8, {3.3, 8.7}},  {9, {3.3, 8.7}},
      {10, {3.3, 8.7}}, {11, {3.3, 8.7}}, {12, {2.8, 7.2}}, {13, {2.2, 5.8}},
      {14, {1.7, 4.3}}, {15, {1.1, 2.9}}};
  // Saved cuts win when RegionCuts.root has them, so the run repeats headless.
  // A strip with no saved cut still prompts, which is how the first pass fills
  // the file. Flip to kTRUE only to deliberately redraw over saved cuts.
  gInstance.STRIP_SUM_SCATTER_CONFIG.REGION_CUT_REDRAW = kFALSE;
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

  std::map<std::pair<Int_t, Int_t>, TString> solaris_map =
      gInstance.channelMap64;

  RunEpoch late;
  late.name = "late";
  late.source = kSolaris;
  late.enabled = kTRUE;
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

  RunEpoch early = late;
  early.name = "early";
  early.enabled = kFALSE;
  early.runs.clear();
  for (Int_t i = 30; i <= 53; i++)
    early.runs.push_back(i);

  //  gInstance.EPOCHS.push_back(early);
  gInstance.EPOCHS.push_back(late);
}

// Static initializer runs before main
struct InitGuard {
  InitGuard() { InitDatasetConfig(); }
};
static InitGuard gInit;

} // namespace Constants
