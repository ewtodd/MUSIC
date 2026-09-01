#include "Constants.hpp"
#include <RtypesCore.h>

namespace Constants {

static DatasetConfig gInstance;

const DatasetConfig &cfg = gInstance;

void InitDatasetConfig() {
  gInstance.USE_SOLARIS_DATA = kFALSE;
  gInstance.COMPASS_BASE_DIR = "/labdata/MUSIC/87Rb/";
  gInstance.RUN_NUMBERS = {16, 20};
  gInstance.N_CHUNKS = -1;
  gInstance.SIM_BEAM_FILE = "traces_87Rb_beam.root";

  gInstance.N_BOARDS = 4;
  gInstance.N_CHANNELS = 16;

  gInstance.REFERENCE_CHANNEL = "NONE";
  gInstance.EVENT_TIME_WINDOW_US = 8;
  gInstance.DEDUP_STRATEGY = kDISCARD;

  gInstance.SKIP_EXISTING = kTRUE;
  gInstance.SAVE_PLOTS = kFALSE;

  gInstance.TIMING_MIN_ENERGY = 300;
  gInstance.TIMING_MAX_ENERGY = 1500;
  gInstance.TIMING_DO_BOARD_SYNC = kFALSE;
  gInstance.TIMING_DO_SORT = kTRUE;

  gInstance.IGNORE_SHORT_STRIPS = kFALSE;
  gInstance.REJECT_FLAGGED_EVENTS = kTRUE;
  gInstance.MAX_FUSED_WORKERS = 16;
  gInstance.MAX_GPU_CONCURRENT_SORTS = 20;

  gInstance.STRIP_SUM_SCATTER_CONFIG.BOTH_MULT_MAX = 3;
  gInstance.STRIP_SUM_SCATTER_CONFIG.BOTH_MULT_COUNT_TO = 16;
  gInstance.STRIP_SUM_SCATTER_CONFIG.ALT_DECODE_REGION_TRACES = kFALSE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.SKIP_SAVGOL_PLOTS = kTRUE;
  // 6 keeps the y-sum span main has always used; the tooling default is 3.
  gInstance.STRIP_SUM_SCATTER_CONFIG.POST_TRIGGER_SUM_STRIPS = 6;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN = 3;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX = 15;
  gInstance.STRIP_SUM_SCATTER_CONFIG.CANDIDATE_REAC_STRIP = 5;
  gInstance.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_X = 3.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_Y = 3.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.X_DISPLAY_MIN = 13;
  gInstance.STRIP_SUM_SCATTER_CONFIG.X_DISPLAY_MAX = 21;
  gInstance.STRIP_SUM_SCATTER_CONFIG.Y_DISPLAY_MIN = 0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.Y_DISPLAY_MAX = 10;
  gInstance.STRIP_SUM_SCATTER_CONFIG.XBINS = 1500;
  gInstance.STRIP_SUM_SCATTER_CONFIG.YBINS = 2000;
  // Saved cuts win when RegionCuts.root has them, so the run repeats headless.
  // A strip with no saved cut still prompts, which is how the first pass fills
  // the file. Flip to kTRUE only to deliberately redraw over saved cuts.
  gInstance.STRIP_SUM_SCATTER_CONFIG.REGION_CUT_REDRAW = kFALSE;
  // Per-strip y zoom. Display-only now (applied with SetRangeUser over the
  // fixed build range), so changing it never forces a refill.
  gInstance.STRIP_SUM_SCATTER_CONFIG.Y_DISPLAY_RANGE = {
      {3, {4, 9}},  {4, {5, 7}},  {5, {5.5, 7.5}}, {6, {5, 7}},  {7, {5, 7}},
      {8, {5, 7}},  {9, {5, 7}},  {10, {5, 7}},    {11, {5, 7}}, {12, {3.5, 6}},
      {13, {2, 5}}, {14, {2, 5}}, {15, {1, 4}}};

  gInstance.STRIP_DE_MIN_NORMED = 0;
  gInstance.STRIP_DE_MAX_NORMED = 4;
  gInstance.STRIP_E_MAX_ADC = 12000;
  gInstance.TOTAL_E_MAX_ADC = 15 * gInstance.STRIP_E_MAX_ADC;

  gInstance.channelMap = {
      {{0, 0}, "Cathode"}, {{0, 1}, ""},       {{0, 2}, "L2"},
      {{0, 3}, ""},        {{0, 4}, "Strip0"}, {{0, 5}, ""},
      {{0, 6}, "L6"},      {{0, 7}, ""},       {{0, 8}, "L1"},
      {{0, 9}, ""},        {{0, 10}, "L10"},   {{0, 11}, ""},
      {{0, 12}, "R2"},     {{0, 13}, "L14"},   {{0, 14}, ""},
      {{0, 15}, "Grid"},

      {{1, 0}, "L3"},      {{1, 1}, ""},       {{1, 2}, "R1"},
      {{1, 3}, ""},        {{1, 4}, "R6"},     {{1, 5}, ""},
      {{1, 6}, "R5"},      {{1, 7}, ""},       {{1, 8}, "L9"},
      {{1, 9}, ""},        {{1, 10}, "R9"},    {{1, 11}, ""},
      {{1, 12}, "R12"},    {{1, 13}, "R13"},   {{1, 14}, ""},
      {{1, 15}, "L15"},

      {{2, 0}, "R4"},      {{2, 1}, ""},       {{2, 2}, "L4"},
      {{2, 3}, ""},        {{2, 4}, "L7"},     {{2, 5}, ""},
      {{2, 6}, "L8"},      {{2, 7}, ""},       {{2, 8}, "R10"},
      {{2, 9}, ""},        {{2, 10}, "L12"},   {{2, 11}, ""},
      {{2, 12}, "L13"},    {{2, 13}, "L16"},   {{2, 14}, ""},
      {{2, 15}, "L11"},

      {{3, 0}, "L5"},      {{3, 1}, ""},       {{3, 2}, "R3"},
      {{3, 3}, "SidE"},    {{3, 4}, "R8"},     {{3, 5}, ""},
      {{3, 6}, "R7"},      {{3, 7}, ""},       {{3, 8}, "R14"},
      {{3, 9}, ""},        {{3, 10}, "R11"},   {{3, 11}, ""},
      {{3, 12}, "R16"},    {{3, 13}, "R15"},   {{3, 14}, ""},
      {{3, 15}, "Strip17"}};
}

// Static initializer runs before main
struct InitGuard {
  InitGuard() { InitDatasetConfig(); }
};
static InitGuard gInit;

} // namespace Constants
