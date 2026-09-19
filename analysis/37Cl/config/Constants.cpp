#include "Constants.hpp"
#include <DedupStrategy.hpp>
#include <RtypesCore.h>
#include <RunEpoch.hpp>

namespace Constants {

static DatasetConfig gInstance;

const DatasetConfig &cfg = gInstance;

void InitDatasetConfig() {
  gInstance.USE_SOLARIS_DATA = kTRUE;
  gInstance.SOL_BASE_DIR = "/labdata/MUSIC/New-37Cl/data_raw/";
  gInstance.SOL_SPLIT_DIR = "/labdata/MUSIC/New-37Cl/data_split/";
  gInstance.SOL_N_SPLIT_WORKERS = 32;
  gInstance.SOL_SPLIT_CHUNK_SECONDS = 120;
  gInstance.COMPASS_BASE_DIR = "/labdata/MUSIC/37Cl/";
  gInstance.N_CHUNKS = -1;

  gInstance.N_BOARDS = 1;
  gInstance.N_CHANNELS = 64;
  gInstance.HAS_CATHODE = kFALSE;
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
  gInstance.STRIP0_MAX_ADC = 1000;
  gInstance.STRIP17_MAX_ADC = 10000;
  gInstance.GRID_MAX_ADC = 10000;
  gInstance.RIGHT_ODD_MAX_ADC = 2500;
  gInstance.RIGHT_EVEN_MAX_ADC = 8000;
  gInstance.LEFT_ODD_MAX_ADC = 8000;
  gInstance.LEFT_EVEN_MAX_ADC = 2500;
  gInstance.STRIP_E_MAX_ADC = 4096.0;
  gInstance.TOTAL_E_MAX_ADC = 60000.0;
  gInstance.STRIP_DE_MIN_NORMED = 0;
  gInstance.STRIP_DE_MAX_NORMED = 4;

  gInstance.REFERENCE_CHANNEL = "Grid";
  gInstance.REFERENCE_CHANNEL_MIN_ADC = 900;
  gInstance.REFERENCE_CHANNEL_MAX_ADC = 3600;
  gInstance.EVENT_TIME_WINDOW_US = 5.0;
  gInstance.DEDUP_STRATEGY = kLARGEST_ENERGY;

  gInstance.PULSE_HISTORY_CORRECTION = kTRUE;
  gInstance.PULSE_HISTORY_GROUPS.long_left.enabled = kTRUE;
  gInstance.PULSE_HISTORY_GROUPS.long_left.kernel = kPulseHistoryBinned;
  gInstance.PULSE_HISTORY_GROUPS.long_right.enabled = kTRUE;
  gInstance.PULSE_HISTORY_GROUPS.long_right.kernel = kPulseHistoryBinned;
  gInstance.PULSE_HISTORY_GROUPS.short_left.enabled = kTRUE;
  gInstance.PULSE_HISTORY_GROUPS.short_left.kernel = kPulseHistoryForm;
  gInstance.PULSE_HISTORY_APPLY_MAX_US = 316.0;
  gInstance.PULSE_HISTORY_AMP_BINS = 4;

  gInstance.MAX_FUSED_WORKERS = 16;
  gInstance.SKIP_EXISTING = kTRUE;
  gInstance.SAVE_SAMPLE_TRACES = 10;

  gInstance.IGNORE_STRIP_0 = kFALSE;
  gInstance.IGNORE_STRIP_17 = kTRUE;
  gInstance.BEAM_GATE_NSIGMA_X = 5;
  gInstance.BEAM_GATE_NSIGMA_Y = 5;
  gInstance.SIM_BEAM_FILE = "traces_37Cl_beam.root";

  gInstance.STRIP_SUM_SCATTER_CONFIG.GATE_STRIP_X = 0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.GATE_STRIP_Y = 1;
  gInstance.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_X = 5.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_Y = 5.0;

  gInstance.STRIP_SUM_SCATTER_CONFIG.PILEUP_NSIGMA = 8.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.PILEUP_MIN_STRIPS = 2;
  gInstance.STRIP_SUM_SCATTER_CONFIG.NOISE_NSIGMA = 12.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.NOISE_MIN_STRIPS = 2;

  gInstance.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN = 2;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX = 8;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REQUIRE_BEAM_UPSTREAM_OF_REAC = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.BEAM_UPSTREAM_NSIGMA = 2.0;

  gInstance.STRIP_SUM_SCATTER_CONFIG.REAC_JUMP_NSIGMA = 2.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.POST_ABOVE_NSIGMA = 1.5;
  gInstance.STRIP_SUM_SCATTER_CONFIG.POST_ABOVE_STRIPS = 3;

  gInstance.STRIP_SUM_SCATTER_CONFIG.SMOOTHNESS_NSIGMA = 6.0;

  gInstance.STRIP_SUM_SCATTER_CONFIG.TAIL_RISE_NSIGMA = 0.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.TAIL_RETURN_NSIGMA = 2.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.TAIL_RERISE_NSIGMA = 4.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REQUIRE_STRIP_16_BELOW_BEAM = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.END_STRIP_NSIGMA = 3.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.TAIL_CLIFF_MAX_FRACTION = 0.5;

  gInstance.STRIP_SUM_SCATTER_CONFIG.X_LO = 1;
  gInstance.STRIP_SUM_SCATTER_CONFIG.X_HI = 13;
  gInstance.STRIP_SUM_SCATTER_CONFIG.POST_TRIGGER_SUM_STRIPS = 4;
  gInstance.STRIP_SUM_SCATTER_CONFIG.POST_WINDOW_LAST_STRIP = 16;

  gInstance.STRIP_SUM_SCATTER_CONFIG.Y_RATIO_TO_UPSTREAM = kTRUE;

  gInstance.STRIP_SUM_SCATTER_CONFIG.AN_REGION_MODE =
      StripSumScatterConfig::AN_REGION_ALL_TAGGED;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REGION_CUT_REDRAW = kFALSE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.X_DISPLAY_MIN = 10;
  gInstance.STRIP_SUM_SCATTER_CONFIG.X_DISPLAY_MAX = 18;
  gInstance.STRIP_SUM_SCATTER_CONFIG.Y_DISPLAY_MIN = 3;
  gInstance.STRIP_SUM_SCATTER_CONFIG.Y_DISPLAY_MAX = 7;
  gInstance.STRIP_SUM_SCATTER_CONFIG.CANDIDATE_REAC_STRIP = 3;

  gInstance.STRIP_SUM_SCATTER_CONFIG.MAX_STRIP_SUM_WORKERS = 8;
  gInstance.STRIP_SUM_SCATTER_CONFIG.SKIP_SAVGOL_PLOTS = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.RERUN_SIM = kFALSE;

  gInstance.STRIP_SUM_SCATTER_CONFIG.CUT_VARIATION_NSIGMA_STEP = 0.2;

  gInstance.CROSS_SECTION_CONFIG.TARGET_GAS = kHELIUM;
  gInstance.CROSS_SECTION_CONFIG.GAS_PRESSURE_TORR = 400.0;
  gInstance.CROSS_SECTION_CONFIG.BEAM_A = 37;
  gInstance.CROSS_SECTION_CONFIG.BEAM_Z = 17;
  gInstance.CROSS_SECTION_CONFIG.BEAM_ELEMENT = "Cl";
  gInstance.CROSS_SECTION_CONFIG.BEAM_SIM_FILE = "traces_37Cl_beam.root";
  gInstance.CROSS_SECTION_CONFIG.XS_STRIP_MIN = 2;
  gInstance.CROSS_SECTION_CONFIG.XS_STRIP_MAX = 8;
  gInstance.CROSS_SECTION_CONFIG.EFFECTIVE_ENERGY = kTRUE;
  gInstance.CROSS_SECTION_CONFIG.PRELIMINARY = kTRUE;
  gInstance.CROSS_SECTION_CONFIG.EPOCHS = {"late"};
  gInstance.CROSS_SECTION_CONFIG.CHANNELS = {
      {"an", "(#alpha, n)", {"n"}, {}, ""}};
  gInstance.CROSS_SECTION_CONFIG.TALYS_MODELS = {
      {"TALYS HF, McFadden-Satchler", {"alphaomp 2"}}};

  RunEpoch early = MakeEpoch("early", RunRange(30, 40));
  early.split_chunk_seconds = 15;
  gInstance.EPOCHS.push_back(early);
  RunEpoch late = MakeEpoch("late", RunRange(97, 137));
  late.pulse_history.long_left.tau_us = 25;
  late.pulse_history.long_right.tau_us = 10;
  late.pulse_history.short_left.tau_us = 25;
  gInstance.EPOCHS.push_back(late);
}

struct InitGuard {
  InitGuard() { InitDatasetConfig(); }
};
static InitGuard gInit;

} // namespace Constants
