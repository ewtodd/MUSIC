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

  gInstance.STRIP_SUM_SCATTER_CONFIG.GATE_STRIP_X = 0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.GATE_STRIP_Y = 1;
  gInstance.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_X = 5.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_Y = 5.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.BOTH_MULT_MAX = 3;
  gInstance.STRIP_SUM_SCATTER_CONFIG.BOTH_MULT_COUNT_TO = 16;
  gInstance.STRIP_SUM_SCATTER_CONFIG.ALT_DECODE_REGION_TRACES = kFALSE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.SKIP_SAVGOL_PLOTS = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN = 2;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX = 12;
  gInstance.STRIP_SUM_SCATTER_CONFIG.CANDIDATE_REAC_STRIP = 2;
  gInstance.STRIP_SUM_SCATTER_CONFIG.POST_TRIGGER_SUM_STRIPS = 6;
  gInstance.STRIP_SUM_SCATTER_CONFIG.POST_WINDOW_LAST_STRIP = 14;
  gInstance.STRIP_SUM_SCATTER_CONFIG.X_DISPLAY_MIN = 13;
  gInstance.STRIP_SUM_SCATTER_CONFIG.X_DISPLAY_MAX = 21;
  gInstance.STRIP_SUM_SCATTER_CONFIG.Y_DISPLAY_MIN = 0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.Y_DISPLAY_MAX = 10;
  gInstance.STRIP_SUM_SCATTER_CONFIG.XBINS = 1500;
  gInstance.STRIP_SUM_SCATTER_CONFIG.YBINS = 2000;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REGION_CUT_REDRAW = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REAC_JUMP_NSIGMA = 1.0;
  // compute-regions (a separate, read-only pass over the scatter cache)
  // fits a bivariate Gaussian mixture per strip and saves these ellipses as
  // the region cuts, overwriting the per-cut files. Hand-drawn cuts stay in
  // the RegionCuts.root and in region_cuts_hand_<date>/.
  gInstance.STRIP_SUM_SCATTER_CONFIG.AN_REGION_NSIGMA = 5.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.AA_REGION_NSIGMA = 5.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REQUIRE_BEAM_UPSTREAM_OF_REAC = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.BEAM_UPSTREAM_NSIGMA = 2.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.Y_DISPLAY_RANGE = {
      {3, {4, 9}},     {4, {5, 7.5}},    {5, {5, 7.5}},    {6, {5, 7.5}},
      {7, {5, 7.5}},   {8, {5, 7.5}},    {9, {4, 6.5}},    {10, {3, 5.5}},
      {11, {2, 4.5}},  {12, {1.2, 3.2}}, {13, {0.5, 1.8}}, {14, {0.5, 1.8}},
      {15, {0.5, 1.8}}};

  gInstance.CROSS_SECTION_CONFIG.TARGET_GAS = kHELIUM;
  gInstance.CROSS_SECTION_CONFIG.GAS_PRESSURE_TORR = 555.0;
  gInstance.CROSS_SECTION_CONFIG.BEAM_A = 87;
  gInstance.CROSS_SECTION_CONFIG.BEAM_Z = 37;
  gInstance.CROSS_SECTION_CONFIG.BEAM_ELEMENT = "Rb";
  // TALYS Hauser-Feshbach for the plot, written by talys-xs. The paper used
  // the Atomki-V2 alpha potential, which TALYS 2.2 does not carry; its
  // default Avrigeanu (2014) potential, alphaomp 6, tracks it closely at
  // these energies (ApJ 983:142 Figure 1) and is the shape the effective
  // energies come from. McFadden-Satchler (alphaomp 2) has the most
  // different energy dependence of TALYS's eight and is the check on them.
  gInstance.CROSS_SECTION_CONFIG.TALYS_MODELS = {
      {"TALYS HF, Avrigeanu", {"alphaomp 6"}},
      {"TALYS HF, McFadden-Satchler", {"alphaomp 2"}}};
  gInstance.CROSS_SECTION_CONFIG.BEAM_SIM_FILE = "traces_87Rb_beam.root";
  gInstance.CROSS_SECTION_CONFIG.XS_STRIP_MIN = 2;
  gInstance.CROSS_SECTION_CONFIG.XS_STRIP_MAX = 11;
  // One channel: the inclusive (a,xn), (a,n) + (a,2n), against ApJ 983:142
  // Table 2 -- effective centre-of-mass energy [MeV] with its +/-
  // uncertainties (the strip's extent about it, asymmetric because the
  // thick-target-yield correction puts E_eff above the midpoint), sigma(a,xn)
  // [mb] and the quadrature sum of its statistical and systematic
  // uncertainties [mb]. Ten strips, 13.01 down to 8.09 MeV: rows are strips
  // 2 to 11 (the published analysis tree has A2..A11).
  gInstance.CROSS_SECTION_CONFIG.CHANNELS = {
      {"an",
       "(#alpha, xn)",
       {"n", "2n"},
       {
           {13.01, 0.26, 0.28, 420.0, 17.0},
           {12.49, 0.24, 0.30, 355.0, 19.0},
           {11.94, 0.25, 0.30, 291.0, 11.0},
           {11.41, 0.23, 0.32, 192.0, 8.0},
           {10.87, 0.22, 0.33, 109.0, 6.0},
           {10.32, 0.22, 0.33, 67.4, 3.8},
           {9.77, 0.21, 0.44, 32.9, 2.0},
           {9.22, 0.19, 0.38, 10.7, 1.0},
           {8.66, 0.18, 0.39, 3.77, 0.66},
           {8.09, 0.18, 0.39, 0.92, 0.31},
       },
       "Foug#grave{e}res et al. (2025)"}};

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
