#include "Constants.hpp"
#include <DedupStrategy.hpp>
#include <RtypesCore.h>
#include <RunEpoch.hpp>

namespace Constants {

static DatasetConfig gInstance;

const DatasetConfig &cfg = gInstance;

void InitDatasetConfig() {
  // Data source: SOLARIS raw files, split into time chunks before building.
  gInstance.USE_SOLARIS_DATA = kTRUE;
  gInstance.SOL_BASE_DIR = "/labdata/MUSIC/New-37Cl/data_raw/";
  gInstance.SOL_SPLIT_DIR = "/labdata/MUSIC/New-37Cl/data_split/";
  gInstance.SOL_N_SPLIT_WORKERS = 32;
  gInstance.SOL_SPLIT_CHUNK_SECONDS = 120;
  gInstance.COMPASS_BASE_DIR = "/labdata/MUSIC/37Cl/";
  gInstance.N_CHUNKS = -1; // all
  gInstance.EPOCHS.push_back(MakeEpoch("late", RunRange(97, 137)));

  // Detector readout: one 64-channel board, channel map, per-channel limits.
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

  // Event building: seeded on the grid, one hit per channel per event, with
  // the long-end pole-zero undershoot of previous pulses removed per subfile.
  gInstance.REFERENCE_CHANNEL = "Grid";
  gInstance.REFERENCE_CHANNEL_MIN_ADC = 900;
  gInstance.REFERENCE_CHANNEL_MAX_ADC = 3600;
  gInstance.EVENT_TIME_WINDOW_US = 5.0;
  gInstance.DEDUP_STRATEGY = kLARGEST_ENERGY;
  gInstance.PULSE_HISTORY_CORRECTION = kTRUE;
  gInstance.PULSE_HISTORY_APPLY_MAX_US = 316.0;
  gInstance.MAX_FUSED_WORKERS = 16;
  gInstance.SKIP_EXISTING = kTRUE;
  gInstance.SAVE_PLOTS = kFALSE;
  gInstance.SAVE_SAMPLE_TRACES = 10;

  // Beam calibration and gating. Strips 0 and 17 are ignored throughout, so
  // the exit gate and the end-strip condition read strip 16.
  gInstance.IGNORE_STRIP_0 = kTRUE;
  gInstance.IGNORE_STRIP_17 = kTRUE;
  gInstance.BEAM_GATE_NSIGMA_X = 3;
  gInstance.BEAM_GATE_NSIGMA_Y = 3;
  gInstance.SIM_BEAM_FILE = "traces_37Cl_beam.root";

  // Strip-sum scatter: event-level cuts before any reaction is asked about.
  gInstance.STRIP_SUM_SCATTER_CONFIG.GATE_STRIP_X = 0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.GATE_STRIP_Y = 1;
  gInstance.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_X = 5.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.GATE_NSIGMA_Y = 5.0;

  // The tag, in the order the conditions are applied. Everything in sigma
  // resolves through the noise measured on the beam (jump sigma ~0.075 per
  // step, strip sigma ~0.065-0.073 at strips 15-16). The 40K residue stops
  // in strips 15-16, so no tail smoothness or zigzag condition is set: any
  // step ceiling there removes the residues. The falling-tail check instead
  // starts at each trace's own post-window peak.
  gInstance.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN = 2;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX = 10;
  gInstance.STRIP_SUM_SCATTER_CONFIG.REQUIRE_BEAM_UPSTREAM_OF_REAC = kTRUE;
  // Post-reaction smoothness at the shared 2.5 sigma to strip 12.
  gInstance.STRIP_SUM_SCATTER_CONFIG.REQUIRE_SMOOTHNESS = kTRUE;
  // Persistence: 3 strips after the reaction each more than 2 sigma of their
  // spread above the beam. The 40K excess is ~0.25, about 4 sigma.
  gInstance.STRIP_SUM_SCATTER_CONFIG.POST_ABOVE_NSIGMA = 2.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.POST_ABOVE_STRIPS = 3;
  // Falling tail after the peak: from the largest deposit in the post window
  // to strip 16 no rise above 2 sigma of a step (0.15). A residue peaks once
  // and declines to its stop; the tags that leaked at reaction strip 6 fell
  // back to the beam and rose again at strips 12-14 before dropping at 16.
  gInstance.STRIP_SUM_SCATTER_CONFIG.TAIL_FALL_AFTER_PEAK = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.TAIL_RISE_NSIGMA = 3.0;
  // No return: after the peak, once back at or below the beam the trace must
  // not rise above 2 sigma of its spread over the beam again. The tags that
  // leaked at reaction strip 6 returned to 1.0 at strips 10-11 and sat at
  // 1.15-1.25 over strips 12-14, a rise too slow for the per-step ceiling.
  gInstance.STRIP_SUM_SCATTER_CONFIG.TAIL_RETURN_NSIGMA = 2.0;
  gInstance.STRIP_SUM_SCATTER_CONFIG.TAIL_RERISE_NSIGMA = 2.0;
  // Below the beam after the reaction, and already on the way down: strips
  // 15 and 16 both more than 1 sigma of their spread under the beam. A
  // residue crosses the beam at strips 10-13 and declines to the exit (the
  // templates at reaction strips 2-6: 0.87 at 15, 0.45-0.55 at 16). Strip 16
  // alone let through a class that sits 10-15% high all the way to strip 15
  // and drops by half in one strip, which dominated the tags from strip 7
  // on; those have strip 15 at 1.1 and fail here.
  gInstance.STRIP_SUM_SCATTER_CONFIG.LATE_STRIP_BELOW_NSIGMA = {{15, 1.0},
                                                                {16, 1.0}};
  gInstance.STRIP_SUM_SCATTER_CONFIG.REQUIRE_STRIP_16_BELOW_BEAM = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.END_STRIP_MAX = 1.0;

  // The plane the regions are fitted in: x sums strips 1-13, y the four
  // strips after the reaction, capped at strip 16. A change here reprojects
  // the cache from the reservoir rather than refilling.
  gInstance.STRIP_SUM_SCATTER_CONFIG.X_LO = 1;
  gInstance.STRIP_SUM_SCATTER_CONFIG.X_HI = 13;
  gInstance.STRIP_SUM_SCATTER_CONFIG.POST_TRIGGER_SUM_STRIPS = 4;
  gInstance.STRIP_SUM_SCATTER_CONFIG.POST_WINDOW_LAST_STRIP = 16;
  // Per-event ratio normalisation of y left off: it moves every region.
  gInstance.STRIP_SUM_SCATTER_CONFIG.Y_RATIO_TO_UPSTREAM = kFALSE;

  // Regions and display.
  gInstance.STRIP_SUM_SCATTER_CONFIG.AN_REGION_MODE =
      StripSumScatterConfig::AN_REGION_ALL_TAGGED;
  // Per-strip override: AN_REGION_ALL_TAGGED counts every event the tag
  // leaves as (a,n) with no fit, for strips where the tag conditions alone
  // isolate the residues, e.g. {{5,
  // StripSumScatterConfig::AN_REGION_ALL_TAGGED}}.
  gInstance.STRIP_SUM_SCATTER_CONFIG.AN_REGION_MODE_STRIPS = {};
  gInstance.STRIP_SUM_SCATTER_CONFIG.REGION_CUT_REDRAW = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.X_DISPLAY_MIN = 11;
  gInstance.STRIP_SUM_SCATTER_CONFIG.X_DISPLAY_MAX = 15;
  gInstance.STRIP_SUM_SCATTER_CONFIG.Y_DISPLAY_MIN = 2;
  gInstance.STRIP_SUM_SCATTER_CONFIG.Y_DISPLAY_MAX = 10;
  gInstance.STRIP_SUM_SCATTER_CONFIG.CANDIDATE_REAC_STRIP = 6;

  // Plotting and workers.
  gInstance.STRIP_SUM_SCATTER_CONFIG.MAX_STRIP_SUM_WORKERS = 8;
  gInstance.STRIP_SUM_SCATTER_CONFIG.SKIP_RUN_PLOTS = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.SKIP_SAVGOL_PLOTS = kTRUE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.PLOT_PARITY_REJECTED_GRID = kFALSE;
  gInstance.STRIP_SUM_SCATTER_CONFIG.RERUN_SIM = kFALSE;

  // Cross section: 37Cl on helium at 400 Torr, the (a,n) channel only, with
  // TALYS Hauser-Feshbach curves for comparison.
  gInstance.CROSS_SECTION_CONFIG.TARGET_GAS = kHELIUM;
  gInstance.CROSS_SECTION_CONFIG.GAS_PRESSURE_TORR = 400.0;
  gInstance.CROSS_SECTION_CONFIG.BEAM_A = 37;
  gInstance.CROSS_SECTION_CONFIG.BEAM_Z = 17;
  gInstance.CROSS_SECTION_CONFIG.BEAM_ELEMENT = "Cl";
  gInstance.CROSS_SECTION_CONFIG.BEAM_SIM_FILE = "traces_37Cl_beam.root";
  gInstance.CROSS_SECTION_CONFIG.XS_STRIP_MIN = 2;
  gInstance.CROSS_SECTION_CONFIG.XS_STRIP_MAX = 8;
  gInstance.CROSS_SECTION_CONFIG.EFFECTIVE_ENERGY = kFALSE;
  gInstance.CROSS_SECTION_CONFIG.EPOCHS = {"late"};
  gInstance.CROSS_SECTION_CONFIG.CHANNELS = {
      {"an", "(#alpha, n)", {"n"}, {}, ""}};
  gInstance.CROSS_SECTION_CONFIG.TALYS_MODELS = {
      {"TALYS HF, Avrigeanu", {"alphaomp 6"}},
      {"TALYS HF, McFadden-Satchler", {"alphaomp 2"}}};
}

// Static initializer runs before main
struct InitGuard {
  InitGuard() { InitDatasetConfig(); }
};
static InitGuard gInit;

} // namespace Constants
