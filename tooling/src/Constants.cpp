#include "Constants.hpp"

void StripSumScatterConfig::SetDefaults() {
  PURE_BEAM_GATE = PURE_BEAM_GATE_S0_S1;

  POST_TRIGGER_SUM_STRIPS = 3;
  MAX_STRIP_SUM_WORKERS = 12;

  REACTION_STRIP_MIN = 2;
  REACTION_STRIP_MAX = 15;

  REQUIRE_SMOOTHNESS_END_STRIP = 12;
  REQUIRE_SMOOTHNESS_MAX_STEP = 1.2;

  REAC_JUMP_MIN = 0.1;
  REAC_JUMP_MAX = 2.0;
  END_STRIP_MAX = 1.0;

  PILEUP_THRESHOLD = 1.75;
  NOISE_THRESHOLD = 0.85;

  REJECT_NOISE = kTRUE;
  NOISE_THRESH_PY = 0.4;
  NOISE_MIN_STRIPS = 1;

  REJECT_PILEUP = kTRUE;
  PILEUP_THRESH_PY = 1.3;
  PILEUP_MIN_STRIPS = 1;

  REGION_CUT_REDRAW = kFALSE;
  BOTH_MULT_MAX = -1;      // disabled by default
  BOTH_MULT_COUNT_TO = 16; // whole trace unless narrowed
  REJECT_OFFBEAM = kFALSE;
  OFFBEAM_DIST = 0.3;
  OFFBEAM_MIN_STRIPS = 4;

  TRIGGER_NSIGMA = 5.0;
  TRIGGER_CFD_FRAC = 0.30;
  PLATEAU_POST = 3;
  CLUSTER_SMOOTH_WINDOW = 1;
  SEED_HALF_BINS = 40;

  SAVGOL_HALF = 2;

  TRACES_PER_CLASS = 40;

  X_LO = 1;
  X_HI = 16;

  GATE_STRIP_X = 1;
  GATE_STRIP_Y = 2;
  GATE_NSIGMA_X = 3.5;
  GATE_NSIGMA_Y = 3.5;
  GATE_MIN = 0.0;
  GATE_MAX = 3.0;
  GATE_BINS = 240;

  X_DISPLAY_MIN = 14;
  X_DISPLAY_MAX = 26;
  XBINS = 1500;
  Y_DISPLAY_MIN = 0;
  Y_DISPLAY_MAX = 20;
  YBINS = 1500;

  SAMPLE_MAX_POINTS = 2000000;

  RERUN_SIM = kFALSE;
  CANDIDATE_REAC_STRIP = 3;

  REQUIRE_SMOOTHNESS = kTRUE;
  REQUIRE_GATE_S3_S4 = kFALSE;
  REQUIRE_GATE_S5_S6 = kFALSE;
  SKIP_SAVGOL_PLOTS = kFALSE;
  REQUIRE_STRIP_16_BELOW_BEAM = kFALSE;
  ALT_DECODE_REGION_TRACES = kFALSE;
}

DatasetConfig::DatasetConfig() {
  // Data source
  USE_SOLARIS_DATA = kFALSE;
  SOL_SPLIT_CHUNK_SECONDS = 30.0;
  SOL_N_SPLIT_WORKERS = 32;
  N_CHUNKS = -1;

  // Hardware layout
  N_BOARDS = 4;
  N_CHANNELS = 16;
  TIMING_REF_BOARD = 0;

  // Timing
  TIMING_MIN_ENERGY = 0;
  TIMING_MAX_ENERGY = 16384;
  TIMING_OVERLAP_MARGIN_S = 5;
  TIMING_THRESH_DT_US = 150.0;
  TIMING_MAX_ABS_SHIFT_S = 0.25;
  TIMING_SHIFT_COARSE_STEP_US = 1000.0;
  TIMING_SHIFT_FINE_STEP_US = 0.5;
  TIMING_SHIFT_FINE_HALF_WIDTH_US = 1000.0;
  TIMING_SHIFT_MIN_NPTS = 10;
  TIMING_SHIFT_MAX_SCAN_CANDIDATES = 2000000;
  TIMING_DO_BOARD_SYNC = kFALSE;
  TIMING_DO_SORT = kFALSE;

  REJECT_FLAGGED_EVENTS = kFALSE;

  IGNORE_SHORT_STRIPS = kFALSE;
  IGNORE_STRIP_0 = kFALSE;
  IGNORE_STRIP_17 = kFALSE;

  HAS_CATHODE = kTRUE;
  HAS_GRID = kTRUE;
  HAS_STRIP0 = kTRUE;
  HAS_STRIP17 = kTRUE;

  SKIP_EXISTING = kTRUE;
  SAVE_PLOTS = kTRUE;

  SKIP_CALIBRATION = kFALSE;

  SAVE_SAMPLE_TRACES = 0;

  MAX_FUSED_WORKERS = 12;

  REFERENCE_CHANNEL = "Grid";
  EVENT_TIME_WINDOW_US = 8.0;
  DEDUP_STRATEGY = kLARGEST_ENERGY;

  USE_GPU_ACCELERATION = kTRUE;
  MAX_GPU_CONCURRENT_SORTS = 20;

  STRIP_DE_OVERVIEW_MIN_NORMED = 0.8;
  STRIP_DE_OVERVIEW_MAX_NORMED = 5;
  STRIP_DE_MIN_NORMED = 0.8;
  STRIP_DE_MAX_NORMED = 1.3;
  CATHODE_E_MAX_NORMED = 300;
  TOTAL_E_MIN_NORMED = 10.0;
  TOTAL_E_MAX_NORMED = 400.0;

  STRIP_SUM_SCATTER_CONFIG.SetDefaults();

  STRIP_E_MIN_ADC = 0.0;
  STRIP_E_MAX_ADC = 4096.0;
  TOTAL_E_MIN_ADC = 0.0;
  TOTAL_E_MAX_ADC = 60000.0;

  GRID_MIN_ADC = 0.0;
  GRID_MAX_ADC = 16384.0;

  REFERENCE_CHANNEL_MIN_ADC = 0.0;
  REFERENCE_CHANNEL_MAX_ADC = 16384.0;

  STRIP0_MAX_ADC = 16384.0;
  STRIP17_MAX_ADC = 16384.0;
  CATHODE_MAX_ADC = 16384.0;

  LEFT_EVEN_MAX_ADC = 16384.0;
  LEFT_ODD_MAX_ADC = 16384.0;
  RIGHT_EVEN_MAX_ADC = 16384.0;
  RIGHT_ODD_MAX_ADC = 16384.0;
}

namespace Constants {

// The epoch currently being processed, or null outside epoch work (and for
// datasets that declare none). Every accessor below falls back to the flat cfg
// block when it is null, so behaviour without epochs is unchanged.
static const RunEpoch *gActiveEpoch = nullptr;

void SetActiveEpoch(const RunEpoch *epoch) { gActiveEpoch = epoch; }
const RunEpoch *GetActiveEpoch() { return gActiveEpoch; }

const RunEpoch *EpochForRun(Int_t run) {
  // Only untagged epochs are addressable by run number. A tagged epoch exists
  // precisely because its run numbers collide with another era's, so answering
  // from the number alone would be a coin flip; those are addressed by tag.
  for (Int_t e = 0; e < Int_t(cfg.EPOCHS.size()); e++) {
    const RunEpoch &ep = cfg.EPOCHS[e];
    if (ep.file_tag.Length() > 0)
      continue;
    for (Int_t r = 0; r < Int_t(ep.runs.size()); r++)
      if (ep.runs[r] == run)
        return &ep;
  }
  return nullptr;
}

Int_t ActiveNBoards() {
  return gActiveEpoch ? gActiveEpoch->n_boards : cfg.N_BOARDS;
}
Int_t ActiveNChannels() {
  return gActiveEpoch ? gActiveEpoch->n_channels : cfg.N_CHANNELS;
}
UShort_t ActiveTimingRefBoard() {
  return gActiveEpoch ? gActiveEpoch->timing_ref_board : cfg.TIMING_REF_BOARD;
}
const std::vector<UShort_t> &ActiveTimingRefBoardChannels() {
  return gActiveEpoch ? gActiveEpoch->timing_ref_board_channels
                      : cfg.TIMING_REF_BOARD_CHANNELS;
}
Bool_t ActiveDoBoardSync() {
  return gActiveEpoch ? gActiveEpoch->do_board_sync : cfg.TIMING_DO_BOARD_SYNC;
}
Bool_t ActiveDoSort() {
  return gActiveEpoch ? gActiveEpoch->do_sort : cfg.TIMING_DO_SORT;
}
Bool_t ActiveHasCathode() {
  return gActiveEpoch ? gActiveEpoch->has_cathode : cfg.HAS_CATHODE;
}
Bool_t ActiveUseSolarisData() {
  return gActiveEpoch ? (gActiveEpoch->source == kSolaris)
                      : cfg.USE_SOLARIS_DATA;
}
Double_t ActiveEventTimeWindowUs() {
  return gActiveEpoch ? gActiveEpoch->event_time_window_us
                      : cfg.EVENT_TIME_WINDOW_US;
}
const TString &ActiveReferenceChannel() {
  return gActiveEpoch ? gActiveEpoch->reference_channel : cfg.REFERENCE_CHANNEL;
}
Double_t ActiveReferenceChannelMinAdc() {
  return gActiveEpoch ? gActiveEpoch->reference_channel_min_adc
                      : cfg.REFERENCE_CHANNEL_MIN_ADC;
}
Double_t ActiveReferenceChannelMaxAdc() {
  return gActiveEpoch ? gActiveEpoch->reference_channel_max_adc
                      : cfg.REFERENCE_CHANNEL_MAX_ADC;
}
DedupStrategy ActiveDedupStrategy() {
  return gActiveEpoch ? gActiveEpoch->dedup_strategy : cfg.DEDUP_STRATEGY;
}
Double_t ActiveStripEMinAdc() {
  return gActiveEpoch ? gActiveEpoch->strip_e_min_adc : cfg.STRIP_E_MIN_ADC;
}
Double_t ActiveStripEMaxAdc() {
  return gActiveEpoch ? gActiveEpoch->strip_e_max_adc : cfg.STRIP_E_MAX_ADC;
}
Double_t ActiveCathodeMaxAdc() {
  return gActiveEpoch ? gActiveEpoch->cathode_max_adc : cfg.CATHODE_MAX_ADC;
}
Double_t ActiveGridMaxAdc() {
  return gActiveEpoch ? gActiveEpoch->grid_max_adc : cfg.GRID_MAX_ADC;
}
Double_t ActiveStrip0MaxAdc() {
  return gActiveEpoch ? gActiveEpoch->strip0_max_adc : cfg.STRIP0_MAX_ADC;
}
Double_t ActiveStrip17MaxAdc() {
  return gActiveEpoch ? gActiveEpoch->strip17_max_adc : cfg.STRIP17_MAX_ADC;
}
Double_t ActiveLeftEvenMaxAdc() {
  return gActiveEpoch ? gActiveEpoch->left_even_max_adc : cfg.LEFT_EVEN_MAX_ADC;
}
Double_t ActiveLeftOddMaxAdc() {
  return gActiveEpoch ? gActiveEpoch->left_odd_max_adc : cfg.LEFT_ODD_MAX_ADC;
}
Double_t ActiveRightEvenMaxAdc() {
  return gActiveEpoch ? gActiveEpoch->right_even_max_adc
                      : cfg.RIGHT_EVEN_MAX_ADC;
}
Double_t ActiveRightOddMaxAdc() {
  return gActiveEpoch ? gActiveEpoch->right_odd_max_adc : cfg.RIGHT_ODD_MAX_ADC;
}

const std::vector<Int_t> &ActiveRunNumbers() {
  return gActiveEpoch ? gActiveEpoch->runs : cfg.RUN_NUMBERS;
}
const TString &ActiveFileTag() {
  static const TString kNone = "";
  return gActiveEpoch ? gActiveEpoch->file_tag : kNone;
}
Int_t ActiveMaxFiles() { return gActiveEpoch ? gActiveEpoch->max_files : -1; }

const std::map<std::pair<Int_t, Int_t>, TString> &ActiveChannelMap() {
  if (gActiveEpoch && !gActiveEpoch->channel_map.empty())
    return gActiveEpoch->channel_map;
  if (!cfg.channelMap64.empty())
    return cfg.channelMap64;
  if (!cfg.channelMap.empty())
    return cfg.channelMap;
  std::cerr << "FATAL: neither channelMap nor channelMap64 is configured."
            << std::endl;
  std::abort();
}

} // namespace Constants
