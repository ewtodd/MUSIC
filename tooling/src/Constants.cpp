#include "Constants.hpp"

Int_t TargetGasA(TargetGas gas) {
  switch (gas) {
  case kHELIUM:
    return 4;
  }
  return 0;
}

Double_t TargetGasAtomsPerMolecule(TargetGas gas) {
  switch (gas) {
  case kHELIUM:
    return 1.0;
  }
  return 0.0;
}

void CrossSectionConfig::SetDefaults() {
  // Deliberately not a working experiment: a cross-section dataset must
  // state its gas and beam; zero pressure makes a forgotten one fail loudly.
  TARGET_GAS = kHELIUM;
  GAS_PRESSURE_TORR = 0.0;
  GAS_PRESSURE_TORR_ERR = 0.0;

  BEAM_A = 0;
  BEAM_Z = 0;
  BEAM_ELEMENT = "";

  TALYS_MODELS.clear();
  // The shared TALYS settings, from the group's 14O(a,p) input; not here:
  // the energies (per file), alphaomp (per model) and the rate options.
  TALYS_COMMON_KEYWORDS = {
      "# level density",
      "ldmodel 1",
      "maxlevelstar 30",
      "# enhanced accuracy",
      "transpower 20",
      "ecisstep 0.02",
      "xseps 1.e-30",
      "transeps 1.e-30",
      "popeps 1.e-30",
      "elow 1.e-06",
      "# width fluctuation correction",
      "widthfluc 20",
      "widthmode 1",
      "# pre-equilibrium off",
      "preequilibrium n",
      "# equidistant excitation-energy binning, enhanced segments",
      "bins 80",
      "equidistant y",
      "segment 4",
  };

  BEAM_SIM_FILE = "";

  XS_STRIP_MIN = 3;
  XS_STRIP_MAX = 15;
  FELDMAN_COUSINS_MAX_COUNT = 50;
  EFFECTIVE_ENERGY = kTRUE;
  PRELIMINARY = kFALSE;
  EPOCHS.clear();

  CHANNELS.clear();
}

void StripSumScatterConfig::SetDefaults() {
  PURE_BEAM_GATE = PURE_BEAM_GATE_S0_S1;
  TAG_RESOLVE = TAG_RESOLVE_FIRST_STRIP;
  PURE_BEAM_NSIGMA = 3.5;

  POST_TRIGGER_SUM_STRIPS = 3;
  POST_WINDOW_LAST_STRIP = 17;
  MAX_STRIP_SUM_WORKERS = 12;

  REACTION_STRIP_MIN = 2;
  REACTION_STRIP_MAX = 15;

  SMOOTHNESS_NSIGMA = 2.5; // the 87Rb macros' 1.2/12 at that noise
  TAIL_FALL_FROM_STRIP = 0;
  TAIL_RISE_NSIGMA = 0.0;
  TAIL_RETURN_NSIGMA = 0.0;
  TAIL_RERISE_NSIGMA = 0.0;
  POST_ABOVE_NSIGMA = 0.0;
  POST_ABOVE_STRIPS = 0;
  TAIL_CLIFF_MAX_FRACTION = 0.0;
  POST_CROSS_MIN_STRIP = 0;
  CUT_VARIATION = kTRUE;
  CUT_VARIATION_NSIGMA_STEP = 3.0; // the resolution, as the 87Rb paper
  CUT_VARIATION_CLIFF_STEP = 0.1;
  Y_RATIO_TO_UPSTREAM = kTRUE;

  REAC_JUMP_NSIGMA = 1.0;
  END_STRIP_NSIGMA = 0.0; // below the beam mean

  PILEUP_NSIGMA = 8.0; // ~1.3 at the 87Rb spread, the old absolute value
  PILEUP_MIN_STRIPS = 1;
  NOISE_NSIGMA = 12.0; // ~0.5 at the 87Rb spread, the old absolute value
  NOISE_MIN_STRIPS = 1;

  REGION_CUT_REDRAW = kFALSE;
  AN_REGION_NSIGMA = 2.0;
  AA_REGION_NSIGMA = 2.0;
  AN_REGION_MODE = AN_REGION_MIXTURE;
  REQUIRE_BEAM_UPSTREAM_OF_REAC = kFALSE;
  BEAM_UPSTREAM_NSIGMA = 2.0;
  BOTH_MULT_MAX = -1;      // disabled by default
  BOTH_MULT_COUNT_TO = 16; // whole trace unless narrowed

  SEED_HALF_BINS = 40;

  TRACES_PER_CLASS = 40;

  X_LO = 1;
  X_HI = 16;

  GATE_STRIP = 1;
  GATE_NSIGMA = 3.5;
  GATE_CENTER = 0.0;
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

  SKIP_SAVGOL_PLOTS = kFALSE;
  PLOT_REGION_MEAN_TRACES = kFALSE;
  PLOT_ADC_TRACES = kFALSE;
  REQUIRE_STRIP_16_BELOW_BEAM = kFALSE;
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

  PULSE_HISTORY_CORRECTION = kFALSE;
  PULSE_HISTORY_MIN_EVENTS = 20000;
  PULSE_HISTORY_BEAM_LO = 0.8;
  PULSE_HISTORY_BEAM_HI = 1.3;
  PULSE_HISTORY_OWN_LO = 0.5;
  PULSE_HISTORY_OWN_HI = 1.6;
  PULSE_HISTORY_APPLY_MAX_US = 316.0;
  PULSE_HISTORY_AMP_BINS = 1;
  PULSE_HISTORY_TRAP_RISE_US = 2.0;
  PULSE_HISTORY_TRAP_FLAT_US = 3.0;
  PULSE_HISTORY_PEAKING_US = 1.5;

  IGNORE_SHORT_STRIPS = kFALSE;
  IGNORE_STRIP_0 = kFALSE;
  IGNORE_STRIP_17 = kFALSE;
  REQUIRE_STRIP_0 = kTRUE;

  HAS_CATHODE = kTRUE;
  HAS_GRID = kTRUE;
  HAS_STRIP0 = kTRUE;
  HAS_STRIP17 = kTRUE;

  SKIP_EXISTING = kTRUE;
  SAVE_FULL_PLOTS = kFALSE;
  PLOT_SAMPLE_FILES = 2;

  SAVE_SAMPLE_TRACES = 0;

  MAX_FUSED_WORKERS = 12;

  REFERENCE_CHANNEL = "Grid";
  EVENT_TIME_WINDOW_US = 8.0;
  SEED_HOLDOFF_US = 0.0;
  SEED_HOLDOFF_MAX_RATIO = 0.5;
  DEDUP_STRATEGY = kLARGEST_ENERGY;

  USE_GPU_ACCELERATION = kTRUE;
  MAX_GPU_CONCURRENT_SORTS = 20;

  BEAM_GATE_NSIGMA = 3.0;
  STRIP_DE_MIN_NORMED = 0.8;
  STRIP_DE_MAX_NORMED = 1.3;

  STRIP_SUM_SCATTER_CONFIG.SetDefaults();
  CROSS_SECTION_CONFIG.SetDefaults();

  STRIP_E_MIN_ADC = 0.0;
  STRIP_E_MAX_ADC = 4096.0;

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

/// The epoch currently being processed, or null outside epoch work (and for
/// datasets that declare none). Every accessor below falls back to the flat cfg
/// block when it is null, so behaviour without epochs is unchanged.
static const RunEpoch *gActiveEpoch = nullptr;

void SetActiveEpoch(const RunEpoch *epoch) { gActiveEpoch = epoch; }
const RunEpoch *GetActiveEpoch() { return gActiveEpoch; }

const RunEpoch *EpochForRun(Int_t run) {
  // Only untagged epochs are addressable by run number; a tagged one's runs
  // collide with another era's (a bare number is a coin flip), so use its tag.
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

RunEpoch MakeEpoch(const TString &name, const std::vector<Int_t> &runs) {
  RunEpoch ep;
  ep.name = name;
  ep.runs = runs;
  ep.source = cfg.USE_SOLARIS_DATA ? kSolaris : kCoMPASS;
  ep.n_boards = cfg.N_BOARDS;
  ep.n_channels = cfg.N_CHANNELS;
  ep.channel_map =
      !cfg.channelMap64.empty() ? cfg.channelMap64 : cfg.channelMap;
  // An epoch declared before the channel map was set would carry the
  // tooling defaults and fail far downstream, so refuse it here.
  if (ep.channel_map.empty()) {
    std::cerr << "FATAL: MakeEpoch(\"" << name
              << "\") called before the channel map was set; declare epochs "
                 "at the end of InitDatasetConfig."
              << std::endl;
    std::abort();
  }
  for (std::map<std::pair<Int_t, Int_t>, TString>::const_iterator it =
           ep.channel_map.begin();
       it != ep.channel_map.end(); ++it)
    if (it->first.first >= ep.n_boards || it->first.second >= ep.n_channels) {
      std::cerr << "FATAL: MakeEpoch(\"" << name << "\"): channel map entry ("
                << it->first.first << ", " << it->first.second
                << ") lies outside N_BOARDS x N_CHANNELS = " << ep.n_boards
                << " x " << ep.n_channels
                << "; set those before declaring the epoch." << std::endl;
      std::abort();
    }
  ep.timing_ref_board = cfg.TIMING_REF_BOARD;
  ep.timing_ref_board_channels = cfg.TIMING_REF_BOARD_CHANNELS;
  ep.do_board_sync = cfg.TIMING_DO_BOARD_SYNC;
  ep.do_sort = cfg.TIMING_DO_SORT;
  ep.event_time_window_us = cfg.EVENT_TIME_WINDOW_US;
  ep.seed_holdoff_us = cfg.SEED_HOLDOFF_US;
  ep.seed_holdoff_max_ratio = cfg.SEED_HOLDOFF_MAX_RATIO;
  ep.reference_channel = cfg.REFERENCE_CHANNEL;
  ep.reference_channel_min_adc = cfg.REFERENCE_CHANNEL_MIN_ADC;
  ep.reference_channel_max_adc = cfg.REFERENCE_CHANNEL_MAX_ADC;
  ep.dedup_strategy = cfg.DEDUP_STRATEGY;
  ep.has_cathode = cfg.HAS_CATHODE;
  ep.split_chunk_seconds = cfg.SOL_SPLIT_CHUNK_SECONDS;
  ep.pulse_history = cfg.PULSE_HISTORY_GROUPS;
  ep.ignore_short_strips = cfg.IGNORE_SHORT_STRIPS;
  ep.strip_e_min_adc = cfg.STRIP_E_MIN_ADC;
  ep.strip_e_max_adc = cfg.STRIP_E_MAX_ADC;
  ep.cathode_max_adc = cfg.CATHODE_MAX_ADC;
  ep.grid_max_adc = cfg.GRID_MAX_ADC;
  ep.strip0_max_adc = cfg.STRIP0_MAX_ADC;
  ep.strip17_max_adc = cfg.STRIP17_MAX_ADC;
  ep.left_even_max_adc = cfg.LEFT_EVEN_MAX_ADC;
  ep.left_odd_max_adc = cfg.LEFT_ODD_MAX_ADC;
  ep.right_even_max_adc = cfg.RIGHT_EVEN_MAX_ADC;
  ep.right_odd_max_adc = cfg.RIGHT_ODD_MAX_ADC;
  return ep;
}

std::vector<Int_t> RunRange(Int_t first, Int_t last) {
  std::vector<Int_t> runs;
  for (Int_t r = first; r <= last; r++)
    runs.push_back(r);
  return runs;
}

/// Runs for a binary that runs with no active epoch: the enabled epochs named
/// in CROSS_SECTION_CONFIG.EPOCHS, or all enabled epochs when that is empty.
/// Computed once; the config does not change after static init.
static const std::vector<Int_t> &EpochUnionRuns() {
  static std::vector<Int_t> runs;
  static Bool_t built = kFALSE;
  if (built)
    return runs;
  built = kTRUE;
  const std::vector<TString> &want = cfg.CROSS_SECTION_CONFIG.EPOCHS;
  for (Int_t w = 0; w < Int_t(want.size()); w++) {
    Bool_t found = kFALSE;
    for (Int_t e = 0; e < Int_t(cfg.EPOCHS.size()) && !found; e++)
      found = cfg.EPOCHS[e].name == want[w];
    if (!found)
      std::cerr << "WARNING: CROSS_SECTION_CONFIG.EPOCHS names \"" << want[w]
                << "\" but no such epoch is declared." << std::endl;
  }
  for (Int_t e = 0; e < Int_t(cfg.EPOCHS.size()); e++) {
    const RunEpoch &ep = cfg.EPOCHS[e];
    if (!ep.enabled)
      continue;
    Bool_t wanted = want.empty();
    for (Int_t w = 0; w < Int_t(want.size()) && !wanted; w++)
      wanted = ep.name == want[w];
    if (wanted)
      runs.insert(runs.end(), ep.runs.begin(), ep.runs.end());
  }
  return runs;
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
Double_t ActiveSeedHoldoffUs() {
  return gActiveEpoch ? gActiveEpoch->seed_holdoff_us : cfg.SEED_HOLDOFF_US;
}
Double_t ActiveSeedHoldoffMaxRatio() {
  return gActiveEpoch ? gActiveEpoch->seed_holdoff_max_ratio
                      : cfg.SEED_HOLDOFF_MAX_RATIO;
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
  if (gActiveEpoch)
    return gActiveEpoch->runs;
  if (cfg.RUN_NUMBERS.empty() && !cfg.EPOCHS.empty())
    return EpochUnionRuns();
  return cfg.RUN_NUMBERS;
}
const TString &ActiveFileTag() {
  static const TString kNone = "";
  return gActiveEpoch ? gActiveEpoch->file_tag : kNone;
}
Int_t ActiveMaxFiles() { return gActiveEpoch ? gActiveEpoch->max_files : -1; }
Double_t ActiveSplitChunkSeconds() {
  return gActiveEpoch ? gActiveEpoch->split_chunk_seconds
                      : cfg.SOL_SPLIT_CHUNK_SECONDS;
}
const PulseHistoryGroups &ActivePulseHistoryGroups() {
  return gActiveEpoch ? gActiveEpoch->pulse_history : cfg.PULSE_HISTORY_GROUPS;
}

static const RunEpoch *AnalysisEpoch() {
  const std::vector<TString> &want = cfg.CROSS_SECTION_CONFIG.EPOCHS;
  const RunEpoch *found = nullptr;
  Int_t n = 0;
  for (Int_t e = 0; e < Int_t(cfg.EPOCHS.size()); e++) {
    const RunEpoch &ep = cfg.EPOCHS[e];
    Bool_t take = kFALSE;
    if (want.empty())
      take = ep.enabled;
    else
      for (Int_t w = 0; w < Int_t(want.size()) && !take; w++)
        take = ep.name == want[w];
    if (take) {
      found = &ep;
      n++;
    }
  }
  return n == 1 ? found : nullptr;
}

void ActivateAnalysisEpoch() {
  const RunEpoch *ep = AnalysisEpoch();
  if (!ep)
    return;
  SetActiveEpoch(ep);
  std::cout << "epoch " << ep->name << " ("
            << (ep->source == kSolaris ? "SOLARIS" : "CoMPASS") << ", "
            << ep->runs.size() << " run(s)"
            << (ep->file_tag.Length() ? ", files tagged " + ep->file_tag : "")
            << ")" << std::endl;
}

Bool_t ActiveIgnoreShortStrips() {
  if (gActiveEpoch)
    return gActiveEpoch->ignore_short_strips;
  if (const RunEpoch *ep = AnalysisEpoch())
    return ep->ignore_short_strips;
  return cfg.IGNORE_SHORT_STRIPS;
}

// The plot sample is a per-thread mark, set by every file-walking tool
// before each file, so the stages need no plumbing to know when to draw.
namespace {
thread_local Bool_t t_plot_sample = kFALSE;
}
Bool_t FileInSample() { return cfg.SAVE_FULL_PLOTS || t_plot_sample; }
Bool_t SavePlots() { return FileInSample(); }
namespace {
// A stream with no buffer discards everything written to it; one per thread
// so the failed-write state it keeps is never shared.
std::ostream &NullStream() {
  static thread_local std::ostream null(nullptr);
  return null;
}
} // namespace
std::ostream &Detail() { return FileInSample() ? std::cout : NullStream(); }
std::ostream &DetailErr() { return FileInSample() ? std::cerr : NullStream(); }
Bool_t InPlotSample(Int_t k) {
  return cfg.SAVE_FULL_PLOTS || k < cfg.PLOT_SAMPLE_FILES;
}
void SetPlotsThisFile(Bool_t on) { t_plot_sample = on; }

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
