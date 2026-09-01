#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include "DedupStrategy.hpp"
#include "RunEpoch.hpp"
#include "SlotLayout.hpp"
#include <Rtypes.h>
#include <RtypesCore.h>
#include <TString.h>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

struct StripSumScatterConfig {
  enum PureBeamGate { PURE_BEAM_GATE_S0_S1, PURE_BEAM_GATE_S1_S2 };
  PureBeamGate PURE_BEAM_GATE;

  Int_t REACTION_STRIP_MIN;
  Int_t REACTION_STRIP_MAX;

  Int_t REQUIRE_SMOOTHNESS_END_STRIP;
  Double_t REQUIRE_SMOOTHNESS_MAX_STEP;

  Double_t REAC_JUMP_MIN;
  Double_t REAC_JUMP_MAX;
  Double_t END_STRIP_MAX;

  Double_t PILEUP_THRESHOLD;
  Double_t NOISE_THRESHOLD;

  Bool_t REJECT_NOISE;
  Double_t NOISE_THRESH_PY;
  Int_t NOISE_MIN_STRIPS;

  Bool_t REJECT_PILEUP;
  Double_t PILEUP_THRESH_PY;
  Int_t PILEUP_MIN_STRIPS;

  // Both-ends multiplicity cut: reject an event when more than MAX strips in
  // 1..COUNT_TO had BOTH ends fire. Read off raw ADC, so it is independent of
  // IGNORE_SHORT_STRIPS. Charge sharing means a displaced track lights the
  // short end on every strip of one parity at once, and tagging on a summed
  // trace then enriches for those -- this removes them. MAX < 0 disables it.
  // COUNT_TO < 16 restricts the count to strips upstream of the reaction, where
  // the sim has no sharing at all, without penalising products downstream.
  Int_t BOTH_MULT_MAX;
  Int_t BOTH_MULT_COUNT_TO;

  Bool_t REJECT_OFFBEAM;
  Double_t OFFBEAM_DIST;
  Int_t OFFBEAM_MIN_STRIPS;

  Double_t TRIGGER_NSIGMA;
  Double_t TRIGGER_CFD_FRAC;
  Int_t PLATEAU_POST;
  Int_t CLUSTER_SMOOTH_WINDOW;
  Int_t SEED_HALF_BINS;

  Int_t SAVGOL_HALF;

  Int_t TRACES_PER_CLASS;

  Int_t X_LO;
  Int_t X_HI;

  Int_t GATE_STRIP_X;
  Int_t GATE_STRIP_Y;
  Double_t GATE_NSIGMA_X;
  Double_t GATE_NSIGMA_Y;
  Double_t GATE_MIN;
  Double_t GATE_MAX;
  Int_t GATE_BINS;

  Double_t XMIN;
  Double_t XMAX;
  Int_t XBINS;
  Double_t YMIN;
  Double_t YMAX;
  Int_t YBINS;

  std::map<Int_t, std::pair<Double_t, Double_t>> Y_RANGE;

  Long64_t SAMPLE_MAX_POINTS;

  Bool_t RERUN_SIM;
  Int_t CANDIDATE_REAC_STRIP;

  Bool_t REQUIRE_SMOOTHNESS;
  Bool_t REQUIRE_GATE_S3_S4;
  Bool_t REQUIRE_GATE_S5_S6;
  Bool_t SKIP_SAVGOL_PLOTS;
  Bool_t REQUIRE_STRIP_16_BELOW_BEAM;

  // Also render the selected region traces under the OTHER decode
  // (long-side-only <-> L+R sum) from the same calibration and the
  // same selected events, into a separate plot plus a text dump.
  Bool_t ALT_DECODE_REGION_TRACES;

  void SetDefaults();
};

class DatasetConfig {
public:
  // Data source
  Bool_t USE_SOLARIS_DATA;
  TString SOL_BASE_DIR;
  TString SOL_SPLIT_DIR;
  Double_t SOL_SPLIT_CHUNK_SECONDS;
  Int_t SOL_N_SPLIT_WORKERS;

  TString COMPASS_BASE_DIR;
  // Flat run list, used when EPOCHS is empty. When EPOCHS is populated the
  // pipeline walks the epochs instead and this is ignored.
  std::vector<Int_t> RUN_NUMBERS;
  std::vector<RunEpoch> EPOCHS;
  Int_t N_CHUNKS;

  TString SIM_BEAM_FILE;

  // Hardware layout
  Int_t N_BOARDS;
  Int_t N_CHANNELS;
  UShort_t TIMING_REF_BOARD;
  std::vector<UShort_t> TIMING_REF_BOARD_CHANNELS;

  // Timing
  Double_t TIMING_MIN_ENERGY;
  Double_t TIMING_MAX_ENERGY;
  Double_t TIMING_OVERLAP_MARGIN_S;
  Double_t TIMING_THRESH_DT_US;
  Double_t TIMING_MAX_ABS_SHIFT_S;
  Double_t TIMING_SHIFT_COARSE_STEP_US;
  Double_t TIMING_SHIFT_FINE_STEP_US;
  Double_t TIMING_SHIFT_FINE_HALF_WIDTH_US;
  Int_t TIMING_SHIFT_MIN_NPTS;
  Int_t TIMING_SHIFT_MAX_SCAN_CANDIDATES;
  Bool_t TIMING_DO_BOARD_SYNC;
  Bool_t TIMING_DO_SORT;

  Bool_t REJECT_FLAGGED_EVENTS;

  Bool_t IGNORE_SHORT_STRIPS;
  Bool_t IGNORE_STRIP_0;
  Bool_t IGNORE_STRIP_17;

  Bool_t HAS_CATHODE;
  Bool_t HAS_GRID;
  Bool_t HAS_STRIP0;
  Bool_t HAS_STRIP17;

  Bool_t SKIP_EXISTING;
  Bool_t SAVE_PLOTS;

  Bool_t SKIP_CALIBRATION;

  // Number of sample traces to save during event build and normed summary
  // passes (0 = disabled). Saved as overlays in events_summary and
  // events_summary_normed alongside the histograms.
  Int_t SAVE_SAMPLE_TRACES;

  Int_t MAX_FUSED_WORKERS;

  TString REFERENCE_CHANNEL;
  Double_t EVENT_TIME_WINDOW_US;
  DedupStrategy DEDUP_STRATEGY;

  Bool_t USE_GPU_ACCELERATION;
  Int_t MAX_GPU_CONCURRENT_SORTS;

  Double_t STRIP_DE_OVERVIEW_MIN_NORMED;
  Double_t STRIP_DE_OVERVIEW_MAX_NORMED;
  Double_t STRIP_DE_MIN_NORMED;
  Double_t STRIP_DE_MAX_NORMED;
  Double_t CATHODE_E_MAX_NORMED;
  Double_t TOTAL_E_MIN_NORMED;
  Double_t TOTAL_E_MAX_NORMED;

  StripSumScatterConfig STRIP_SUM_SCATTER_CONFIG;

  Double_t STRIP_E_MIN_ADC;
  Double_t STRIP_E_MAX_ADC;
  Double_t TOTAL_E_MIN_ADC;
  Double_t TOTAL_E_MAX_ADC;

  Double_t GRID_MIN_ADC;
  Double_t GRID_MAX_ADC;

  Double_t REFERENCE_CHANNEL_MIN_ADC;
  Double_t REFERENCE_CHANNEL_MAX_ADC;

  Double_t STRIP0_MAX_ADC;
  Double_t STRIP17_MAX_ADC;
  Double_t CATHODE_MAX_ADC;

  Double_t LEFT_EVEN_MAX_ADC;
  Double_t LEFT_ODD_MAX_ADC;
  Double_t RIGHT_EVEN_MAX_ADC;
  Double_t RIGHT_ODD_MAX_ADC;

  std::map<std::pair<Int_t, Int_t>, TString> channelMap;
  std::map<std::pair<Int_t, Int_t>, TString> channelMap64;
  std::map<std::pair<Int_t, Int_t>, Long64_t> ttfOffsetPs;

  DatasetConfig();
};

namespace Constants {
extern const DatasetConfig &cfg;

// Epoch selection. Pipeline sets the active epoch around each epoch's work;
// every Active*() below reads it when set and falls back to the flat cfg block
// otherwise, so a dataset that declares no epochs behaves exactly as before.
void SetActiveEpoch(const RunEpoch *epoch);
const RunEpoch *GetActiveEpoch();
// Epoch owning a run number, or null when no UNTAGGED epoch declares it. A
// tagged epoch reuses another era's run numbers and is addressed by tag, never
// by number alone.
const RunEpoch *EpochForRun(Int_t run);
// Output-name prefix of the active epoch ("" when none, so existing filenames
// are unchanged).
const TString &ActiveFileTag();

Int_t ActiveNBoards();
Int_t ActiveNChannels();
UShort_t ActiveTimingRefBoard();
const std::vector<UShort_t> &ActiveTimingRefBoardChannels();
Bool_t ActiveDoBoardSync();
Bool_t ActiveDoSort();
Bool_t ActiveHasCathode();
Bool_t ActiveUseSolarisData();
Double_t ActiveEventTimeWindowUs();
const TString &ActiveReferenceChannel();
Double_t ActiveReferenceChannelMinAdc();
Double_t ActiveReferenceChannelMaxAdc();
DedupStrategy ActiveDedupStrategy();
Double_t ActiveStripEMinAdc();
Double_t ActiveStripEMaxAdc();
Double_t ActiveCathodeMaxAdc();
Double_t ActiveGridMaxAdc();
Double_t ActiveStrip0MaxAdc();
Double_t ActiveStrip17MaxAdc();
Double_t ActiveLeftEvenMaxAdc();
Double_t ActiveLeftOddMaxAdc();
Double_t ActiveRightEvenMaxAdc();
Double_t ActiveRightOddMaxAdc();

// Runs of the active epoch, else the flat RUN_NUMBERS. ActiveMaxFiles caps
// subfiles per run for the active epoch (-1 = all, and -1 without an epoch).
const std::vector<Int_t> &ActiveRunNumbers();
Int_t ActiveMaxFiles();

// Returns the active epoch's channel map when an epoch is set, else
// channelMap64 if populated, else channelMap. Aborts if none is configured.
const std::map<std::pair<Int_t, Int_t>, TString> &ActiveChannelMap();
} // namespace Constants

#endif
