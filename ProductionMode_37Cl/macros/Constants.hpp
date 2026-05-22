#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <map>
#include <utility>
#include <vector>

namespace Paths {
inline TString ProjectRootOf(const char *file) {
  TString path = file;
  path.ReplaceAll("/./", "/");
  TString macros_dir = gSystem->DirName(path);
  return gSystem->DirName(macros_dir);
}
} // namespace Paths

namespace Constants {

inline const std::vector<Int_t> RUN_NUMBERS = {16};

const Int_t N_FILES = -1;

inline const TString COMPASS_BASE_DIR = "/home/e-work/LabData/MUSIC/37Cl/";

const Int_t N_BOARDS = 4;
const Int_t N_CHANNELS = 16;

const UShort_t REF_BOARD = 1;
inline const std::vector<UShort_t> BOARD_CHANNELS = {12, 0, 0, 0};

const Double_t MIN_ENERGY = 300;
const Double_t MAX_ENERGY = 1500;
const Double_t OVERLAP_MARGIN_S = 1.0;
const Double_t THRESH_DT_US = 175.0;
const Double_t TIMING_MAX_ABS_SHIFT_S = 1.5;

const Int_t MAX_TRACE_SAVES = 10;

const Bool_t REJECT_FLAGGED_EVENTS = kFALSE;

// Post EventBuilder only!
const Bool_t IGNORE_SHORT_STRIPS = kFALSE;
const Bool_t IGNORE_STRIP_0 = kFALSE;
const Bool_t IGNORE_STRIP_17 = kFALSE;

const Bool_t SKIP_EXISTING = kTRUE;
const Bool_t RUN_TRACES = kTRUE;

const Bool_t BEAM_CAL_REDRAW_CUTS = kFALSE;

const Bool_t SAVE_PLOTS = kTRUE;

const Int_t MAX_FUSED_WORKERS = 16;
const Int_t MAX_DIAGNOSE_WORKERS = 32;

const Bool_t USE_NEAREST_TO_GRID = kFALSE;
const Bool_t SAVE_PER_CHANNEL_TIMESTAMPS_FLAGS = kFALSE;

const Bool_t USE_GPU_ACCELERATION = kTRUE;
const Int_t MAX_GPU_CONCURRENT_SORTS = 5;

inline Char_t LongAnodeSide(Int_t strip) {
  if (strip < 1 || strip > 16)
    return ' ';
  return (strip % 2 == 1) ? 'L' : 'R';
}

inline Bool_t IsLongSide(Int_t strip, Char_t side) {
  return strip >= 1 && strip <= 16 && side == LongAnodeSide(strip);
}

inline Bool_t IsShortSide(Int_t strip, Char_t side) {
  return strip >= 1 && strip <= 16 && (side == 'L' || side == 'R') &&
         side != LongAnodeSide(strip);
}

inline Int_t FirstStrip() { return IGNORE_STRIP_0 ? 1 : 0; }
inline Int_t LastStrip() { return IGNORE_STRIP_17 ? 16 : 17; }
inline Int_t NumActiveStrips() { return LastStrip() - FirstStrip() + 1; }
inline Bool_t IsActiveStrip(Int_t strip) {
  return strip >= FirstStrip() && strip <= LastStrip();
}

const Long64_t BASELINE_MIN_ENTRIES = 500;
const Int_t BASELINE_HIST_BINS = 4096;
const Double_t BASELINE_HIST_MAX_ADC = 16384.0;

const Long64_t REGION_TRACES_MAX_INDIV = 10000;

const Double_t STRIP_E_MIN_MEV = -0.2;
const Double_t STRIP_E_MAX_MEV = 5.0;
const Double_t TOTAL_E_MIN_MEV = 10;
const Double_t TOTAL_E_MAX_MEV = 60.0;

const Double_t STRIP_E_MIN_ADC = 0.0;
const Double_t STRIP_E_MAX_ADC = 5500.0;
const Double_t TOTAL_E_MIN_ADC = 0.0;
const Double_t TOTAL_E_MAX_ADC = 60000.0;

inline const std::map<std::pair<Int_t, Int_t>, TString> channelMap = {
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

// Per-channel CoMPASS TTF (trapezoidal trigger filter) delay relative to the
// 896 ns baseline used by most channels. Subtract this from a channel's raw
// Timestamp to put every channel on the same timing baseline. The board-level
// beam-pattern shift is computed at second-scale resolution and does not see
// these ns-scale offsets, so this correction is independent and complementary.
// Pattern (from CoMPASS settings.xml across all runs): every named board-2
// channel runs at TTF=992 ns; every other named channel runs at TTF=896 ns.
inline const std::map<std::pair<Int_t, Int_t>, Long64_t> ttfOffsetPs = {
    {{2, 0}, 96000},  {{2, 2}, 96000}, {{2, 4}, 96000},  {{2, 5}, 96000},
    {{2, 6}, 96000},  {{2, 8}, 96000}, {{2, 10}, 96000}, {{2, 12}, 96000},
    {{2, 13}, 96000}, {{2, 15}, 96000}};

inline Long64_t LookupTTFOffsetPs(Int_t board, Int_t channel) {
  std::map<std::pair<Int_t, Int_t>, Long64_t>::const_iterator it =
      ttfOffsetPs.find(std::pair<Int_t, Int_t>(board, channel));
  return (it == ttfOffsetPs.end()) ? 0 : it->second;
}
} // namespace Constants

#endif
