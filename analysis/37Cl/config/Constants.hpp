#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include "ChannelTiming.hpp"
#include "Paths.hpp"
#include "SlotLayout.hpp"
#include <Rtypes.h>
#include <TString.h>
#include <TSystem.h>
#include <map>
#include <utility>
#include <vector>

namespace Constants {

const std::vector<Int_t> RUN_NUMBERS = {
    37}; // {12, 13, 14, 15, 16, 17, 20, 37};

const TString COMPASS_BASE_DIR = "/home/e-work/LabData/MUSIC/37Cl/";
const Int_t N_FILES = 5;

const TString SIM_BEAM_FILE = "traces_37Cl_beam.root";

const Int_t N_BOARDS = 4;
const Int_t N_CHANNELS = 16;

const UShort_t TIMING_REF_BOARD = 1;
const std::vector<UShort_t> TIMING_REF_BOARD_CHANNELS = {12, 0, 0, 0};
const Double_t TIMING_MIN_ENERGY = 300;
const Double_t TIMING_MAX_ENERGY = 1500;
const Double_t TIMING_OVERLAP_MARGIN_S = 1.0;
const Double_t TIMING_THRESH_DT_US = 175.0;
const Double_t TIMING_MAX_ABS_SHIFT_S = 1.5;
const Double_t TIMING_MIN_NSD2_GAIN = 1.10;

const Int_t MAX_TRACE_SAVES = 10;

const Bool_t REJECT_FLAGGED_EVENTS = kTRUE;

// Post EventBuilder + CalibrateBeam only!
const Bool_t IGNORE_SHORT_STRIPS = kFALSE;
const Bool_t IGNORE_STRIP_0 = kFALSE;
const Bool_t IGNORE_STRIP_17 = kFALSE;

const Bool_t SKIP_EXISTING = kTRUE;
const Bool_t RUN_TRACES = kTRUE;
const Bool_t SAVE_PLOTS = kTRUE;

const Bool_t SKIP_CALIBRATION = kFALSE;

const Int_t MAX_FUSED_WORKERS = 16;
const Int_t MAX_DIAGNOSE_WORKERS = 32;

// 37Cl recorded grid events
const Bool_t USE_NEAREST_TO_GRID = kTRUE;
const Bool_t USE_TIME_WINDOW_EVENTS = kFALSE;
const Double_t EVENT_TIME_WINDOW_US = 8.0;

const Bool_t USE_GPU_ACCELERATION = kTRUE;
const Int_t MAX_GPU_CONCURRENT_SORTS = 5;

const Double_t STRIP_E_MIN_MEV = -0.2;
const Double_t STRIP_E_MAX_MEV = 10.0;
const Double_t CATHODE_E_MAX_MEV = 300;
const Double_t TOTAL_E_MIN_MEV = 10.0;
const Double_t TOTAL_E_MAX_MEV = 150.0;

const Double_t NORM_MUSIC_MEV = 2.30;

const Int_t REACTION_STRIP_MIN = 7;
const Int_t REACTION_STRIP_MAX = 15;

const Double_t STRIP_E_MIN_ADC = 0.0;
const Double_t STRIP_E_MAX_ADC = 5500.0;
const Double_t TOTAL_E_MIN_ADC = 0.0;
const Double_t TOTAL_E_MAX_ADC = 60000.0;

const std::map<std::pair<Int_t, Int_t>, TString> channelMap = {
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

// Per-channel TTF timing offsets (ps). Board 2 was misconfigured with a 992 ns
// TTF delay vs the 896 ns baseline on every other board, so every named board-2
// channel needs a +96 ns (96000 ps) correction subtracted from its raw
// timestamp. See ChannelTiming.hpp.
const std::map<std::pair<Int_t, Int_t>, Long64_t> ttfOffsetPs = {
    {{2, 0}, 96000},  {{2, 2}, 96000}, {{2, 4}, 96000},  {{2, 5}, 96000},
    {{2, 6}, 96000},  {{2, 8}, 96000}, {{2, 10}, 96000}, {{2, 12}, 96000},
    {{2, 13}, 96000}, {{2, 15}, 96000}};

} // namespace Constants

#endif
