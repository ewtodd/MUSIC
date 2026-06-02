#ifndef SI_CALIB_CONSTANTS_HPP
#define SI_CALIB_CONSTANTS_HPP

#include "Constants.hpp"
#include <Rtypes.h>
#include <TString.h>
#include <map>
#include <vector>

// Per-dataset constants for the Si-detector calibration / stopping-power
// analysis. This is dataset config (a data blob, like Constants.hpp); it lives
// on the -I analysis/<DATASET>/config path and is resolved per dataset exactly
// as Constants.hpp is. Datasets without a Si-calibration run set simply omit
// this file -- the Makefile then drops the Si binaries from the build.
//
// Kept in namespace SiCalib so it does not collide with namespace Constants.
namespace SiCalib {

const Bool_t INCLUDE_PPAC = kFALSE;

// Energy (MeV) that the beam has when it reaches the Si detector for each
// zero-pressure anchor run.  PPAC_OUT = PPAC not in beam; PPAC_IN = PPAC
// (4x1.5 um mylar) in beam.  All PPAC_IN values are model averages over models
// 1-4 (Ziegler, ATIMA12LS, ATIMA12NoLS, ATIMA14Weick) from
// StoppingPowerDiagnostic.
//
// Run 21:  TOF = 91.87 MeV, no Ti windows.
const Double_t E_RUN_21_PPAC_OUT = 91.87;
const Double_t E_RUN_21_PPAC_IN = E_RUN_21_PPAC_OUT - 18.1162;
// Run 22:  TOF = 71.35 MeV, no Ti windows.
const Double_t E_RUN_22_PPAC_OUT = 71.35;
const Double_t E_RUN_22_PPAC_IN = E_RUN_22_PPAC_OUT - 19.5146;
// Run 23:  TOF = 71.35 MeV, front Ti window only (0.9 mg/cm^2).
const Double_t E_RUN_23_PPAC_OUT = 58.6312;
const Double_t E_RUN_23_PPAC_IN = 38.3174;
// Run 24:  TOF = 92.00 MeV, front Ti window only (0.9 mg/cm^2).
const Double_t E_RUN_24_PPAC_OUT = 80.0511;
const Double_t E_RUN_24_PPAC_IN = 61.2780;
// Run 25:  TOF = 92.00 MeV, both Ti windows (0.9 + 1.3 mg/cm^2),
const Double_t E_RUN_25_PPAC_OUT = 62.5413;
const Double_t E_RUN_25_PPAC_IN = 42.3088;

inline const std::vector<Int_t> RUN_NUMBERS = {21, 22, 23, 24, 25, 26, 27, 28,
                                               29, 30, 31, 32, 33, 34, 35, 36};

// Absolute path to a Si run's ROOT file, under this dataset's CoMPASS run data
// directory (shared with the main analysis config,
// Constants::COMPASS_BASE_DIR).
inline TString GetRunFilepath(Int_t run) {
  return Form("%s/run_%d/ROOT/DataR_run_%d.root",
              Constants::COMPASS_BASE_DIR.Data(), run, run);
}

const UShort_t SI_DETECTOR_CHANNEL = 11;

const Double_t ADC_MIN = 0.0;
const Double_t ADC_MAX = 16384.0;
const Double_t BIN_WIDTH_ADC = 25.0;
const Int_t N_ADC_BINS = (Int_t)((ADC_MAX - ADC_MIN) / BIN_WIDTH_ADC);

inline const std::map<Int_t, Double_t> MU_GUESSES = {
    {21, 15300}, {22, 10000}, {23, 9500}, {24, 13200}, {25, 10150}, {26, 9200},
    {27, 8500},  {28, 7800},  {29, 6800}, {30, 6000},  {31, 5100},  {32, 4200},
    {33, 3400},  {34, 2520},  {35, 1760}, {36, 1500}};

inline const std::map<Int_t, Double_t> GAS_PRESSURE_TORR = {
    {25, 0},   {26, 20},  {27, 40},  {28, 60},  {29, 80},  {30, 100},
    {31, 120}, {32, 140}, {33, 160}, {34, 180}, {35, 200}, {36, 220}};

const Double_t BEAM_ENERGY_MEV = 92.00;
const Double_t BEAM_ENERGY_ERR_MEV = 0.05;

// Mass number of the projectile (beam energy per nucleon = E / BEAM_A).
const Int_t BEAM_A = 37;

const Double_t FIT_HALF_WIDTH_ADC = 2000.0;

// --- Stopping-power / energy-loss geometry (the experimental Si-calib setup) -
// MUSIC active-gas length along the beam (cm).
const Double_t DETECTOR_LENGTH_CM = 35.0;
// Integration step through the gas (um).
const Double_t GAS_SEGMENT_MICRON = 100.0;
// PPAC: four single mylar foils in the beam.
const Double_t PPAC_LAYER_MICRON = 1.5;
const Int_t PPAC_N_LAYERS = 4;
// Ti entrance / exit window thicknesses (mg/cm^2).
const Double_t TI_FRONT_MG_CM2 = 0.9;
const Double_t TI_EXIT_MG_CM2 = 1.3;
// Run-25 survival threshold used by the stopping-power diagnostic (MeV).
const Double_t STOPPING_THRESHOLD_MEV = 71.35;

// Mulgin gain-group boundary: runs <= this are group A, the rest group B (the
// gain shift occurred between runs 22 and 23).
const Int_t GAIN_GROUP_BOUNDARY_RUN = 22;

// One Mulgin calibration anchor run: its TOF beam energy, which windows were in
// the beam, the measured Si ADC centroid/sigma, and its gain group (0=A, 1=B).
struct AnchorRun {
  Int_t run;
  Double_t tof_MeV;
  Bool_t front_Ti; // 0.9 mg/cm^2 Ti entrance window
  Bool_t exit_Ti;  // 1.3 mg/cm^2 Ti exit window
  Double_t adc_centroid;
  Double_t adc_sigma;
  Int_t group; // 0 = A, 1 = B
};

// The 5 anchor runs the Mulgin two-parameter fit is constrained on.
inline const std::vector<AnchorRun> ANCHOR_RUNS = {
    {21, 91.87, kFALSE, kFALSE, 15332.8, 48.8, 0},
    {22, 71.35, kFALSE, kFALSE, 10061.9, 42.4, 0},
    {23, 71.35, kTRUE, kFALSE, 9604.6, 120.4, 1},
    {24, 92.00, kTRUE, kFALSE, 13344.3, 119.4, 1},
    {25, 92.00, kTRUE, kTRUE, 10189.3, 140.8, 1}};

inline const TString FIT_RESULTS_FILE = "SiFits.root";
inline const TString DIAGNOSTIC_RESULTS_FILE = "Diagnostic.root";
inline const TString CALIBRATION_RESULTS_FILE = "SiCalibration.root";
inline const TString MULGIN_CALIBRATION_FILE = "MulginCalibration.root";

} // namespace SiCalib

#endif
