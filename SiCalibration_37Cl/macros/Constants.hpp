#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <map>
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

const Bool_t INCLUDE_PPAC = kTRUE;

// Energy (MeV) that the 37Cl beam has when it reaches the Si detector for
// each zero-pressure anchor run.  PPAC_OUT = PPAC not in beam; PPAC_IN =
// PPAC (4×1.5 μm mylar) in beam.  All PPAC_IN values are model averages
// over models 1–4 (Ziegler, ATIMA12LS, ATIMA12NoLS, ATIMA14Weick) from
// StoppingPowerDiagnostic.
//
// Run 21:  TOF = 91.87 MeV, no Ti windows.
const Double_t E_RUN_21_PPAC_OUT = 91.87;
const Double_t E_RUN_21_PPAC_IN = E_RUN_21_PPAC_OUT - 18.1162;
// Run 22:  TOF = 71.35 MeV, no Ti windows.
const Double_t E_RUN_22_PPAC_OUT = 71.35;
const Double_t E_RUN_22_PPAC_IN = E_RUN_22_PPAC_OUT - 19.5146;
// Run 23:  TOF = 71.35 MeV, front Ti window only (0.9 mg/cm²).
const Double_t E_RUN_23_PPAC_OUT = 58.6312;
const Double_t E_RUN_23_PPAC_IN = 38.3174;
// Run 24:  TOF = 92.00 MeV, front Ti window only (0.9 mg/cm²).
const Double_t E_RUN_24_PPAC_OUT = 80.0511;
const Double_t E_RUN_24_PPAC_IN = 61.2780;
// Run 25:  TOF = 92.00 MeV, both Ti windows (0.9 + 1.3 mg/cm²),
const Double_t E_RUN_25_PPAC_OUT = 62.5413;
const Double_t E_RUN_25_PPAC_IN = 42.3088;

inline const std::vector<Int_t> RUN_NUMBERS = {21, 22, 23, 24, 25, 26, 27, 28,
                                               29, 30, 31, 32, 33, 34, 35, 36};

inline const TString DATA_BASE_PATH = "/home/e-work/LabData/MUSIC/37Cl";

inline TString GetRunFilepath(Int_t run) {
  return Form("%s/run_%d/ROOT/DataR_run_%d.root", DATA_BASE_PATH.Data(), run,
              run);
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

const Double_t FIT_HALF_WIDTH_ADC = 2000.0;

inline const TString FIT_RESULTS_FILE = "SiFits.root";
inline const TString DIAGNOSTIC_RESULTS_FILE = "Diagnostic.root";
inline const TString CALIBRATION_RESULTS_FILE = "SiCalibration.root";
inline const TString MULGIN_CALIBRATION_FILE = "MulginCalibration.root";

} // namespace Constants

#endif
