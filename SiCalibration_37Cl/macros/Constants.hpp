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

// Gas pressure for each pressure-scan run (Torr).
inline const std::map<Int_t, Double_t> GAS_PRESSURE_TORR = {
    {25, 0},   {26, 20},  {27, 40},  {28, 60},  {29, 80},  {30, 100},
    {31, 120}, {32, 140}, {33, 160}, {34, 180}, {35, 200}, {36, 220}};

const Double_t E_RUN_21_MEV = 91.87;
const Double_t E_RUN_22_MEV = 71.35;

const Double_t BEAM_ENERGY_MEV = 92.00;
const Double_t BEAM_ENERGY_ERR_MEV = 0.05;

const Double_t FIT_HALF_WIDTH_ADC = 2000.0;

inline const TString FIT_RESULTS_FILE = "SiCalibration_FitResults.root";
inline const TString ANALYSIS_RESULTS_FILE = "SiCalibration_Results.root";

} // namespace Constants

#endif
