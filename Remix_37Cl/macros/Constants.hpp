#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

namespace Paths {
inline TString ProjectRootOf(const char *file) {
  TString path = file;
  path.ReplaceAll("/./", "/");
  TString macros_dir = gSystem->DirName(path);
  return gSystem->DirName(macros_dir);
}
} // namespace Paths

namespace Constants {

const Int_t N_STRIPS = 18;
const Int_t MAX_TRACE_SAVES = 10;

// Histogram axis ranges in MeV (sim is calibrated truth; no ADC).
const Double_t STRIP_E_MIN_MEV = -0.2;
const Double_t STRIP_E_MAX_MEV = 5.0;
const Double_t TOTAL_E_MIN_MEV = -0.5;
const Double_t TOTAL_E_MAX_MEV = 50.0;

} // namespace Constants

#endif
