#ifndef REMIX_CONSTANTS_HPP
#define REMIX_CONSTANTS_HPP

#include <Rtypes.h>
#include <TString.h>

// Sim-side ("Remix") constants and strip helpers, in namespace Remix to avoid
// clashing with the dataset config namespace Constants. These describe the
// simulation output (calibrated MeV truth; no ADC), distinct from the
// experimental analysis/<iso>/config/Constants.hpp.
namespace Remix {

const Int_t MAX_TRACE_SAVES = 10;

const Bool_t IGNORE_SHORT_STRIPS = kFALSE;
const Bool_t IGNORE_STRIP_0 = kFALSE;
const Bool_t IGNORE_STRIP_17 = kFALSE;

// Histogram axis ranges in MeV (sim is calibrated truth; no ADC).
const Double_t STRIP_E_MIN_MEV = -0.2;
const Double_t STRIP_E_MAX_MEV = 5.0;
const Double_t TOTAL_E_MIN_MEV = -0.5;
const Double_t TOTAL_E_MAX_MEV = 50.0;

inline Int_t FirstStrip() { return IGNORE_STRIP_0 ? 1 : 0; }
inline Int_t LastStrip() { return IGNORE_STRIP_17 ? 16 : 17; }

} // namespace Remix

#endif
