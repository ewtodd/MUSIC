#ifndef REMIX_CONSTANTS_HPP
#define REMIX_CONSTANTS_HPP

#include <Rtypes.h>
#include <TString.h>

namespace Remix {

const Int_t MAX_TRACE_SAVES = 10;

const Bool_t IGNORE_SHORT_STRIPS = kFALSE;
const Bool_t IGNORE_STRIP_0 = kFALSE;
const Bool_t IGNORE_STRIP_17 = kFALSE;

// Histogram axis ranges (a.u., normalized truth; no ADC).
const Double_t STRIP_E_MIN_NORMED = -0.2;
const Double_t STRIP_E_MAX_NORMED = 5.0;
const Double_t TOTAL_E_MIN_NORMED = -0.5;
const Double_t TOTAL_E_MAX_NORMED = 50.0;

inline Int_t FirstStrip() { return IGNORE_STRIP_0 ? 1 : 0; }
inline Int_t LastStrip() { return IGNORE_STRIP_17 ? 16 : 17; }

} // namespace Remix

#endif
