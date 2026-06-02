#ifndef REMIX_CONSTANTS_HPP
#define REMIX_CONSTANTS_HPP

#include <Rtypes.h>
#include <TString.h>

// Sim-side ("Remix") constants and strip helpers, in namespace Remix to avoid
// clashing with the dataset config namespace Constants. These describe the
// simulation output (calibrated MeV truth; no ADC), distinct from the
// experimental analysis/<iso>/config/Constants.hpp.
namespace Remix {

const Int_t N_STRIPS = 18;
const Int_t MAX_TRACE_SAVES = 10;

const Bool_t IGNORE_SHORT_STRIPS = kFALSE;
const Bool_t IGNORE_STRIP_0 = kFALSE;
const Bool_t IGNORE_STRIP_17 = kFALSE;

// Histogram axis ranges in MeV (sim is calibrated truth; no ADC).
const Double_t STRIP_E_MIN_MEV = -0.2;
const Double_t STRIP_E_MAX_MEV = 5.0;
const Double_t TOTAL_E_MIN_MEV = -0.5;
const Double_t TOTAL_E_MAX_MEV = 50.0;

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

} // namespace Remix

#endif
