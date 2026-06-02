#ifndef CHANNEL_TIMING_HPP
#define CHANNEL_TIMING_HPP

#include <Rtypes.h>

// Per-channel TTF (trapezoidal trigger filter) timing offset lookup. The
// offsets themselves live in the per-dataset `ttfOffsetPs` map in Constants.hpp
// (empty = no correction, e.g. 87Rb); this tooling function reads that map.
// Timing subtracts the returned ps offset from each hit's raw timestamp to put
// all channels on a common timing baseline.
namespace Constants {

Long64_t LookupTTFOffsetPs(Int_t board, Int_t channel);

} // namespace Constants

#endif
