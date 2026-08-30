#ifndef CHANNEL_TIMING_HPP
#define CHANNEL_TIMING_HPP

#include "Constants.hpp"
#include <Rtypes.h>
#include <map>
#include <utility>

// Per-channel TTF (trapezoidal trigger filter) timing offset lookup.
// The offsets live in Constants::cfg.ttfOffsetPs (empty = no correction).
// Returns 0 for channels not listed.
namespace Constants {

inline Long64_t LookupTTFOffsetPs(Int_t board, Int_t channel) {
  std::map<std::pair<Int_t, Int_t>, Long64_t>::const_iterator it =
      cfg.ttfOffsetPs.find(std::pair<Int_t, Int_t>(board, channel));
  return (it == cfg.ttfOffsetPs.end()) ? 0 : it->second;
}

} // namespace Constants

#endif
