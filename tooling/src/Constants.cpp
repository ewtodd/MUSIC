#include "Constants.hpp"

namespace Constants {

// Per-channel TTF timing offset (ps) from the per-dataset ttfOffsetPs map.
// Returns 0 for channels not listed (e.g. every channel when the map is empty,
// as for 87Rb). For 37Cl this corrects board 2's misconfigured TTF delay.
Long64_t LookupTTFOffsetPs(Int_t board, Int_t channel) {
  std::map<std::pair<Int_t, Int_t>, Long64_t>::const_iterator it =
      ttfOffsetPs.find(std::pair<Int_t, Int_t>(board, channel));
  return (it == ttfOffsetPs.end()) ? 0 : it->second;
}

} // namespace Constants
