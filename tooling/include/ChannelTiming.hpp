#ifndef CHANNEL_TIMING_HPP
#define CHANNEL_TIMING_HPP

#include "Constants.hpp"
#include <Rtypes.h>
#include <map>
#include <utility>

/**
 * @file ChannelTiming.hpp
 * @brief Per-channel trapezoidal-trigger-filter timing offsets.
 */
namespace Constants {

/**
 * @brief Timing offset to apply to one channel's TTF timestamp.
 *
 * The trapezoidal trigger filter introduces a channel-dependent delay, so
 * timestamps from different boards and channels need aligning before events can
 * be built across them. The offsets are a per-dataset measurement, held in
 * `Constants::cfg.ttfOffsetPs`.
 *
 * @param board   Digitiser board id.
 * @param channel Channel index on that board.
 *
 * @return The offset in picoseconds, or `0` for a channel with no entry. An
 *         empty `ttfOffsetPs` therefore disables the correction entirely rather
 *         than failing.
 */
inline Long64_t LookupTTFOffsetPs(Int_t board, Int_t channel) {
  std::map<std::pair<Int_t, Int_t>, Long64_t>::const_iterator it =
      cfg.ttfOffsetPs.find(std::pair<Int_t, Int_t>(board, channel));
  return (it == cfg.ttfOffsetPs.end()) ? 0 : it->second;
}

} // namespace Constants

#endif
