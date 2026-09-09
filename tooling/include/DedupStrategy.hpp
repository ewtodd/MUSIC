#ifndef DEDUPSTRATEGY_HPP
#define DEDUPSTRATEGY_HPP

/**
 * @brief How to resolve more than one hit on the same channel within an event.
 *
 * The first four values pick a survivor and keep the event. #kDISCARD is not a
 * selection at all — it throws the whole event away.
 */
enum DedupStrategy {
  /// Keep the hit whose timestamp is nearest the first grid hit.
  kCLOSEST_TO_FIRST_GRID,
  kSMALLEST_ENERGY,  ///< Keep the lowest-energy hit.
  kLARGEST_ENERGY,   ///< Keep the highest-energy hit.
  kLATEST_TIMESTAMP, ///< Keep the last hit to arrive.
  /// Discard the entire event if any anode channel has a repeated hit. Mirrors
  /// the upstream builder, which sets `event_rej` on a second hit to `de_l` or
  /// `de_r` and never fills that event.
  kDISCARD
};

#endif
