#ifndef DEDUPSTRATEGY_HPP
#define DEDUPSTRATEGY_HPP

// Dedup strategy for multi-hit resolution within an event.
enum DedupStrategy {
  kCLOSEST_TO_FIRST_GRID,
  kSMALLEST_ENERGY,
  kLARGEST_ENERGY,
  kLATEST_TIMESTAMP,
  // Not a selection: any event with a repeated hit on an anode channel is
  // thrown away whole. Mirrors the upstream builder, which sets event_rej on a
  // second hit to de_l/de_r and never fills that event.
  kDISCARD
};

#endif
