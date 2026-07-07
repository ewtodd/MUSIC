#ifndef DEDUPSTRATEGY_HPP
#define DEDUPSTRATEGY_HPP

// Dedup strategy for multi-hit resolution within an event.
enum DedupStrategy {
  kCLOSEST_TO_FIRST_GRID,
  kCLOSEST_TO_LAST_GRID,
  kSMALLEST_ENERGY,
  kLARGEST_ENERGY,
  kEARLIEST_TIMESTAMP,
  kLATEST_TIMESTAMP
};

#endif
