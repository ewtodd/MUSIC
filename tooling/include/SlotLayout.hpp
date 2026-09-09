#ifndef SLOT_LAYOUT_HPP
#define SLOT_LAYOUT_HPP

#include <Rtypes.h>

/**
 * @file SlotLayout.hpp
 * @brief Fixed array indices for the per-event channel arrays.
 *
 * These positions are set by the tooling's data model — the 36-element
 * `EventState` and `PerChannelData` arrays that EventBuilder fills — and not by
 * any dataset. That is why they live here rather than in the per-dataset
 * configuration: changing one changes the meaning of every event array in the
 * tooling, for every experiment.
 */
namespace Constants {

/// First anode strip, occupying the low end of the array.
const Int_t ARR_SLOT_STRIP_0 = 0;
/// Last anode strip. Strips run contiguously from #ARR_SLOT_STRIP_0 to here.
const Int_t ARR_SLOT_STRIP_17 = 33;
const Int_t ARR_SLOT_CATHODE = 34; ///< Cathode signal.
const Int_t ARR_SLOT_GRID = 35;    ///< Frisch grid signal.

/// Length of the per-event channel arrays. Indices run 0 to #ARR_SLOT_GRID
/// inclusive, so this is one past the last valid slot.
const Int_t N_ARR_SLOTS = ARR_SLOT_GRID + 1;

} // namespace Constants

#endif
