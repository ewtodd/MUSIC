#ifndef SLOT_LAYOUT_HPP
#define SLOT_LAYOUT_HPP

#include <Rtypes.h>

// Structural slot layout for the EventState / PerChannelData 36-element arrays.
// These indices are fixed by the tooling's data model (see EventBuilder), not
// by the dataset, so they live with the tooling rather than per-dataset config.
namespace Constants {

const Int_t ARR_SLOT_STRIP_0 = 0;
const Int_t ARR_SLOT_STRIP_17 = 33;
const Int_t ARR_SLOT_CATHODE = 34;
const Int_t ARR_SLOT_GRID = 35;

// Number of slots: indices run 0..ARR_SLOT_GRID, so the EventState /
// PerChannelData arrays are this long.
const Int_t N_ARR_SLOTS = ARR_SLOT_GRID + 1;

} // namespace Constants

#endif
