#ifndef GPU_ACCEL_HPP
#define GPU_ACCEL_HPP

#include <TROOT.h>

namespace GpuAccel {

typedef int (*SortFunc)(void *, long long);

Bool_t Init();
Bool_t Available();
SortFunc GetSort();

Bool_t TryAcquireSortSlot();
void ReleaseSortSlot();

} // namespace GpuAccel

#endif
