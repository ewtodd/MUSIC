#ifndef GPU_ACCEL_HPP
#define GPU_ACCEL_HPP

#include "Constants.hpp"
#include <Rtypes.h>
#include <dlfcn.h>
#include <iostream>
#include <mutex>

class GpuAccel {
public:
  typedef int (*SortFunc)(void *, long long);

  static Bool_t Init();
  static Bool_t Available();
  static SortFunc GetSort();

  static Bool_t TryAcquireSortSlot();
  static void ReleaseSortSlot();
};

#endif
