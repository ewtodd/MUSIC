#include "GpuAccel.hpp"
#include "Constants.hpp"
#include <dlfcn.h>
#include <iostream>
#include <mutex>

namespace GpuAccel {

static SortFunc g_sort = nullptr;
static Bool_t g_available = kFALSE;
static Bool_t g_initialized = kFALSE;

Bool_t Init() {
  if (g_initialized)
    return g_available;
  g_initialized = kTRUE;

  if (!Constants::USE_GPU_ACCELERATION)
    return kFALSE;

  void *lib = dlopen("libgpuaccel.so", RTLD_LAZY);
  if (!lib) {
    std::cerr << "[GPU] Library not found, using CPU fallback: " << dlerror()
              << std::endl;
    return kFALSE;
  }

  g_sort = (SortFunc)dlsym(lib, "gpu_sort_hits_by_timestamp");
  g_available = (g_sort != nullptr) ? kTRUE : kFALSE;

  if (g_available)
    std::cout << "[GPU] Acceleration loaded successfully." << std::endl;
  else
    std::cerr << "[GPU] Symbols missing, using CPU fallback." << std::endl;

  return g_available;
}

Bool_t Available() { return g_available; }
SortFunc GetSort() { return g_sort; }

static std::mutex g_sort_mutex;
static Int_t g_sort_active = 0;

Bool_t TryAcquireSortSlot() {
  std::lock_guard<std::mutex> lk(g_sort_mutex);
  if (g_sort_active >= Constants::MAX_GPU_CONCURRENT_SORTS)
    return kFALSE;
  g_sort_active++;
  return kTRUE;
}

void ReleaseSortSlot() {
  std::lock_guard<std::mutex> lk(g_sort_mutex);
  g_sort_active--;
}

} // namespace GpuAccel
