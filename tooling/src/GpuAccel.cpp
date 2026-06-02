#include "GpuAccel.hpp"

GpuAccel::SortFunc g_sort = nullptr;
Bool_t g_available = kFALSE;
Bool_t g_initialized = kFALSE;
std::mutex g_sort_mutex;
Int_t g_sort_active = 0;

Bool_t GpuAccel::Init() {
  if (g_initialized)
    return g_available;
  g_initialized = kTRUE;

  if (!Constants::USE_GPU_ACCELERATION)
    return kFALSE;

  // Absolute path to the tooling GPU lib, injected at build by the Makefile
  // (it lives with the tooling, not the dataset). Falls back to a bare name
  // resolved via LD_LIBRARY_PATH if the macro is somehow undefined.
#ifndef MUSIC_GPU_LIB
#define MUSIC_GPU_LIB "libgpuaccel.so"
#endif
  void *lib = dlopen(MUSIC_GPU_LIB, RTLD_LAZY);
  if (!lib) {
    std::cerr << "[GPU] Library not found, using CPU fallback: " << dlerror()
              << std::endl;
    return kFALSE;
  }

  g_sort = (GpuAccel::SortFunc)dlsym(lib, "gpu_sort_hits_by_timestamp");
  g_available = (g_sort != nullptr) ? kTRUE : kFALSE;

  if (g_available)
    std::cout << "[GPU] Acceleration loaded successfully." << std::endl;
  else
    std::cerr << "[GPU] Symbols missing, using CPU fallback." << std::endl;

  return g_available;
}

Bool_t GpuAccel::Available() { return g_available; }

GpuAccel::SortFunc GpuAccel::GetSort() { return g_sort; }

Bool_t GpuAccel::TryAcquireSortSlot() {
  std::lock_guard<std::mutex> lk(g_sort_mutex);
  if (g_sort_active >= Constants::MAX_GPU_CONCURRENT_SORTS)
    return kFALSE;
  g_sort_active++;
  return kTRUE;
}

void GpuAccel::ReleaseSortSlot() {
  std::lock_guard<std::mutex> lk(g_sort_mutex);
  g_sort_active--;
}
