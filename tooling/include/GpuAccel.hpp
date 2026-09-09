#ifndef GPU_ACCEL_HPP
#define GPU_ACCEL_HPP

#include "Constants.hpp"
#include <Rtypes.h>
#include <dlfcn.h>
#include <iostream>
#include <mutex>

/**
 * @brief Runtime access to the CUDA timestamp-sort kernel.
 *
 * The kernel lives in `libgpuaccel.so`, built from `tooling/gpu/` and loaded
 * with `dlopen` at first use rather than linked. That keeps the tooling
 * runnable on a machine with no CUDA runtime: if the library or its symbol is
 * missing, every entry point here degrades to reporting unavailable and the
 * caller takes the CPU path.
 *
 * All state is process-global and initialised once.
 */
class GpuAccel {
public:
  /// @brief Signature of the sort entry point: buffer, element count, status.
  typedef int (*SortFunc)(void *, long long);

  /**
   * @brief Load the GPU library and resolve the sort symbol. Idempotent.
   *
   * Does nothing and reports unavailable when
   * `Constants::cfg.USE_GPU_ACCELERATION` is off. The library path is baked in
   * at build time by the Makefile, falling back to a bare `libgpuaccel.so`
   * resolved through `LD_LIBRARY_PATH`.
   *
   * @return `kTRUE` if the kernel is usable. Repeat calls return the first
   *         result without retrying the load.
   *
   * @note Diagnostics go to stdout on success and stderr on failure; a failed
   *       load is reported, not raised.
   */
  static Bool_t Init();

  /// @brief Whether the kernel loaded successfully.
  /// @return `kFALSE` before Init() has run, and whenever the load failed.
  static Bool_t Available();

  /// @brief The resolved sort entry point.
  /// @return The function pointer, or null if unavailable. Check Available()
  ///         rather than testing this.
  static SortFunc GetSort();

  /**
   * @brief Claim one of the concurrent GPU sort slots.
   *
   * The pipeline processes subfiles in parallel, so without a cap they would
   * all dispatch to the GPU at once and exhaust device memory.
   * `Constants::cfg.MAX_GPU_CONCURRENT_SORTS` bounds how many may be in flight.
   *
   * @return `kTRUE` if a slot was taken, in which case the caller **must** pair
   *         it with ReleaseSortSlot(). `kFALSE` if all slots are busy, in which
   *         case the caller should take the CPU path rather than wait.
   *
   * @warning Non-blocking, and not RAII. An early return between acquire and
   *          release leaks a slot for the lifetime of the process.
   */
  static Bool_t TryAcquireSortSlot();

  /// @brief Give back a slot claimed by TryAcquireSortSlot().
  /// @warning Call exactly once per successful acquire. It is not guarded
  ///          against an unmatched call, which would let the cap drift upward.
  static void ReleaseSortSlot();
};

#endif
