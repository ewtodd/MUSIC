#include "gpu_interface.h"

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>

#include <cstdint>
#include <cstdio>

// Mirrors layout of RawHit in include/BinaryUtils.hpp.
// The host struct uses ROOT types (UShort_t, ULong64_t, UInt_t) but the
// memory layout under the standard ABI is:
//   offset 0:  board    (2 bytes)
//   offset 2:  channel  (2 bytes)
//   offset 4:  energy   (2 bytes)
//   offset 6:  padding  (2 bytes for ULong64_t alignment)
//   offset 8:  timestamp (8 bytes)
//   offset 16: flags    (4 bytes)
//   offset 20: padding  (4 bytes for struct alignment)
// Total: 24 bytes.
struct RawHitGPU {
  uint16_t board;
  uint16_t channel;
  uint16_t energy;
  uint16_t _pad0;
  uint64_t timestamp;
  uint32_t flags;
  uint32_t _pad1;
};

static_assert(sizeof(RawHitGPU) == 24,
              "RawHitGPU must be 24 bytes to match host RawHit layout");
static_assert(offsetof(RawHitGPU, timestamp) == 8,
              "timestamp must sit at byte offset 8");

struct TimestampLess {
  __host__ __device__ bool operator()(const RawHitGPU &a,
                                      const RawHitGPU &b) const {
    return a.timestamp < b.timestamp;
  }
};

extern "C" int gpu_sort_hits_by_timestamp(void *hits, long long n_hits) {
  if (n_hits <= 0)
    return 0;

  RawHitGPU *host_hits = static_cast<RawHitGPU *>(hits);

  try {
    thrust::device_vector<RawHitGPU> d_hits(host_hits, host_hits + n_hits);

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
      fprintf(stderr, "[GPU sort] H2D copy failed: %s\n",
              cudaGetErrorString(err));
      return 1;
    }

    thrust::sort(d_hits.begin(), d_hits.end(), TimestampLess());

    err = cudaGetLastError();
    if (err != cudaSuccess) {
      fprintf(stderr, "[GPU sort] sort kernel failed: %s\n",
              cudaGetErrorString(err));
      return 2;
    }

    thrust::copy(d_hits.begin(), d_hits.end(), host_hits);

    err = cudaGetLastError();
    if (err != cudaSuccess) {
      fprintf(stderr, "[GPU sort] D2H copy failed: %s\n",
              cudaGetErrorString(err));
      return 3;
    }
  } catch (const std::exception &e) {
    fprintf(stderr, "[GPU sort] exception: %s\n", e.what());
    return 4;
  } catch (...) {
    fprintf(stderr, "[GPU sort] unknown exception\n");
    return 5;
  }

  return 0;
}
