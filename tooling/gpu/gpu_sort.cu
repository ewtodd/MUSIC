#include "gpu_interface.h"

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

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

// Sort hits by timestamp. Rather than shipping the full 24-byte structs to the
// device and sorting them with a comparator (which forces thrust into a
// comparator merge sort needing ~2x the data in scratch), we ship only the
// 8-byte uint64 timestamp keys plus a 4-byte uint32 index and run a key/value
// sort -- a primitive-key radix sort. The device footprint per call drops from
// ~48 B/hit to ~24 B/hit (so more sorts fit concurrently) and H2D/D2H traffic
// drops from 24 down to 8 (down) + 4 (back) B/hit. The resulting permutation is
// applied to the host array by a gather; the full structs never touch the GPU.
extern "C" int gpu_sort_hits_by_timestamp(void *hits, long long n_hits) {
  if (n_hits <= 0)
    return 0;

  // Indices are uint32 to keep the device payload small; a single subfile's hit
  // count is far below 2^32, but guard anyway and let the caller fall back to
  // the CPU sort if that ever stops holding.
  if (static_cast<unsigned long long>(n_hits) > 0xFFFFFFFFULL) {
    fprintf(stderr, "[GPU sort] n_hits %lld exceeds uint32 index range\n",
            n_hits);
    return 6;
  }

  RawHitGPU *host_hits = static_cast<RawHitGPU *>(hits);

  try {
    thrust::host_vector<uint64_t> h_keys(n_hits);
    thrust::host_vector<uint32_t> h_idx(n_hits);
    for (long long i = 0; i < n_hits; i++) {
      h_keys[i] = host_hits[i].timestamp;
      h_idx[i] = static_cast<uint32_t>(i);
    }

    thrust::device_vector<uint64_t> d_keys = h_keys;
    thrust::device_vector<uint32_t> d_idx = h_idx;

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
      fprintf(stderr, "[GPU sort] H2D copy failed: %s\n",
              cudaGetErrorString(err));
      return 1;
    }

    thrust::sort_by_key(d_keys.begin(), d_keys.end(), d_idx.begin());

    err = cudaGetLastError();
    if (err != cudaSuccess) {
      fprintf(stderr, "[GPU sort] sort kernel failed: %s\n",
              cudaGetErrorString(err));
      return 2;
    }

    thrust::copy(d_idx.begin(), d_idx.end(), h_idx.begin());

    err = cudaGetLastError();
    if (err != cudaSuccess) {
      fprintf(stderr, "[GPU sort] D2H copy failed: %s\n",
              cudaGetErrorString(err));
      return 3;
    }

    // Gather the hits into sorted order on the host.
    std::vector<RawHitGPU> sorted(n_hits);
    for (long long i = 0; i < n_hits; i++)
      sorted[i] = host_hits[h_idx[i]];
    std::memcpy(host_hits, sorted.data(),
                static_cast<size_t>(n_hits) * sizeof(RawHitGPU));
  } catch (const std::exception &e) {
    fprintf(stderr, "[GPU sort] exception: %s\n", e.what());
    return 4;
  } catch (...) {
    fprintf(stderr, "[GPU sort] unknown exception\n");
    return 5;
  }

  return 0;
}
