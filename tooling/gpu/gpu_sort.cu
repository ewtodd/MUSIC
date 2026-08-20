#include "gpu_interface.h"

#include <cuda_runtime.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>

#include <cstdint>
#include <cstring>
#include <iostream>
#include <vector>

// Mirrors layout of RawHit in include/BinaryUtils.hpp: 24 bytes, timestamp at offset 8.
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

// Ship only 8-byte timestamp keys + 4-byte uint32 indices for a key/value radix sort;
// drops device footprint from ~48 to ~24 B/hit and avoids thrust comparator merge.
extern "C" int gpu_sort_hits_by_timestamp(void *hits, long long n_hits) {
  if (n_hits <= 0)
    return 0;

// uint32 indices suffice: a single subfile's hit count is far below 2^32.
  if (static_cast<unsigned long long>(n_hits) > 0xFFFFFFFFULL) {
    std::cerr << "[GPU sort] n_hits " << n_hits
              << " exceeds uint32 index range" << std::endl;
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
      std::cerr << "[GPU sort] H2D copy failed: " << cudaGetErrorString(err)
                << std::endl;
      return 1;
    }

    thrust::sort_by_key(d_keys.begin(), d_keys.end(), d_idx.begin());

    err = cudaGetLastError();
    if (err != cudaSuccess) {
      std::cerr << "[GPU sort] sort kernel failed: " << cudaGetErrorString(err)
                << std::endl;
      return 2;
    }

    thrust::copy(d_idx.begin(), d_idx.end(), h_idx.begin());

    err = cudaGetLastError();
    if (err != cudaSuccess) {
      std::cerr << "[GPU sort] D2H copy failed: " << cudaGetErrorString(err)
                << std::endl;
      return 3;
    }

    // Gather the hits into sorted order on the host.
    std::vector<RawHitGPU> sorted(n_hits);
    for (long long i = 0; i < n_hits; i++)
      sorted[i] = host_hits[h_idx[i]];
    std::memcpy(host_hits, sorted.data(),
                static_cast<size_t>(n_hits) * sizeof(RawHitGPU));
  } catch (const std::exception &e) {
    std::cerr << "[GPU sort] exception: " << e.what() << std::endl;
    return 4;
  } catch (...) {
    std::cerr << "[GPU sort] unknown exception" << std::endl;
    return 5;
  }

  return 0;
}
