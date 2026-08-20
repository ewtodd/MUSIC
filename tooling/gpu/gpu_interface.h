#ifndef GPU_INTERFACE_H
#define GPU_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

// Sort RawHit-layout structs by timestamp (uint64 at byte offset 8).
int gpu_sort_hits_by_timestamp(void *hits, long long n_hits);

#ifdef __cplusplus
}
#endif

#endif
