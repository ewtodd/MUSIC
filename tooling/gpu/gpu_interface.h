#ifndef GPU_INTERFACE_H
#define GPU_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

/* Sort an array of RawHit-layout structs by timestamp (uint64 at byte
   offset 8). hits: pointer to contiguous array of structs. n_hits: number
   of elements. Returns 0 on success, nonzero on failure. */
int gpu_sort_hits_by_timestamp(void *hits, long long n_hits);

#ifdef __cplusplus
}
#endif

#endif
