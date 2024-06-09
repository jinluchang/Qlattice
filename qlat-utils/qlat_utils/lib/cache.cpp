#include <qlat-utils/cache.h>

namespace qlat
{  //

void clear_mem_cache()
{
  TIMER_VERBOSE_FLOPS("clear_mem_cache");
  Long total_bytes = 0;
  total_bytes += get_mem_cache(false).mem_cache_size;
  total_bytes += get_mem_cache(true).mem_cache_size;
  get_mem_cache(false).gc();
  get_mem_cache(true).gc();
  displayln_info(
      0, fname + ssprintf(": %ld bytes (%.3f GB) freed.", total_bytes,
                          (double)total_bytes / (1024.0 * 1024.0 * 1024.0)));
  timer.flops += total_bytes;
}

void clear_all_caches()
// (1) Clear all cache using CacheBase (or Cache which internally uses
// CacheBase).
//
// (2) Clear memory cache (for vector and vector_acc memory freed but not
// cached).
{
  TIMER_VERBOSE("clear_all_caches");
  for (auto it = get_all_caches().cbegin(); it != get_all_caches().cend();
       ++it) {
    (*it)->clear();
  }
  clear_mem_cache();
}

}  // namespace qlat
