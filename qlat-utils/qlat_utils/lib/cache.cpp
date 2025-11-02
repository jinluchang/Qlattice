#include <qlat-utils/cache.h>

namespace qlat
{  //

void clear_mem_cache()
{
  TIMER_VERBOSE_FLOPS("clear_mem_cache");
  Long total_bytes = 0;
  for (Int i = 0; i < static_cast<Int>(MemType::SIZE); ++i) {
    MemType mem_type = static_cast<MemType>(i);
    const Long bytes = get_mem_cache(mem_type).mem_cache_size;
    total_bytes += bytes;
    get_mem_cache(mem_type).gc();
    displayln_info(
        0, fname + ssprintf(": mem_type=%s %ld bytes (%.3f GB) freed.",
                            show(mem_type).c_str(), bytes,
                            (RealD)bytes / (1024.0 * 1024.0 * 1024.0)));
  }
  displayln_info(
      0, fname + ssprintf(": total %ld bytes (%.3f GB) freed.", total_bytes,
                          (RealD)total_bytes / (1024.0 * 1024.0 * 1024.0)));
  timer.flops += total_bytes;
}

void clear_all_caches()
// (1) Clear all cache using CacheBase (or Cache which internally uses
// CacheBase).
//
// (2) Clear memory cache (for vector memory freed but not cached).
{
  TIMER_VERBOSE("clear_all_caches");
  for (auto it = get_all_caches().cbegin(); it != get_all_caches().cend();
       ++it) {
    (*it)->clear();
  }
  clear_mem_cache();
}

}  // namespace qlat
