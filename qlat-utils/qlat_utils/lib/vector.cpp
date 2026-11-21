#include <qlat-utils/vector.h>

#include <cstdint>

namespace qlat
{  //

static void* alloc_mem_aligned(const Long size, const MemType mem_type,
                               const Long min_alignment)
// Allocate `size` bytes aligned to at least `min_alignment` bytes (usually
// alignof(M) of the element type).  The configured alignment
// (get_alignment()) is raised to `min_alignment` when it is smaller, so an
// env override (q_alignment / q_alignment_comm) can never under-align an
// over-aligned element type.  std::aligned_alloc() also requires the size to
// be a multiple of the alignment, so the size is rounded up.
{
  void* ptr = NULL;
#if defined QLAT_NO_ALIGNED_ALLOC
  (void)min_alignment;
  ptr = malloc(size);
#else
  Long alignment = get_effective_alignment(mem_type, min_alignment);
  while (true) {
    const Long alloc_size = ((size + alignment - 1) / alignment) * alignment;
    ptr = std::aligned_alloc(alignment, alloc_size);
    if (ptr != NULL or alignment <= min_alignment) {
      break;
    }
    alignment /= 2;
    if (alignment < min_alignment) {
      alignment = min_alignment;
    }
  }
  if (ptr != NULL) {
    qassert(((uintptr_t)ptr) % min_alignment == 0);
  }
#endif
  return ptr;
}

static void* alloc_mem_alloc_impl(const Long size, const MemType mem_type,
                                  const Long min_alignment)
{
  TIMER_FLOPS("alloc_mem_alloc");
  timer.flops += size;
  if (size == 0) {
    return NULL;
  }
  static MemoryStats& ms = get_mem_stats();
  ms.alloc[static_cast<Int>(mem_type)] += size;
  void* ptr = NULL;
#ifdef QLAT_USE_ACC
  const MemType eff_mem_type = get_eff_mem_type(mem_type);
  if (eff_mem_type == MemType::Cpu) {
    TIMER_FLOPS("alloc_mem_alloc(MemType::Cpu)");
    timer.flops += size;
    ptr = alloc_mem_aligned(size, mem_type, min_alignment);
    qassert(ptr != NULL);
    set_mem(ptr, 0, size, mem_type);
  } else if (eff_mem_type == MemType::Acc) {
    TIMER_FLOPS("alloc_mem_alloc(MemType::Acc)");
    timer.flops += size;
    qacc_Error err = qacc_GetLastError();
    if (qacc_Success != err) {
      qerr(fname + ssprintf(": ACC error '%s' before qacc_MallocManaged.",
                            qacc_GetErrorString(err)));
    }
    err = qacc_Malloc(&ptr, size);
    if (qacc_Success != err) {
      qerr(fname + ssprintf(": ACC error '%s', size=%ld, ptr=%lX.",
                            qacc_GetErrorString(err), size, ptr));
    }
    qassert(ptr != NULL);
    set_mem(ptr, 0, size, mem_type);
  } else if (eff_mem_type == MemType::Uvm) {
    TIMER_FLOPS("alloc_mem_alloc(MemType::Uvm)");
    timer.flops += size;
    qacc_Error err = qacc_GetLastError();
    if (qacc_Success != err) {
      qerr(fname + ssprintf(": ACC error '%s' before qacc_MallocManaged.",
                            qacc_GetErrorString(err)));
    }
    err = qacc_MallocManaged(&ptr, size);
    if (qacc_Success != err) {
      qerr(fname + ssprintf(": ACC error '%s', size=%ld, ptr=%lX.",
                            qacc_GetErrorString(err), size, ptr));
    }
    qassert(ptr != NULL);
    set_mem(ptr, 0, size, mem_type);
  } else {
    qassert(false);
  }
#else
  ptr = alloc_mem_aligned(size, mem_type, min_alignment);
  qassert_info(ptr != NULL,
               { printf("alloc_mem_alloc: %lX %ld\n", (Long)ptr, size); });
  set_mem(ptr, 0, size, mem_type);
#endif
  if (ptr != NULL) {
    qassert(((uintptr_t)ptr) % min_alignment == 0);
  }
  return ptr;
}

void* alloc_mem_alloc(const Long size, const MemType mem_type)
{
  return alloc_mem_alloc_impl(size, mem_type, 1);
}

void free_mem_free(void* ptr, const Long size, const MemType mem_type)
{
  TIMER_FLOPS("free_mem_free");
  timer.flops += size;
  if (size == 0) {
    return;
  }
  static MemoryStats& ms = get_mem_stats();
  ms.alloc[static_cast<Int>(mem_type)] -= size;
#ifdef QLAT_USE_ACC
  const MemType eff_mem_type = get_eff_mem_type(mem_type);
  if (eff_mem_type == MemType::Cpu) {
    TIMER_FLOPS("free_mem_free(MemType::Cpu)");
    timer.flops += size;
    free(ptr);
  } else if (eff_mem_type == MemType::Acc) {
    TIMER_FLOPS("free_mem_free(MemType::Acc)");
    timer.flops += size;
    qacc_Error err = qacc_Free(ptr);
    if (qacc_Success != err) {
      if (qacc_ErrorUnloading != err) {
        qerr(fname + ssprintf(": ACC error '%s' (%d) after qacc_Free.",
                              qacc_GetErrorString(err), err));
      }
    }
  } else if (eff_mem_type == MemType::Uvm) {
    TIMER_FLOPS("free_mem_free(MemType::Uvm)");
    timer.flops += size;
    qacc_Error err = qacc_Free(ptr);
    if (qacc_Success != err) {
      if (qacc_ErrorUnloading != err) {
        qerr(fname + ssprintf(": ACC error '%s' (%d) after qacc_Free.",
                              qacc_GetErrorString(err), err));
      }
    }
  } else {
    qassert(false);
  }
#else
  free(ptr);
#endif
}

void MemCache::init(const MemType mem_type_, const Long mem_cache_max_size_)
{
  mem_type = mem_type_;
  mem_cache_max_size = mem_cache_max_size_;
  if (mem_cache_max_size < 0) {
    mem_cache_max_size = get_mem_cache_max_size(mem_type);
  }
  mem_cache_size = 0;
}

void MemCache::add(void* ptr, const Long size, const Long alignment)
{
  qassert(size > 0);
  mem_cache_size += size;
  static MemoryStats& ms = get_mem_stats();
  ms.cache[static_cast<Int>(mem_type)] += size;
  db.insert(std::pair<Long, void*>(size, ptr));
  db_alignment[ptr] = alignment;
  if (mem_cache_size > mem_cache_max_size) {
    gc();
  }
}

void* MemCache::del(const Long size, const Long min_alignment)
{
  const auto range = db.equal_range(size);
  for (auto iter = range.first; iter != range.second; ++iter) {
    void* ptr = iter->second;
    const auto alignment_iter = db_alignment.find(ptr);
    const Long alignment =
        (alignment_iter == db_alignment.end()) ? 1 : alignment_iter->second;
    if (alignment < min_alignment) {
      continue;
    }
    mem_cache_size -= size;
    static MemoryStats& ms = get_mem_stats();
    ms.cache[static_cast<Int>(mem_type)] -= size;
    db.erase(iter);
    db_alignment.erase(ptr);
    return ptr;
  }
  return NULL;
}

void MemCache::gc()
{
  if (mem_cache_size == 0) {
    return;
  }
  TIMER_FLOPS("MemCache::gc()");
  timer.flops += mem_cache_size;
  static MemoryStats& ms = get_mem_stats();
  ms.cache[static_cast<Int>(mem_type)] -= mem_cache_size;
  for (auto iter = db.cbegin(); iter != db.cend(); ++iter) {
    const Long size = iter->first;
    void* ptr = iter->second;
    qassert(ptr != NULL);
    free_mem_free(ptr, size, mem_type);
    mem_cache_size -= size;
  }
  qassert(mem_cache_size == 0);
  db.clear();
  db_alignment.clear();
}

std::vector<MemCache> mk_mem_cache_vec()
{
  std::vector<MemCache> ret;
  const Int size = static_cast<Int>(MemType::SIZE);
  ret.resize(size);
  for (Int i = 0; i < size; ++i) {
    MemType mem_type = static_cast<MemType>(i);
    ret[i].init(mem_type);
  }
  return ret;
}

void* alloc_mem(const Long min_size, const MemType mem_type)
{
  return alloc_mem(min_size, mem_type, 1);
}

void* alloc_mem(const Long min_size, const MemType mem_type,
                const Long min_alignment)
{
  if (min_size <= 0) {
    return NULL;
  }
  TIMER_FLOPS("alloc_mem");
  timer.flops += min_size;
  const Long chunk_size = get_mem_chunk_size(mem_type);
  const Long size = get_chunked_mem_size(chunk_size, min_size);
  MemCache& cache = get_mem_cache(mem_type);
  void* ptr = cache.del(size, min_alignment);
  if (NULL != ptr) {
    return ptr;
  }
  if (size + cache.mem_cache_size > cache.mem_cache_max_size) {
    cache.gc();
  }
  static MemoryStats& ms = get_mem_stats();
  if (size + ms.total() > get_alloc_mem_max_size()) {
    displayln_info(
        fname + ssprintf(": alloc %.3lf (GB) memory (current total %.3lf (GB))",
                         (RealD)min_size / (1024.0 * 1024.0 * 1024.0),
                         (RealD)ms.total() / (1024.0 * 1024.0 * 1024.0)));
    clear_mem_cache();
    displayln_info(
        fname + ssprintf(": after clear mem_cache (current total %.3lf (GB))",
                         (RealD)ms.total() / (1024.0 * 1024.0 * 1024.0)));
  }
  {
    TIMER_FLOPS("alloc_mem-alloc");
    timer.flops += min_size;
    void* ptr = NULL;
    ptr = alloc_mem_alloc_impl(size, mem_type, min_alignment);
    qassert(ptr == NULL or ((uintptr_t)ptr) % min_alignment == 0);
    return ptr;
  }
}

void free_mem(void* ptr, const Long min_size, const MemType mem_type)
{
  free_mem(ptr, min_size, mem_type, 1);
}

void free_mem(void* ptr, const Long min_size, const MemType mem_type,
              const Long min_alignment)
{
  TIMER_FLOPS("free_mem");
  timer.flops += min_size;
  const Long chunk_size = get_mem_chunk_size(mem_type);
  const Long size = get_chunked_mem_size(chunk_size, min_size);
  MemCache& cache = get_mem_cache(mem_type);
  cache.add(ptr, size, min_alignment);
}

void copy_mem(void* dst, const MemType mem_type_dst, const void* src,
              const MemType mem_type_src, const Long size)
{
#ifdef QLAT_USE_ACC
  const MemType eff_mem_type_dst = get_eff_mem_type(mem_type_dst);
  const MemType eff_mem_type_src = get_eff_mem_type(mem_type_src);
  qacc_Error err = qacc_ErrorUnknown;
  if (eff_mem_type_src == MemType::Uvm or eff_mem_type_dst == MemType::Uvm) {
    bool with_cpu = false;
    if (size < 4096) {
      if (eff_mem_type_src == MemType::Uvm and
          eff_mem_type_dst == MemType::Uvm) {
        with_cpu = true;
      }
      if (eff_mem_type_src == MemType::Cpu or
          eff_mem_type_dst == MemType::Cpu) {
        with_cpu = true;
      }
    }
    if (with_cpu == true) {
      TIMER_FLOPS("copy_mem(Cpu<-Cpu)");
      timer.flops += size;
      std::memcpy(dst, src, size);
      err = qacc_Success;
    }
    if (with_cpu == false) {
      TIMER_FLOPS("copy_mem(Uvm)");
      timer.flops += size;
      err = qacc_Memcpy(dst, src, size, qacc_MemcpyDefault);
    }
  } else if (eff_mem_type_src == MemType::Cpu) {
    if (eff_mem_type_dst == MemType::Cpu) {
      TIMER_FLOPS("copy_mem(Cpu<-Cpu)");
      timer.flops += size;
      std::memcpy(dst, src, size);
      err = qacc_Success;
    } else if (eff_mem_type_dst == MemType::Acc) {
      TIMER_FLOPS("copy_mem(Acc<-Cpu)");
      timer.flops += size;
      err = qacc_Memcpy(dst, src, size, qacc_MemcpyHostToDevice);
    } else {
      qassert(false);
    }
  } else if (eff_mem_type_src == MemType::Acc) {
    if (eff_mem_type_dst == MemType::Cpu) {
      TIMER_FLOPS("copy_mem(Cpu<-Acc)");
      timer.flops += size;
      err = qacc_Memcpy(dst, src, size, qacc_MemcpyDeviceToHost);
    } else if (eff_mem_type_dst == MemType::Acc) {
      TIMER_FLOPS("copy_mem(Acc<-Acc)");
      timer.flops += size;
      err = qacc_Memcpy(dst, src, size, qacc_MemcpyDeviceToDevice);
    } else {
      qassert(false);
    }
  } else {
    qassert(false);
  }
  qacc_ErrCheck(qacc_DeviceSynchronize());  // HIP XNACK need barrier...
  if (qacc_Success != err) {
    qerr(ssprintf("mem copy : ACC error '%s' (%d) after qacc_Malloc.",
                  qacc_GetErrorString(err), err));
  }
#else
  TIMER_FLOPS("copy_mem");
  timer.flops += size;
  (void)mem_type_dst;
  (void)mem_type_src;
  std::memcpy(dst, src, size);
#endif
}

void set_mem(void* ptr, const Int v, const Long size, const MemType mem_type)
{
#ifdef QLAT_USE_ACC
  const MemType eff_mem_type = get_eff_mem_type(mem_type);
  if (eff_mem_type == MemType::Cpu) {
    memset(ptr, v, size);
  } else {
    qacc_Error err = qacc_ErrorUnknown;
    err = qacc_Memset(ptr, v, size);
    if (qacc_Success != err) {
      qerr(ssprintf("set_zero(vector): ACC error '%s' (%d) after qacc_Memset.",
                    qacc_GetErrorString(err), err));
    }
  }
#else
  (void)mem_type;
  memset(ptr, v, size);
#endif
}

}  // namespace qlat
