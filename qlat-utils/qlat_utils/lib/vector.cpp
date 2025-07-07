#include <qlat-utils/vector.h>

namespace qlat
{  //

std::string show(const MemType mem_type)
{
  if (mem_type == MemType::Cpu) {
    return "cpu";
  } else if (mem_type == MemType::Acc) {
    return "acc";
  } else if (mem_type == MemType::Comm) {
    return "comm";
  } else if (mem_type == MemType::Uvm) {
    return "uvm";
  } else {
    qassert(false);
    return "";
  }
}

MemType read_mem_type(const std::string& mem_type_str)
{
  if (mem_type_str == "cpu") {
    return MemType::Cpu;
  } else if (mem_type_str == "acc") {
    return MemType::Acc;
  } else if (mem_type_str == "comm") {
    return MemType::Comm;
  } else if (mem_type_str == "uvm") {
    return MemType::Uvm;
  } else {
    qassert(false);
    return MemType::Cpu;
  }
}

void* alloc_mem_aligned(const Long size)
{
#if defined QLAT_NO_ALIGNED_ALLOC
  return malloc(size);
#else
  const Long alignment = get_alignment();
  return std::aligned_alloc(alignment, size);
#endif
}

void* alloc_mem_alloc(const Long size, const MemType mem_type)
{
  TIMER_FLOPS("alloc_mem_alloc");
  timer.flops += size;
  static MemoryStats& ms = get_mem_stats();
  ms.alloc[static_cast<Int>(mem_type)] += size;
  if (mem_type == MemType::Cpu or mem_type == MemType::Comm) {
    return alloc_mem_aligned(size);
  } else if (mem_type == MemType::Acc) {
#ifdef QLAT_USE_ACC
    void* ptr = NULL;
    qacc_Error err = qacc_GetLastError();
    if (qacc_Success != err) {
      qerr(fname + ssprintf(": Cuda error '%s' before qacc_MallocManaged.",
                            qacc_GetErrorString(err)));
    }
    err = qacc_Malloc(&ptr, size);
    if (qacc_Success != err) {
      qerr(fname + ssprintf(": Cuda error '%s', size=%ld, ptr=%lX.",
                            qacc_GetErrorString(err), size, ptr));
    }
    return ptr;
#else
    return alloc_mem_aligned(size);
#endif
  } else if (mem_type == MemType::Uvm) {
#ifdef QLAT_USE_ACC
    void* ptr = NULL;
    qacc_Error err = qacc_GetLastError();
    if (qacc_Success != err) {
      qerr(fname + ssprintf(": Cuda error '%s' before qacc_MallocManaged.",
                            qacc_GetErrorString(err)));
    }
    err = qacc_MallocManaged(&ptr, size);
    if (qacc_Success != err) {
      qerr(fname + ssprintf(": Cuda error '%s', size=%ld, ptr=%lX.",
                            qacc_GetErrorString(err), size, ptr));
    }
    return ptr;
#else
    return alloc_mem_aligned(size);
#endif
  } else {
    qassert(false);
    return NULL;
  }
}

void free_mem_free(void* ptr, const Long size, const MemType mem_type)
{
  TIMER_FLOPS("free_mem_free");
  timer.flops += size;
  static MemoryStats& ms = get_mem_stats();
  ms.alloc[static_cast<Int>(mem_type)] -= size;
  if (mem_type == MemType::Cpu or mem_type == MemType::Comm) {
    free(ptr);
  } else if (mem_type == MemType::Acc or mem_type == MemType::Uvm) {
#ifdef QLAT_USE_ACC
    qacc_Error err = qacc_Free(ptr);
    if (qacc_Success != err) {
      if (qacc_ErrorCudartUnloading != err) {
        qerr(fname + ssprintf(": Cuda error '%s' (%d) after qacc_Free.",
                              qacc_GetErrorString(err), err));
      }
    }
#else
    free(ptr);
#endif
  } else {
    qassert(false);
  }
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

void MemCache::add(void* ptr, const Long size)
{
  qassert(size > 0);
  mem_cache_size += size;
  static MemoryStats& ms = get_mem_stats();
  ms.cache[static_cast<Int>(mem_type)] += size;
  std::pair<Long, void*> p(size, ptr);
  db.insert(p);
  if (mem_cache_size > mem_cache_max_size) {
    gc();
  }
}

void* MemCache::del(const Long size)
{
  auto iter = db.find(size);
  if (iter == db.end()) {
    return NULL;
  } else {
    mem_cache_size -= size;
    static MemoryStats& ms = get_mem_stats();
    ms.cache[static_cast<Int>(mem_type)] -= size;
    void* ptr = iter->second;
    db.erase(iter);
    return ptr;
  }
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
  if (min_size <= 0) {
    return NULL;
  }
  TIMER_FLOPS("alloc_mem");
  timer.flops += min_size;
  const Long alignment = get_alignment();
  const Long size = get_aligned_mem_size(alignment, min_size);
  MemCache& cache = get_mem_cache(mem_type);
  void* ptr = cache.del(size);
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
                         (double)min_size / (1024.0 * 1024.0 * 1024.0),
                         (double)ms.total() / (1024.0 * 1024.0 * 1024.0)));
    clear_mem_cache();
    displayln_info(
        fname + ssprintf(": after clear mem_cache (current total %.3lf (GB))",
                         (double)ms.total() / (1024.0 * 1024.0 * 1024.0)));
  }
  {
    TIMER_FLOPS("alloc_mem-alloc");
    timer.flops += min_size;
    void* ptr = NULL;
    ptr = alloc_mem_alloc(size, mem_type);
    memset(ptr, 0, size);
    return ptr;
  }
}

void free_mem(void* ptr, const Long min_size, const MemType mem_type)
{
  TIMER_FLOPS("free_mem");
  timer.flops += min_size;
  const Long alignment = get_alignment();
  const Long size = get_aligned_mem_size(alignment, min_size);
  MemCache& cache = get_mem_cache(mem_type);
  cache.add(ptr, size);
}

}  // namespace qlat
