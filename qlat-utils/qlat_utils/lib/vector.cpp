#include <qlat-utils/vector.h>

namespace qlat
{  //

static void* alloc_mem_aligned(const Long size, const MemType mem_type)
{
  void* ptr = NULL;
#if defined QLAT_NO_ALIGNED_ALLOC
  ptr = malloc(size);
#else
  Long alignment = get_alignment(mem_type);
  ptr = std::aligned_alloc(alignment, size);
  while (ptr == NULL) {
    alignment /= 2;
    if (alignment == 0) {
      ptr = malloc(size);
      break;
    } else {
      ptr = std::aligned_alloc(alignment, size);
    }
  }
#endif
  return ptr;
}

void* alloc_mem_alloc(const Long size, const MemType mem_type)
{
  TIMER_FLOPS("alloc_mem_alloc");
  timer.flops += size;
  static MemoryStats& ms = get_mem_stats();
  ms.alloc[static_cast<Int>(mem_type)] += size;
  void* ptr = NULL;
#ifdef QLAT_USE_ACC
  const MemType eff_mem_type = get_eff_mem_type(mem_type);
  if (eff_mem_type == MemType::Cpu) {
    TIMER_FLOPS("alloc_mem_alloc(MemType::Cpu)");
    timer.flops += size;
    ptr = alloc_mem_aligned(size, mem_type);
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
  ptr = alloc_mem_aligned(size, mem_type);
  qassert(ptr != NULL);
  set_mem(ptr, 0, size, mem_type);
#endif
  return ptr;
}

void free_mem_free(void* ptr, const Long size, const MemType mem_type)
{
  TIMER_FLOPS("free_mem_free");
  timer.flops += size;
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
  const Long chunk_size = get_mem_chunk_size(mem_type);
  const Long size = get_chunked_mem_size(chunk_size, min_size);
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
    ptr = alloc_mem_alloc(size, mem_type);
    return ptr;
  }
}

void free_mem(void* ptr, const Long min_size, const MemType mem_type)
{
  TIMER_FLOPS("free_mem");
  timer.flops += min_size;
  const Long chunk_size = get_mem_chunk_size(mem_type);
  const Long size = get_chunked_mem_size(chunk_size, min_size);
  MemCache& cache = get_mem_cache(mem_type);
  cache.add(ptr, size);
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
    if(size < 4096){
      if(eff_mem_type_src == MemType::Uvm and eff_mem_type_dst == MemType::Uvm){with_cpu = true;}
      if(eff_mem_type_src == MemType::Cpu or  eff_mem_type_dst == MemType::Cpu){with_cpu = true;}
    }   

    if(with_cpu == true){
      TIMER_FLOPS("copy_mem(Cpu<-Cpu)");
      timer.flops += size;

      std::memcpy(dst, src, size);
      err = qacc_Success;
    }   
    if(with_cpu == false)
    {   
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
