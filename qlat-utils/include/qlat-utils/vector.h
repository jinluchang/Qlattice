#pragma once

#include <qlat-utils/config.h>
#include <qlat-utils/handle.h>
#include <qlat-utils/qutils.h>

#include <cstdlib>
#include <cstring>
#include <unordered_map>

#ifdef QLAT_USE_MALLOC_STATS
#include <malloc.h>
#endif

namespace qlat
{  //

API inline size_t& get_alignment()
// qlat parameter
//
// Should NOT change in the middle of the run.
{
  static size_t alignment = 256;
  return alignment;
}

inline size_t get_aligned_mem_size(const size_t alignment, const long min_size)
{
  const long n_elem = 1 + (min_size - 1) / alignment;
  const size_t size = n_elem * alignment;
  return size;
}

API inline size_t& get_mem_cache_max_size(const bool is_acc = false)
// qlat parameter
{
  static size_t max_size =
      get_env_long_default("q_mem_cache_max_size", 512) * 1024L * 1024L;
  static size_t max_size_acc =
      get_env_long_default("q_mem_cache_acc_max_size",
                           max_size / (1024L * 1024L)) *
      1024L * 1024L;
  if (is_acc) {
    return max_size_acc;
  } else {
    return max_size;
  }
}

struct API MemCache {
  bool is_acc;
  size_t mem_cache_size;
  std::unordered_multimap<size_t, void*> db;
  //
  MemCache(const bool is_acc_ = false)
  {
    is_acc = is_acc_;
    mem_cache_size = 0;
  }
  ~MemCache() { gc(); }
  //
  void add(void* ptr, const size_t size)
  {
    qassert(size > 0);
    mem_cache_size += size;
    std::pair<size_t, void*> p(size, ptr);
    db.insert(p);
    if (mem_cache_size > get_mem_cache_max_size(is_acc)) {
      gc();
    }
  }
  //
  void* del(const size_t size)
  {
    auto iter = db.find(size);
    if (iter == db.end()) {
      return NULL;
    } else {
      mem_cache_size -= size;
      void* ptr = iter->second;
      db.erase(iter);
      return ptr;
    }
  }
  //
  void gc()
  {
    if (mem_cache_size == 0) {
      return;
    }
    TIMER_FLOPS("MemCache::gc()");
    timer.flops += mem_cache_size;
    for (auto iter = db.cbegin(); iter != db.cend(); ++iter) {
      void* ptr = iter->second;
      qassert(ptr != NULL);
#ifdef QLAT_USE_ACC
      if (is_acc) {
        cudaError err = cudaFree(ptr);
        if (cudaSuccess != err) {
          displayln(fname + ssprintf(": Cuda error %s after cudaFree.",
                                     cudaGetErrorString(err)));
          qassert(err == cudaSuccess);
        }
      } else {
        free(ptr);
      }
#else
      free(ptr);
#endif
    }
    mem_cache_size = 0;
    db.clear();
  }
};

API inline MemCache& get_mem_cache(const bool is_acc = false)
{
  static MemCache cache(false);
  static MemCache cache_acc(true);
  if (is_acc) {
    return cache_acc;
  } else {
    return cache;
  }
}

inline void clear_mem_cache()
{
  TIMER_VERBOSE_FLOPS("clear_mem_cache");
  long total_bytes = 0;
  total_bytes += get_mem_cache(false).mem_cache_size;
  total_bytes += get_mem_cache(true).mem_cache_size;
  get_mem_cache(false).gc();
  get_mem_cache(true).gc();
  displayln_info(
      0, fname + ssprintf(": %ld bytes (%.3f GB) freed.", total_bytes,
                          (double)total_bytes / (1024.0 * 1024.0 * 1024.0)));
  timer.flops += total_bytes;
}

inline void* alloc_mem_alloc_no_acc(const long size)
{
  const size_t alignment = get_alignment();
#if defined QLAT_NO_ALIGNED_ALLOC
  return malloc(size);
#else
  return aligned_alloc(alignment, size);
#endif
}

inline void* alloc_mem(const long min_size, const bool is_acc = false)
{
  if (min_size <= 0) {
    return NULL;
  }
  TIMER_FLOPS("alloc_mem");
  timer.flops += min_size;
  const size_t alignment = get_alignment();
  const size_t size = get_aligned_mem_size(alignment, min_size);
  MemCache& cache = get_mem_cache(is_acc);
  void* ptr = cache.del(size);
  if (NULL != ptr) {
    return ptr;
  }
  {
    TIMER_FLOPS("alloc_mem-alloc");
    timer.flops += min_size;
    void* ptr = NULL;
#ifdef QLAT_USE_ACC
    if (is_acc) {
      cudaError err = cudaGetLastError();
      if (cudaSuccess != err) {
        displayln(fname + ssprintf(": Cuda error %s before cudaMallocManaged.",
                                   cudaGetErrorString(err)));
        qassert(err == cudaSuccess);
      }
      err = cudaMallocManaged(&ptr, size);
      if (cudaSuccess != err) {
        displayln(fname +
                  ssprintf(": Cuda error %s, min_size=%ld, size=%ld, ptr=%lX.",
                           cudaGetErrorString(err), min_size, size, ptr));
        usleep((useconds_t)(10.0 * 1.0e6));
        qassert(err == cudaSuccess);
      }
    } else {
      ptr = alloc_mem_alloc_no_acc(size);
    }
#else
    ptr = alloc_mem_alloc_no_acc(size);
#endif
    memset(ptr, 0, size);
    return ptr;
  }
}

inline void free_mem(void* ptr, const long min_size, const bool is_acc = false)
{
  TIMER_FLOPS("free_mem");
  timer.flops += min_size;
  const size_t alignment = get_alignment();
  const size_t size = get_aligned_mem_size(alignment, min_size);
  MemCache& cache = get_mem_cache(is_acc);
  cache.add(ptr, size);
}

inline void displayln_malloc_stats()
{
#ifdef QLAT_USE_MALLOC_STATS
  malloc_stats();
#endif
}

template <class M>
struct API vector {
  // Avoid copy constructor when possible
  // (it is likely not what you think it is)
  //
  bool is_copy;  // do not free memory if is_copy=true
  bool is_acc;   // if place data on cudaMallocManaged memory (default false)
  Vector<M> v;
  //
  vector()
  {
    qassert(v.p == NULL);
    is_copy = false;
    is_acc = false;
  }
  vector(const vector<M>& vp)
  {
#ifndef QLAT_USE_ACC
    qassert(false);
#endif
    is_copy = true;
    qassert(vp.is_acc);
    is_acc = vp.is_acc;
    v = vp.v;
  }
  vector(vector<M>&& vp) noexcept
  {
    is_copy = vp.is_copy;
    is_acc = vp.is_acc;
    v = vp.v;
    vp.is_copy = true;
  }
  vector(const long size)
  {
    qassert(v.p == NULL);
    is_copy = false;
    is_acc = false;
    resize(size);
  }
  vector(const long size, const M& x)
  {
    qassert(v.p == NULL);
    is_copy = false;
    is_acc = false;
    resize(size, x);
  }
  vector(const std::vector<M>& vp)
  {
    is_copy = false;
    is_acc = false;
    *this = vp;
  }
  //
  ~vector()
  {
    if (not is_copy) {
      clear();
    }
  }
  //
  void init()
  {
    if (not is_copy) {
      clear();
    }
    is_copy = false;
  }
  //
  void clear()
  {
    qassert(not is_copy);
    if (v.p != NULL) {
      free_mem(v.p, v.n * sizeof(M), is_acc);
    }
    v = Vector<M>();
    qassert(v.p == NULL);
  }
  //
  void set_acc(const bool is_acc_)
  {
    qassert(not is_copy);
    if (is_acc == is_acc_) {
      return;
    }
    if (NULL == v.p) {
      is_acc = is_acc_;
      return;
    }
    vector<M> vec;
    vec.set_acc(is_acc_);
    vec = *this;
    swap(vec);
    qassert(is_acc == is_acc_);
  }
  //
  qacc void swap(vector<M>& x)
  {
    qassert(not is_copy);
    qassert(not x.is_copy);
    bool tb = x.is_acc;
    x.is_acc = is_acc;
    is_acc = tb;
    Vector<M> t = v;
    v = x.v;
    x.v = t;
  }
  //
  void resize(const long size)
  {
    qassert(not is_copy);
    qassert(0 <= size);
    if (v.p == NULL) {
      v.p = (M*)alloc_mem(size * sizeof(M), is_acc);
      v.n = size;
    } else {
      vector<M> vp;
      vp.set_acc(is_acc);
      vp.v = v;
      v.p = (M*)alloc_mem(size * sizeof(M), is_acc);
      v.n = size;
      if (size <= vp.v.n) {
        std::memcpy((void*)v.p, (void*)vp.v.p, size * sizeof(M));
      } else {
        std::memcpy((void*)v.p, (void*)vp.v.p, vp.v.n * sizeof(M));
      }
    }
  }
  void resize(const long size, const M& x)
  {
    qassert(not is_copy);
    qassert(0 <= size);
    if (v.p == NULL) {
      v.p = (M*)alloc_mem(size * sizeof(M), is_acc);
      v.n = size;
      for (long i = 0; i < v.n; ++i) {
        v[i] = x;
      }
    } else {
      vector<M> vp;
      vp.set_acc(is_acc);
      vp.v = v;
      v.p = (M*)alloc_mem(size * sizeof(M), is_acc);
      v.n = size;
      if (size <= vp.v.n) {
        std::memcpy(v.p, vp.v.p, size * sizeof(M));
      } else {
        std::memcpy(v.p, vp.v.p, vp.v.n * sizeof(M));
        for (long i = size; i < v.n; ++i) {
          v[i] = x;
        }
      }
    }
  }
  //
  const vector<M>& operator=(const vector<M>& vp)
  {
    qassert(not is_copy);
    clear();
    resize(vp.size());
    for (long i = 0; i < v.n; ++i) {
      v[i] = vp[i];
    }
    return *this;
  }
  const vector<M>& operator=(const std::vector<M>& vp)
  {
    qassert(not is_copy);
    clear();
    resize(vp.size());
    for (long i = 0; i < v.n; ++i) {
      v[i] = vp[i];
    }
    return *this;
  }
  //
  qacc const M& operator[](const long i) const { return v[i]; }
  qacc M& operator[](const long i) { return v[i]; }
  //
  qacc long size() const { return v.size(); }
  //
  qacc M* data() { return v.data(); }
  qacc const M* data() const { return v.data(); }
};

template <class M>
struct API vector_acc : vector<M> {
  // Avoid copy constructor when possible
  // (it is likely not what you think it is)
  //
  using vector<M>::v;
  using vector<M>::is_copy;
  using vector<M>::is_acc;
  using vector<M>::resize;
  using vector<M>::operator=;
  //
  vector_acc()
  {
    qassert(v.p == NULL);
    is_copy = false;
    is_acc = true;
  }
  vector_acc(const long size)
  {
    qassert(v.p == NULL);
    is_copy = false;
    is_acc = true;
    resize(size);
  }
  vector_acc(const long size, const M& x)
  {
    qassert(v.p == NULL);
    is_copy = false;
    is_acc = true;
    resize(size, x);
  }
  vector_acc(const std::vector<M>& vp)
  {
    is_copy = false;
    is_acc = true;
    *this = vp;
  }
  vector_acc(const vector<M>& vp)
  {
#ifndef QLAT_USE_ACC
    qassert(false);
#endif
    is_copy = true;
    qassert(vp.is_acc);
    is_acc = vp.is_acc;
    v = vp.v;
  }
  vector_acc(vector<M>&& vp) noexcept
  {
    qassert(vp.is_acc);
    is_copy = vp.is_copy;
    is_acc = vp.is_acc;
    v = vp.v;
    vp.is_copy = true;
  }
  ~vector_acc() {}
  //
  const vector_acc<M>& operator=(const vector<M>& vp)
  {
    this->vector<M>::operator=(vp);
    return *this;
  }
  const vector_acc<M>& operator=(const std::vector<M>& vp)
  {
    this->vector<M>::operator=(vp);
    return *this;
  }
};

template <class M>
void clear(vector<M>& v)
{
  v.clear();
}

template <class M>
qacc void qswap(vector<M>& v1, vector<M>& v2)
{
  v1.swap(v2);
}

template <class M>
qacc Vector<M> get_data(const vector<M>& v)
{
  return v.v;
}

template <class M>
qacc void set_zero(vector<M>& v)
{
  set_zero(v.v);
}

template <class T>
qacc double qnorm(const vector<T>& mm)
{
  double sum = 0.0;
  for (long i = 0; i < mm.size(); ++i) {
    sum += qnorm(mm[i]);
  }
  return sum;
}

template <class M>
struct API box {
  //
  // like a one element vector
  //
  // Avoid copy constructor when possible
  // (it is likely not be what you think it is)
  //
  bool is_copy;  // do not free memory if is_copy=true
  bool is_acc;   // if place data on cudaMallocManaged memory (default false)
  Handle<M> v;
  //
  box()
  {
    qassert(v.p == NULL);
    is_copy = false;
    is_acc = false;
  }
  box(const box<M>& vp)
  {
#ifndef QLAT_USE_ACC
    qassert(false);
#endif
    is_copy = true;
    qassert(vp.is_acc);
    is_acc = vp.is_acc;
    v = vp.v;
  }
  box(box<M>&& vp) noexcept
  {
    is_copy = vp.is_copy;
    is_acc = vp.is_acc;
    v = vp.v;
    vp.is_copy = true;
  }
  box(const M& x)
  {
    qassert(v.p == NULL);
    is_copy = false;
    is_acc = false;
    set(x);
  }
  //
  ~box()
  {
    if (not is_copy) {
      clear();
    }
  }
  //
  void init()
  {
    if (not is_copy) {
      clear();
    }
    is_copy = false;
  }
  //
  void clear()
  {
    qassert(not is_copy);
    if (v.p != NULL) {
      free_mem(v.p, sizeof(M), is_acc);
    }
    v = Handle<M>();
    qassert(v.p == NULL);
  }
  //
  void set_acc(const bool is_acc_)
  {
    qassert(not is_copy);
    if (is_acc == is_acc_) {
      return;
    }
    if (NULL == v.p) {
      is_acc = is_acc_;
      return;
    }
    box<M> b;
    b.set_acc(is_acc_);
    b = *this;
    swap(b);
    qassert(is_acc == is_acc_);
  }
  //
  qacc void swap(box<M>& x)
  {
    qassert(not is_copy);
    qassert(not x.is_copy);
    bool tb = x.is_acc;
    x.is_acc = is_acc;
    is_acc = tb;
    Handle<M> t = v;
    v = x.v;
    x.v = t;
  }
  //
  void set(const M& x)
  {
    qassert(not is_copy);
    if (v.p == NULL) {
      v.p = (M*)alloc_mem(sizeof(M), is_acc);
    }
    v() = x;
  }
  //
  const box<M>& operator=(const box<M>& vp)
  {
    qassert(not is_copy);
    set(vp());
    return *this;
  }
  //
  qacc const M& operator()() const { return v(); }
  qacc M& operator()() { return v(); }
  //
  qacc bool null() const { return v.null(); }
};

template <class M>
struct API box_acc : box<M> {
  //
  // like a one element vector_acc
  //
  // Avoid copy constructor when possible
  // (it is likely not be what you think it is)
  //
  using box<M>::v;
  using box<M>::is_copy;
  using box<M>::is_acc;
  using box<M>::set;
  //
  box_acc()
  {
    qassert(v.p == NULL);
    is_copy = false;
    is_acc = true;
  }
  box_acc(const M& x)
  {
    qassert(v.p == NULL);
    is_copy = false;
    is_acc = true;
    set(x);
  }
};

template <class M>
void clear(box<M>& v)
{
  v.clear();
}

template <class M>
qacc void qswap(box<M>& v1, box<M>& v2)
{
  v1.swap(v2);
}

template <class M>
qacc Vector<M> get_data(const box<M>& v)
{
  if (not v.null()) {
    return get_data_one_elem(v());
  } else {
    return Vector<M>();
  }
}

template <class M>
qacc void set_zero(box<M>& v)
{
  if (not v.null()) {
    set_zero(v());
  }
}

}  // namespace qlat
