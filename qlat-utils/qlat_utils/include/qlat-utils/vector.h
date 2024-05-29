#pragma once

#include <qlat-utils/config.h>
#include <qlat-utils/handle.h>
#include <qlat-utils/utils.h>
#include <unistd.h>

#include <cstdlib>
#include <cstring>
#include <unordered_map>

#ifdef QLAT_USE_MALLOC_STATS
#include <malloc.h>
#endif

namespace qlat
{  //

API inline Long& get_alignment()
// qlat parameter
//
// Should NOT change in the middle of the run.
{
  static Long alignment = 256;
  return alignment;
}

inline Long get_aligned_mem_size(const Long alignment, const Long min_size)
{
  const Long n_elem = 1 + (min_size - 1) / alignment;
  const Long size = n_elem * alignment;
  return size;
}

API inline Long& get_mem_cache_max_size(const bool is_acc = false)
// qlat parameter
// unit in MB
{
  static Long max_size =
      get_env_long_default("q_mem_cache_max_size", 512) * 1024L * 1024L;
  static Long max_size_acc = get_env_long_default("q_mem_cache_acc_max_size",
                                                  max_size / (1024L * 1024L)) *
                             1024L * 1024L;
  if (is_acc) {
    return max_size_acc;
  } else {
    return max_size;
  }
}

API inline Long& get_alloc_mem_max_size()
// qlat parameter
// unit in MB
{
  static Long max_size =
      get_env_long_default("q_alloc_mem_max_size", 256L * 1024L) * 1024L *
      1024L;
  return max_size;
}

struct MemoryStats {
  Long alloc;
  Long alloc_acc;
  Long cache;
  Long cache_acc;
  //
  MemoryStats() { init(); }
  //
  void init()
  {
    alloc = 0;
    alloc_acc = 0;
    cache = 0;
    cache_acc = 0;
  }
  //
  Long total()
  // total memory include the size of the cache
  {
    return alloc + alloc_acc;
  }
};

API inline MemoryStats& get_mem_stats()
{
  static MemoryStats ms;
  return ms;
}

inline void* alloc_mem_alloc(const Long size)
{
  TIMER_FLOPS("alloc_mem_alloc");
  timer.flops += size;
  static MemoryStats& ms = get_mem_stats();
  ms.alloc += size;
#if defined QLAT_NO_ALIGNED_ALLOC
  return malloc(size);
#else
  const Long alignment = get_alignment();
  return aligned_alloc(alignment, size);
#endif
}

inline void* alloc_mem_alloc_acc(const Long size)
{
#ifdef QLAT_USE_ACC
  TIMER_FLOPS("alloc_mem_alloc_acc");
  timer.flops += size;
  static MemoryStats& ms = get_mem_stats();
  ms.alloc_acc += size;
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
  return alloc_mem_alloc(size);
#endif
}

inline void free_mem_free(void* ptr, const Long size)
{
  TIMER_FLOPS("free_mem_free");
  timer.flops += size;
  static MemoryStats& ms = get_mem_stats();
  ms.alloc -= size;
  free(ptr);
}

inline void free_mem_free_acc(void* ptr, const Long size)
{
#ifdef QLAT_USE_ACC
  TIMER_FLOPS("free_mem_free_acc");
  timer.flops += size;
  static MemoryStats& ms = get_mem_stats();
  ms.alloc_acc -= size;
  qacc_Error err = qacc_Free(ptr);
  if (qacc_Success != err) {
    if (qacc_ErrorCudartUnloading != err) {
      qerr(fname + ssprintf(": Cuda error '%s' (%d) after qacc_Free.",
                            qacc_GetErrorString(err), err));
    }
  }
#else
  free_mem_free(ptr, size);
#endif
}

struct API MemCache {
  bool is_acc;
  Long mem_cache_size;
  Long mem_cache_max_size;
  std::unordered_multimap<Long, void*> db;
  //
  MemCache(const bool is_acc_ = false, const Long mem_cache_max_size_ = -1)
  {
    is_acc = is_acc_;
    mem_cache_max_size = mem_cache_max_size_;
    if (mem_cache_max_size < 0) {
      mem_cache_max_size = get_mem_cache_max_size(is_acc);
    }
    mem_cache_size = 0;
  }
  ~MemCache() { gc(); }
  //
  void add(void* ptr, const Long size)
  {
    qassert(size > 0);
    mem_cache_size += size;
    static MemoryStats& ms = get_mem_stats();
    if (is_acc) {
      ms.cache_acc += size;
    } else {
      ms.cache += size;
    }
    std::pair<Long, void*> p(size, ptr);
    db.insert(p);
    if (mem_cache_size > mem_cache_max_size) {
      gc();
    }
  }
  //
  void* del(const Long size)
  {
    auto iter = db.find(size);
    if (iter == db.end()) {
      return NULL;
    } else {
      mem_cache_size -= size;
      static MemoryStats& ms = get_mem_stats();
      if (is_acc) {
        ms.cache_acc -= size;
      } else {
        ms.cache -= size;
      }
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
    static MemoryStats& ms = get_mem_stats();
    if (is_acc) {
      ms.cache_acc -= mem_cache_size;
    } else {
      ms.cache -= mem_cache_size;
    }
    for (auto iter = db.cbegin(); iter != db.cend(); ++iter) {
      const Long size = iter->first;
      void* ptr = iter->second;
      qassert(ptr != NULL);
      if (is_acc) {
        free_mem_free_acc(ptr, size);
      } else {
        free_mem_free(ptr, size);
      }
      mem_cache_size -= size;
    }
    qassert(mem_cache_size == 0);
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

inline void* alloc_mem(const Long min_size, const bool is_acc = false)
{
  if (min_size <= 0) {
    return NULL;
  }
  TIMER_FLOPS("alloc_mem");
  timer.flops += min_size;
  const Long alignment = get_alignment();
  const Long size = get_aligned_mem_size(alignment, min_size);
  MemCache& cache = get_mem_cache(is_acc);
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
    if (is_acc) {
      ptr = alloc_mem_alloc_acc(size);
    } else {
      ptr = alloc_mem_alloc(size);
    }
    memset(ptr, 0, size);
    return ptr;
  }
}

inline void free_mem(void* ptr, const Long min_size, const bool is_acc = false)
{
  TIMER_FLOPS("free_mem");
  timer.flops += min_size;
  const Long alignment = get_alignment();
  const Long size = get_aligned_mem_size(alignment, min_size);
  MemCache& cache = get_mem_cache(is_acc);
  cache.add(ptr, size);
}

inline void displayln_malloc_stats()
{
#ifdef QLAT_USE_MALLOC_STATS
  malloc_stats();
#endif
}

// --------------------

template <class M>
struct API vector {
  // Avoid copy constructor when possible
  // (it is likely not what you think it is)
  // Only used in qacc macros, or if it is already a copy.
  //
  bool is_copy;  // do not free memory if is_copy=true
  bool is_acc;   // if place data on qacc_MallocManaged memory (default false)
  Vector<M> v;
  //
  vector()
  {
    // TIMER("vector::vector()")
    qassert(v.p == NULL);
    is_copy = false;
    is_acc = false;
  }
  vector(const vector<M>& vp)
  {
    // TIMER("vector::vector(&)")
#ifndef QLAT_USE_ACC
    qassert(vp.is_copy);
#endif
    is_copy = true;
    qassert(vp.is_acc);
    is_acc = vp.is_acc;
    v = vp.v;
  }
  vector(vector<M>&& vp) noexcept
  {
    // TIMER("vector::vector(&&)")
    is_copy = vp.is_copy;
    is_acc = vp.is_acc;
    v = vp.v;
    vp.is_copy = true;
  }
  vector(const Long size)
  {
    // TIMER("vector::vector(size)")
    qassert(v.p == NULL);
    is_copy = false;
    is_acc = false;
    resize(size);
  }
  vector(const Long size, const M& x)
  {
    // TIMER("vector::vector(size,x)")
    qassert(v.p == NULL);
    is_copy = false;
    is_acc = false;
    resize(size, x);
  }
  vector(const std::vector<M>& vp)
  {
    // TIMER("vector::vector(std::vector&)")
    is_copy = false;
    is_acc = false;
    *this = vp;
  }
  //
  ~vector()
  {
    // TIMER("vector::~vector()")
    if (not is_copy) {
      clear();
    }
  }
  //
  void init()
  // does not change is_acc
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
  void set_view(const Vector<M>& vec)
  // does not change is_acc
  {
    init();
    is_copy = true;
    v = vec;
  }
  void set_view(const vector<M>& vec)
  {
    init();
    is_copy = true;
    is_acc = vec.is_acc;
    v = vec.v;
  }
  //
  template <class N>
  void set_view_cast(const vector<N>& vec)
  {
    init();
    is_copy = true;
    is_acc = vec.is_acc;
    v.set_cast(vec.v);
  }
  //
  void resize(const Long size)
  {
    qassert(not is_copy);
    qassert(0 <= size);
    if (v.p == NULL) {
      v.p = (M*)alloc_mem(size * sizeof(M), is_acc);
      v.n = size;
    } else if (v.n != size) {
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
  void resize(const Long size, const M& x)
  {
    qassert(not is_copy);
    qassert(0 <= size);
    if (v.p == NULL) {
      v.p = (M*)alloc_mem(size * sizeof(M), is_acc);
      v.n = size;
      for (Long i = 0; i < v.n; ++i) {
        v[i] = x;
      }
    } else if (v.n != size) {
      vector<M> vp;
      vp.set_acc(is_acc);
      vp.v = v;
      v.p = (M*)alloc_mem(size * sizeof(M), is_acc);
      v.n = size;
      if (size <= vp.v.n) {
        std::memcpy(v.p, vp.v.p, v.n * sizeof(M));
      } else {
        std::memcpy(v.p, vp.v.p, vp.v.n * sizeof(M));
        for (Long i = vp.v.n; i < v.n; ++i) {
          v[i] = x;
        }
      }
    }
  }
  //
  vector<M>& operator=(const vector<M>& vp)
  {
    // TIMER("vector::operator=(&)");
    qassert(not is_copy);
    clear();
    resize(vp.size());
    for (Long i = 0; i < v.n; ++i) {
      v[i] = vp[i];
    }
    return *this;
  }
  vector<M>& operator=(vector<M>&& vp) noexcept
  {
    // TIMER("vector::operator=(&&)");
    is_copy = vp.is_copy;
    is_acc = vp.is_acc;
    v = vp.v;
    vp.is_copy = true;
    return *this;
  }
  vector<M>& operator=(const std::vector<M>& vp)
  {
    // TIMER("vector::operator=(std::vector&)");
    qassert(not is_copy);
    clear();
    resize(vp.size());
    for (Long i = 0; i < v.n; ++i) {
      v[i] = vp[i];
    }
    return *this;
  }
  //
  qacc const M& operator[](const Long i) const { return v[i]; }
  qacc M& operator[](const Long i) { return v[i]; }
  //
  qacc Long size() const { return v.size(); }
  //
  qacc M* data() { return v.data(); }
  qacc const M* data() const { return v.data(); }
};

template <class M>
struct API vector_acc : vector<M> {
  // Avoid copy constructor when possible
  // (it is likely not what you think it is)
  // Only used in qacc macros.
  //
  using vector<M>::v;
  using vector<M>::is_copy;
  using vector<M>::is_acc;
  using vector<M>::resize;
  //
  vector_acc()
  {
    // TIMER("vector_acc::vector_acc()");
    qassert(v.p == NULL);
    is_copy = false;
    is_acc = true;
  }
  vector_acc(const vector_acc<M>& vp)
  {
    // TIMER("vector_acc::vector_acc(&)");
#ifndef QLAT_USE_ACC
    qassert(vp.is_copy);
#endif
    is_copy = true;
    qassert(vp.is_acc);
    is_acc = vp.is_acc;
    v = vp.v;
  }
  vector_acc(vector_acc<M>&& vp) noexcept
  {
    // TIMER("vector_acc::vector_acc(&&)")
    // qassert(vp.is_acc);
    is_copy = vp.is_copy;
    is_acc = vp.is_acc;
    v = vp.v;
    vp.is_copy = true;
  }
  vector_acc(const Long size)
  {
    // TIMER("vector_acc::vector_acc(size)");
    qassert(v.p == NULL);
    is_copy = false;
    is_acc = true;
    resize(size);
  }
  vector_acc(const Long size, const M& x)
  {
    // TIMER("vector_acc::vector_acc(size,x)");
    qassert(v.p == NULL);
    is_copy = false;
    is_acc = true;
    resize(size, x);
  }
  vector_acc(const std::vector<M>& vp)
  {
    // TIMER("vector_acc::vector_acc(std::vector&)");
    is_copy = false;
    is_acc = true;
    *this = vp;
  }
  //
  vector_acc<M>& operator=(const vector_acc<M>& vp)
  {
    // TIMER("vector_acc::operator=(&)");
    vector<M>::operator=(vp);
    return *this;
  }
  vector_acc<M>& operator=(const vector<M>& vp)
  {
    // TIMER("vector_acc::operator=(&)");
    vector<M>::operator=(vp);
    return *this;
  }
  vector_acc<M>& operator=(vector_acc<M>&& vp) noexcept
  {
    // TIMER("vector_acc::operator=(&&)");
    vector<M>::operator=(std::move(vp));
    return *this;
  }
  vector_acc<M>& operator=(vector<M>&& vp) noexcept
  {
    // TIMER("vector_acc::operator=(&&)");
    vector<M>::operator=(std::move(vp));
    return *this;
  }
  vector_acc<M>& operator=(const std::vector<M>& vp)
  {
    // TIMER("vector_acc::operator=(std::vector&)");
    vector<M>::operator=(vp);
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

// --------------------

template <class M>
qacc bool operator==(const vector<M>& v1, const vector<M>& v2)
{
  if (v1.size() != v2.size()) {
    return false;
  }
  const int cmp = std::memcmp(v1.data(), v2.data(), v1.size() * sizeof(M));
  return cmp == 0;
}

template <class M>
qacc bool operator==(const vector_acc<M>& v1, const vector_acc<M>& v2)
{
  if (v1.size() != v2.size()) {
    return false;
  }
  const int cmp = std::memcmp(v1.data(), v2.data(), v1.size() * sizeof(M));
  return cmp == 0;
}

template <class M>
qacc bool operator!=(const vector<M>& v1, const vector<M>& v2)
{
  return not(v1 == v2);
}

template <class M>
qacc bool operator!=(const vector_acc<M>& v1, const vector_acc<M>& v2)
{
  return not(v1 == v2);
}

// --------------------

template <class M>
struct IsDataVectorType<vector<M>> {
  using DataType = M;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsDataValueType<DataType>::ElementaryType;
  static constexpr bool value = is_data_value_type<DataType>();
};

template <class M>
struct IsDataVectorType<vector_acc<M>> {
  using DataType = M;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsDataValueType<DataType>::ElementaryType;
  static constexpr bool value = is_data_value_type<DataType>();
};

template <class M>
qacc Vector<M> get_data(const vector<M>& v)
{
  return v.v;
}

template <class M>
qacc Vector<M> get_data(const vector_acc<M>& v)
{
  return get_data((const vector<M>&)v);
}

// --------------------

template <class M>
struct API box {
  //
  // like a one element vector
  //
  // Avoid copy constructor when possible
  // (it is likely not be what you think it is)
  // Only used in qacc macros, or if it is already a copy.
  //
  bool is_copy;  // do not free memory if is_copy=true
  bool is_acc;   // if place data on qacc_MallocManaged memory (default false)
  Handle<M> v;
  //
  box()
  {
    // TIMER("box::box()");
    qassert(v.p == NULL);
    is_copy = false;
    is_acc = false;
  }
  box(const box<M>& vp)
  {
    // TIMER("box::box(&)");
#ifndef QLAT_USE_ACC
    qassert(vp.is_copy);
#endif
    is_copy = true;
    qassert(vp.is_acc);
    is_acc = vp.is_acc;
    v = vp.v;
  }
  box(box<M>&& vp) noexcept
  {
    // TIMER("box::box(&&)");
    is_copy = vp.is_copy;
    is_acc = vp.is_acc;
    v = vp.v;
    vp.is_copy = true;
  }
  box(const M& x)
  {
    // TIMER("box::box(x)");
    qassert(v.p == NULL);
    is_copy = false;
    is_acc = false;
    set(x);
  }
  //
  ~box()
  {
    // TIMER("box::~box()");
    if (not is_copy) {
      clear();
    }
  }
  //
  void init()
  // does not change is_acc
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
  void set_view(const M& x)
  // does not change is_acc
  {
    init();
    is_copy = true;
    v.p = &x;
  }
  void set_view(const Handle<M> h)
  // does not change is_acc
  {
    init();
    is_copy = true;
    v = h;
  }
  void set_view(const box<M>& b)
  {
    init();
    is_copy = true;
    is_acc = b.is_acc;
    v = b.v;
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
  box<M>& operator=(const box<M>& vp)
  {
    // TIMER("box::operator=(&)");
    qassert(not is_copy);
    set(vp());
    return *this;
  }
  box<M>& operator=(box<M>&& vp) noexcept
  {
    // TIMER("box::operator=(&&)");
    is_copy = vp.is_copy;
    is_acc = vp.is_acc;
    v = vp.v;
    vp.is_copy = true;
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
  // Only used in qacc macros, or if it is already a copy.
  //
  using box<M>::v;
  using box<M>::is_copy;
  using box<M>::is_acc;
  using box<M>::set;
  //
  box_acc()
  {
    // TIMER("box_acc::box_acc()");
    qassert(v.p == NULL);
    is_copy = false;
    is_acc = true;
  }
  box_acc(const box_acc<M>& vp)
  {
    // TIMER("box::box(&)");
#ifndef QLAT_USE_ACC
    qassert(vp.is_copy);
#endif
    is_copy = true;
    qassert(vp.is_acc);
    is_acc = vp.is_acc;
    v = vp.v;
  }
  box_acc(box_acc<M>&& vp) noexcept
  {
    // TIMER("box_acc::box_acc(&&)");
    // qassert(vp.is_acc);
    is_copy = vp.is_copy;
    is_acc = vp.is_acc;
    v = vp.v;
    vp.is_copy = true;
  }
  box_acc(const M& x)
  {
    // TIMER("box_acc::box_acc(x)");
    qassert(v.p == NULL);
    is_copy = false;
    is_acc = true;
    set(x);
  }
  //
  box_acc<M>& operator=(const box_acc<M>& vp)
  {
    // TIMER("box_acc::operator=(&)");
    box<M>::operator=(vp);
    return *this;
  }
  box_acc<M>& operator=(const box<M>& vp)
  {
    // TIMER("box_acc::operator=(&)");
    box<M>::operator=(vp);
    return *this;
  }
  box_acc<M>& operator=(box_acc<M>&& vp) noexcept
  {
    // TIMER("box_acc::operator=(&&)");
    box<M>::operator=(std::move(vp));
    return *this;
  }
  box_acc<M>& operator=(box<M>&& vp) noexcept
  {
    // TIMER("box_acc::operator=(&&)");
    box<M>::operator=(std::move(vp));
    return *this;
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

// --------------------

template <class M>
struct IsDataVectorType<box<M>> {
  using DataType = M;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsDataValueType<DataType>::ElementaryType;
  static constexpr bool value = is_data_value_type<DataType>();
};

template <class M>
struct IsDataVectorType<box_acc<M>> {
  using DataType = M;
  using BasicDataType = typename IsDataValueType<DataType>::BasicDataType;
  using ElementaryType = typename IsDataValueType<DataType>::ElementaryType;
  static constexpr bool value = is_data_value_type<DataType>();
};

template <class M>
qacc Vector<M> get_data(const box<M>& v)
{
  if (not v.null()) {
    return get_data(v());
  } else {
    return Vector<M>();
  }
}

template <class M>
qacc Vector<M> get_data(const box_acc<M>& v)
{
  return get_data((const box<M>&)v);
}

}  // namespace qlat
