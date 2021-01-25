#pragma once

#include <qutils/qutils-vec.h>
#include <qutils/qutils.h>

#include <cstdlib>
#include <cstring>
#include <unordered_map>

namespace qlat
{  //

inline size_t& get_alignment()
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

inline size_t get_mem_cache_max_size_default()
{
  TIMER_VERBOSE("get_mem_cache_max_size_default");
  const std::string n1 = get_env("Q_MEM_CACHE_MAX_SIZE");
  if (n1 != "") {
    const int n = read_long(n1);
    displayln_info(
        fname +
        ssprintf(": get_mem_cache_max_size() = %d MB via Q_MEM_CACHE_MAX_SIZE.",
                 n));
    return n;
  }
  const std::string n2 = get_env("q_mem_cache_max_size");
  if (n2 != "") {
    const int n = read_long(n2);
    displayln_info(
        fname +
        ssprintf(": get_mem_cache_max_size() = %d MB via q_mem_cache_max_size.",
                 n));
    return n;
  }
  displayln_info(fname + ssprintf(": get_mem_cache_max_size() = %d MB.", 512));
  return 512;
}

inline size_t& get_mem_cache_max_size()
// qlat parameter
{
  static size_t max_size = get_mem_cache_max_size_default() * 1024L * 1024L;
  return max_size;
}

struct MemCache {
  size_t mem_cache_size;
  std::unordered_multimap<size_t, void*> db;
  //
  MemCache() { mem_cache_size = 0; }
  ~MemCache() { gc(); }
  //
  void add(void* ptr, const size_t size)
  {
    mem_cache_size += size;
    std::pair<size_t, void*> p(size, ptr);
    db.insert(p);
    if (mem_cache_size > get_mem_cache_max_size()) {
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
    TIMER_VERBOSE_FLOPS("MemCache::gc()");
    timer.flops += mem_cache_size;
    for (auto iter = db.cbegin(); iter != db.cend(); ++iter) {
      void* ptr = iter->second;
#ifdef QLAT_USE_ACC
      cudaError_t code = cudaFree(ptr);
      qassert(code == cudaSuccess);
#else
      free(ptr);
#endif
    }
    mem_cache_size = 0;
    db.clear();
  }
};

inline MemCache& get_mem_cache()
{
  static MemCache cache;
  return cache;
}

inline void* alloc_mem(const long min_size)
{
  if (min_size <= 0) {
    return NULL;
  }
  TIMER_FLOPS("alloc_mem");
  timer.flops += min_size;
  const size_t alignment = get_alignment();
  const size_t size = get_aligned_mem_size(alignment, min_size);
  MemCache& cache = get_mem_cache();
  void* ptr = cache.del(size);
  if (NULL != ptr) {
    return ptr;
  }
  {
    TIMER_FLOPS("alloc_mem-alloc");
    timer.flops += min_size;
#ifdef QLAT_USE_ACC
    cudaError err = cudaGetLastError();
    if (cudaSuccess != err) {
      displayln(fname + ssprintf(": Cuda error %s before cudaMallocManaged.",
                                 cudaGetErrorString(err)));
      qassert(err == cudaSuccess);
    }
    void* ptr = NULL;
    err = cudaMallocManaged(&ptr, size);
    if (cudaSuccess != err) {
      displayln(fname +
                ssprintf(": Cuda error %s, min_size=%ld, size=%ld, ptr=%lX.",
                         cudaGetErrorString(err), min_size, size, ptr));
      usleep((useconds_t)(10.0 * 1.0e6));
      qassert(err == cudaSuccess);
    }
#else
    void* ptr = aligned_alloc(alignment, size);
#endif
    memset(ptr, 0, size);
    return ptr;
  }
}

inline void free_mem(void* ptr, const long min_size)
{
  TIMER_FLOPS("free_mem");
  timer.flops += min_size;
  const size_t alignment = get_alignment();
  const size_t size = get_aligned_mem_size(alignment, min_size);
  MemCache& cache = get_mem_cache();
  cache.add(ptr, size);
}

template <class M>
struct vector {
  // Avoid copy constructor when possible
  // (it is likely not be what you think it is)
  //
  bool is_copy;  // do not free memory if is_copy=true
  Vector<M> v;
  //
  vector()
  {
    qassert(v.p == NULL);
    is_copy = false;
  }
  vector(const vector<M>& vp)
  {
    is_copy = true;
    this->v = vp.v;
  }
  vector(const long size)
  {
    qassert(v.p == NULL);
    is_copy = false;
    resize(size);
  }
  vector(const long size, const M& x)
  {
    qassert(v.p == NULL);
    is_copy = false;
    resize(size, x);
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
      free_mem(v.p, v.n * sizeof(M));
    }
    v = Vector<M>();
    qassert(v.p == NULL);
  }
  //
  qacc void swap(vector<M>& x)
  {
    qassert(not is_copy);
    qassert(not x.is_copy);
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
      v.p = (M*)alloc_mem(size * sizeof(M));
      v.n = size;
    } else {
      vector<M> vp;
      vp.v = v;
      v.p = (M*)alloc_mem(size * sizeof(M));
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
      v.p = (M*)alloc_mem(size * sizeof(M));
      v.n = size;
      for (long i = 0; i < v.n; ++i) {
        v[i] = x;
      }
    } else {
      vector<M> vp;
      vp.v = v;
      v.p = (M*)alloc_mem(size * sizeof(M));
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
  for (size_t i = 0; i < mm.size(); ++i) {
    sum += qnorm(mm[i]);
  }
  return sum;
}

template <class M>
struct box {
  //
  // like a one element vector
  //
  // Avoid copy constructor when possible
  // (it is likely not be what you think it is)
  //
  bool is_copy;  // do not free memory if is_copy=true
  Handle<M> v;
  //
  box()
  {
    qassert(v.p == NULL);
    is_copy = false;
  }
  box(const box<M>& vp)
  {
    is_copy = true;
    this->v = vp.v;
  }
  box(const M& x)
  {
    qassert(v.p == NULL);
    is_copy = false;
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
      free_mem(v.p, sizeof(M));
    }
    v.init();
    qassert(v.p == NULL);
  }
  //
  qacc void swap(box<M>& x)
  {
    qassert(not is_copy);
    qassert(not x.is_copy);
    Handle<M> t = v;
    v = x.v;
    x.v = t;
  }
  //
  void set(const M& x)
  {
    qassert(not is_copy);
    if (v.p == NULL) {
      v.p = (M*)alloc_mem(sizeof(M));
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
  qacc const M& operator()() const
  {
    return v();
  }
  qacc M& operator()()
  {
    return v();
  }
  //
  qacc bool null() const { return v.null(); }
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
