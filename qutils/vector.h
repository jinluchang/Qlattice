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
  static size_t alignment = 16 * 1024;
  return alignment;
}

inline size_t get_aligned_mem_size(const size_t alignment, const long min_size)
{
  const long n_elem = 1 + (min_size - 1) / alignment;
  const size_t size = n_elem * alignment;
  return size;
}

inline size_t& get_mem_cache_max_size()
// qlat parameter
{
  static size_t max_size = 512L * 1024L * 1024L;
  return max_size;
}

struct MemCache {
  size_t mem_cache_size;
  std::unordered_multimap<size_t, void*> db;
  //
  MemCache() {
    mem_cache_size = 0;
  }
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
      free(iter->second);
    }
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
    void* ptr = aligned_alloc(alignment, size);
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
  Vector<M> v;
  //
  vector() { qassert(v.p == NULL); }
  vector(const vector<M>& vp)
  {
    qassert(v.p == NULL);
    *this = vp;
  }
  vector(const long size)
  {
    qassert(v.p == NULL);
    resize(size);
  }
  vector(const long size, const M& x)
  {
    qassert(v.p == NULL);
    resize(size, x);
  }
  //
  ~vector() { clear(); }
  //
  void clear()
  {
    if (v.p != NULL) {
      free_mem(v.p, v.n * sizeof(M));
    }
    v = Vector<M>();
    qassert(v.p == NULL);
  }
  //
  void swap(vector<M>& x)
  {
    Vector<M> t = v;
    v = x.v;
    x.v = t;
  }
  //
  void resize(const long size)
  {
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
        std::memcpy(v.p, vp.v.p, size * sizeof(M));
      } else {
        std::memcpy(v.p, vp.v.p, vp.v.n * sizeof(M));
      }
    }
  }
  void resize(const long size, const M& x)
  {
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
    clear();
    resize(vp.size());
    std::memcpy(v.p, vp.v.p, v.n * sizeof(M));
    return *this;
  }
  //
  const M& operator[](long i) const { return v[i]; }
  M& operator[](long i) { return v[i]; }
  //
  long size() const { return v.size(); }
  //
  M* data() { return v.data(); }
  const M* data() const { return v.data(); }
};

template <class M>
void clear(vector<M>& v)
{
  v.clear();
}

template <class M>
void qswap(vector<M>& v1, vector<M>& v2)
{
  v1.swap(v2);
}

template <class M>
Vector<M> get_data(const vector<M>& v)
{
  return v.v;
}

template <class M>
void set_zero(vector<M>& v)
{
  set_zero(v.v);
}

}  // namespace qlat
