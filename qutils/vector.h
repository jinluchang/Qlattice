#pragma once

#include <qutils/qutils-vec.h>
#include <qutils/qutils.h>

#include <cstdlib>
#include <cstring>

namespace qlat
{  //

inline size_t& get_alignment()
// qlat parameter
{
  static size_t alignment = 16 * 1024;
  return alignment;
}

inline void* alloc_mem(const long min_size)
{
  const size_t alignment = get_alignment();
  const long n_elem = 1 + (min_size - 1) / alignment;
  const size_t size = n_elem * alignment;
  return aligned_alloc(alignment, size);
}

inline void free_mem(void* ptr) { free(ptr); }

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
      free_mem(v.p);
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
  const M& operator[](long i) const
  {
    return v[i];
  }
  M& operator[](long i)
  {
    return v[i];
  }
  //
  long size() const { return v.n; }
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
