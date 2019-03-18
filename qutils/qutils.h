#pragma once

#ifndef OLD_CPP
#include <array>
#endif

#include <vector>
#include <complex>
#include <cassert>

#include "show.h"

// #define SKIP_ASSERT

#ifdef SKIP_ASSERT
#define qassert(x) assert(true)
#else
#define qassert(x)                        \
  {                                       \
    if (not(x)) {                         \
      displayln("qassert failed: " #x);   \
      usleep((useconds_t)(10.0 * 1.0e6)); \
      assert(false);                      \
    }                                     \
  }
#endif

namespace qlat
{  //

typedef std::complex<double> Complex;

typedef std::complex<float> ComplexF;

const double PI = 3.141592653589793;

const Complex ii(0, 1);

template <class T>
T sqr(const T& x)
{
  return x * x;
}

template <class M>
void clear(std::vector<M>& vec)
{
  std::vector<M> empty;
  swap(empty, vec);
}

template <class M>
std::array<M, 0> make_array()
{
  std::array<M, 0> arr;
  return arr;
}

template <class M>
std::array<M, 1> make_array(const M& x)
{
  std::array<M, 1> arr;
  arr[0] = x;
  return arr;
}

template <class M>
std::array<M, 2> make_array(const M& x, const M& x1)
{
  std::array<M, 2> arr;
  arr[0] = x;
  arr[1] = x1;
  return arr;
}

template <class M>
std::array<M, 3> make_array(const M& x, const M& x1, const M& x2)
{
  std::array<M, 3> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  return arr;
}

template <class M>
std::array<M, 4> make_array(const M& x, const M& x1, const M& x2, const M& x3)
{
  std::array<M, 4> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  return arr;
}

template <class M>
std::array<M, 5> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                            const M& x4)
{
  std::array<M, 5> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  arr[4] = x4;
  return arr;
}

template <class M>
std::array<M, 6> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                            const M& x4, const M& x5)
{
  std::array<M, 6> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  arr[4] = x4;
  arr[5] = x5;
  return arr;
}

template <class M>
std::array<M, 7> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                            const M& x4, const M& x5, const M& x6)
{
  std::array<M, 7> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  arr[4] = x4;
  arr[5] = x5;
  arr[6] = x6;
  return arr;
}

template <class M>
std::array<M, 8> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                            const M& x4, const M& x5, const M& x6, const M& x7)
{
  std::array<M, 8> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  arr[4] = x4;
  arr[5] = x5;
  arr[6] = x6;
  arr[7] = x7;
  return arr;
}

template <class M>
std::array<M, 9> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                            const M& x4, const M& x5, const M& x6, const M& x7,
                            const M& x8)
{
  std::array<M, 9> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  arr[4] = x4;
  arr[5] = x5;
  arr[6] = x6;
  arr[7] = x7;
  arr[8] = x8;
  return arr;
}

template <class M>
std::array<M, 10> make_array(const M& x, const M& x1, const M& x2, const M& x3,
                             const M& x4, const M& x5, const M& x6, const M& x7,
                             const M& x8, const M& x9)
{
  std::array<M, 10> arr;
  arr[0] = x;
  arr[1] = x1;
  arr[2] = x2;
  arr[3] = x3;
  arr[4] = x4;
  arr[5] = x5;
  arr[6] = x6;
  arr[7] = x7;
  arr[8] = x8;
  arr[9] = x9;
  return arr;
}

template <class M>
struct Handle {
  M* p;
  //
  Handle<M>() { init(); }
  Handle<M>(M& obj) { init(obj); }
  //
  void init() { p = NULL; }
  void init(M& obj) { p = (M*)&obj; }
  //
  bool null() const { return p == NULL; }
  //
  M& operator()() const
  {
    qassert(NULL != p);
    return *p;
  }
};

template <class M>
struct ConstHandle {
  const M* p;
  //
  ConstHandle<M>() { init(); }
  ConstHandle<M>(const M& obj) { init(obj); }
  ConstHandle<M>(const Handle<M>& h) { init(h()); }
  //
  void init() { p = NULL; }
  void init(const M& obj) { p = (M*)&obj; }
  //
  bool null() const { return p == NULL; }
  //
  const M& operator()() const
  {
    qassert(NULL != p);
    return *p;
  }
};

template <class M>
struct Vector {
  M* p;
  long n;
  //
  Vector<M>()
  {
    p = NULL;
    n = 0;
  }
  Vector<M>(const Vector<M>& vec)
  {
    p = vec.p;
    n = vec.n;
  }
  template <int N>
  Vector<M>(const std::array<M, N>& arr)
  {
    p = (M*)arr.data();
    n = arr.size();
  }
  Vector<M>(const std::vector<M>& vec)
  {
    p = (M*)vec.data();
    n = vec.size();
  }
  Vector<M>(const M* p_, const long n_)
  {
    p = (M*)p_;
    n = n_;
  }
  Vector<M>(const M& x)
  {
    p = (M*)&x;
    n = 1;
  }
  //
  const M& operator[](long i) const
  {
    qassert(0 <= i && i < n);
    return p[i];
  }
  M& operator[](long i)
  {
    qassert(0 <= i && i < n);
    return p[i];
  }
  //
  M* data() { return p; }
  const M* data() const { return p; }
  //
  long size() const { return n; }
  //
  long data_size() const { return n * sizeof(M); }
  //
  const Vector<M>& operator=(const Vector<M>& v)
  {
    n = v.n;
    p = v.p;
    return *this;
  }
};

inline uint16_t flip_endian_16(uint16_t x) { return ((x >> 8)) | ((x << 8)); }

inline uint32_t flip_endian_32(uint32_t x)
{
  return ((x >> 24)) | ((x >> 8) & 0x0000FF00) | ((x << 8) & 0x00FF0000) |
         ((x << 24));
}

inline uint64_t flip_endian_64(uint64_t x)
{
  return ((x >> 56)) | ((x >> 40) & 0xFF00) | ((x >> 24) & 0xFF0000) |
         ((x >> 8) & 0xFF000000) | ((x << 8) & 0xFF00000000) |
         ((x << 24) & 0xFF0000000000) | ((x << 40) & 0xFF000000000000) |
         ((x << 56));
}

inline void flip_endian_16(void* str, const size_t len)
{
  qassert(0 == len % 2);
  uint16_t* p = (uint16_t*)str;
  for (size_t i = 0; i < len / 2; ++i) {
    p[i] = flip_endian_16(p[i]);
  }
}

inline void flip_endian_32(void* str, const size_t len)
{
  qassert(0 == len % 4);
  uint32_t* p = (uint32_t*)str;
  for (size_t i = 0; i < len / 4; ++i) {
    p[i] = flip_endian_32(p[i]);
  }
}

inline void flip_endian_64(void* str, const size_t len)
{
  qassert(0 == len % 8);
  uint64_t* p = (uint64_t*)str;
  for (size_t i = 0; i < len / 8; ++i) {
    p[i] = flip_endian_64(p[i]);
  }
}

inline bool is_big_endian()
{
#if defined(__BYTE_ORDER) && (__BYTE_ORDER != 0) && \
    (__BYTE_ORDER == __BIG_ENDIAN)
  return true;
#else
  return false;
#endif
}

inline bool is_little_endian() { return not is_big_endian(); }

inline void to_from_little_endian_16(void* str, const size_t len)
{
  qassert(0 == len % 2);
  if (is_big_endian()) {
    flip_endian_16(str, len);
  }
}

inline void to_from_little_endian_32(void* str, const size_t len)
{
  qassert(0 == len % 4);
  if (is_big_endian()) {
    flip_endian_32(str, len);
  }
}

inline void to_from_little_endian_64(void* str, const size_t len)
{
  qassert(0 == len % 8);
  if (is_big_endian()) {
    flip_endian_64(str, len);
  }
}

inline void to_from_big_endian_16(void* str, const size_t len)
{
  qassert(0 == len % 2);
  if (is_little_endian()) {
    flip_endian_16(str, len);
  }
}

inline void to_from_big_endian_32(void* str, const size_t len)
{
  qassert(0 == len % 4);
  if (is_little_endian()) {
    flip_endian_32(str, len);
  }
}

inline void to_from_big_endian_64(void* str, const size_t len)
{
  qassert(0 == len % 8);
  if (is_little_endian()) {
    flip_endian_64(str, len);
  }
}

template <class M>
void to_from_little_endian_16(Vector<M> v)
{
  to_from_little_endian_16((void*)v.data(), v.data_size());
}

template <class M>
void to_from_little_endian_32(Vector<M> v)
{
  to_from_little_endian_32((void*)v.data(), v.data_size());
}

template <class M>
void to_from_little_endian_64(Vector<M> v)
{
  to_from_little_endian_64((void*)v.data(), v.data_size());
}

template <class M>
void to_from_big_endian_16(Vector<M> v)
{
  to_from_big_endian_16((void*)v.data(), v.data_size());
}

template <class M>
void to_from_big_endian_32(Vector<M> v)
{
  to_from_big_endian_32((void*)v.data(), v.data_size());
}

template <class M>
void to_from_big_endian_64(Vector<M> v)
{
  to_from_big_endian_64((void*)v.data(), v.data_size());
}

inline void from_big_endian_32(char* str, const size_t len)
// obsolete
{
  to_from_big_endian_32(str, len);
}

inline void from_big_endian_64(char* str, const size_t len)
// obsolete
{
  to_from_big_endian_64(str, len);
}

}  // namespace qlat

#ifndef USE_NAMESPACE
using namespace qlat;
#endif
