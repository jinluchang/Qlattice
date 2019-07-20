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
#define qassert(x)                            \
  {                                           \
    if (not(x)) {                             \
      qlat::displayln("qassert failed: " #x); \
      usleep((useconds_t)(10.0 * 1.0e6));     \
      assert(false);                          \
    }                                         \
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

inline void set_zero(double& x) { x = 0; }

inline void set_zero(Complex& x) { x = 0; }

template <class M, unsigned long N>
void set_zero(std::array<M, N>& arr)
{
  long size = N * sizeof(M);
  std::memset(arr.data(), 0, size);
}

inline void set_unit(double& x, const double& coef = 1.0) { x = coef; }

inline void set_unit(Complex& x, const Complex& coef = 1.0) { x = coef; }

template <class M>
void set_zero(std::vector<M>& vec)
{
  long size = vec.size() * sizeof(M);
  std::memset((void*)vec.data(), 0, size);
}

inline double qnorm(const double& x) { return x * x; }

inline double qnorm(const Complex& x) { return std::norm(x); }

template <class T, size_t N>
double qnorm(const std::array<T, N>& mm)
{
  double sum = 0.0;
  for (size_t i = 0; i < N; ++i) {
    sum += qnorm(mm[i]);
  }
  return sum;
}

template <class Vec>
bool is_equal_vec(const Vec& v1, const Vec& v2)
{
  const bool b = v1.size() == v2.size();
  if (not b) {
    return false;
  } else {
    const long s = v1.size();
    for (long i = 0; i < s; ++i) {
      if (not(v1[i] == v2[i])) {
        return false;
      }
    }
    return true;
  }
}

template <class M>
bool operator==(const std::vector<M>& v1, const std::vector<M>& v2)
{
  return is_equal_vec(v1, v2);
}

template <class M, int N>
bool operator==(const std::array<M, N>& v1, const std::array<M, N>& v2)
{
  return is_equal_vec(v1, v2);
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

template <typename F>
double simpson(const F& f, const double a, const double b)
{
  return 1.0 / 6.0 * (f(a) + 4 * f(0.5 * (a + b)) + f(b)) * (b - a);
}

inline int adaptive_simpson_min_level()
{
  static int level = 6;
  return level;
}

inline int adaptive_simpson_max_level()
{
  static int level = 20;
  return level;
}

template <typename F>
double adaptive_simpson_level(const F& f, const double a, const double b,
                              const double eps, const int level)
{
  const double w = simpson(f, a, b);
  const double l = simpson(f, a, 0.5 * (a + b));
  const double r = simpson(f, 0.5 * (a + b), b);
  const double error = 1.0 / 15.0 * (l + r - w);
  const double result = l + r + error;
  if ((level >= adaptive_simpson_min_level() and std::abs(error) <= eps) or
      level >= adaptive_simpson_max_level()) {
    return result;
  } else {
    return adaptive_simpson_level(f, a, 0.5 * (a + b), eps / 2.0, level + 1) +
           adaptive_simpson_level(f, 0.5 * (a + b), b, eps / 2.0, level + 1);
  }
  return result;
}

template <typename F>
struct AdaptiveSimpsonToInf
{
  ConstHandle<F> f;
  double start;
  //
  double operator()(const double x) const
  {
    if (x == 1.0) {
      return 0.0;
    } else {
      return f()(start + x / (1.0 - x)) / sqr(1.0 - x);
    }
  }
};

template <typename F>
struct AdaptiveSimpsonFromInf
{
  ConstHandle<F> f;
  double end;
  //
  double operator()(const double x) const
  {
    if (x == 1.0) {
      return 0.0;
    } else {
      return f()(end - x / (1.0 - x)) / sqr(1.0 - x);
    }
  }
};

template <typename F>
double adaptive_simpson(const F& f, const double a, const double b,
                        const double eps)
{
  if (b < a) {
    return -adaptive_simpson(f, b, a, eps);
  } else if (a == b) {
    return 0.0;
  } else {
    const double inf = 1.0 / 0.0;
    if (a == -inf and b == inf) {
      return adaptive_simpson(f, a, 0.0, eps / 2.0) +
             adaptive_simpson(f, 0.0, b, eps / 2.0);
    } else if (b == inf) {
      AdaptiveSimpsonToInf<F> ff;
      ff.f.init(f);
      ff.start = a;
      return adaptive_simpson_level(ff, 0.0, 1.0, eps, 0);
    } else if (a == -inf) {
      AdaptiveSimpsonFromInf<F> ff;
      ff.f.init(f);
      ff.end = b;
      return adaptive_simpson_level(ff, 0.0, 1.0, eps, 0);
    } else {
      return adaptive_simpson_level(f, a, b, eps, 0);
    }
  }
}

inline void split_work(long& start, long& size, const long total, const long num_worker, const long id_worker)
{
  const long size_max = (total - 1) / num_worker + 1;
  start = std::min(id_worker * size_max, total);
  const long stop = std::min(start + size_max, total);
  size = stop - start;
}

inline long find_worker(const long idx, const long total, const long num_worker)
{
  const long size_max = (total - 1) / num_worker + 1;
  return idx / size_max;
}

}  // namespace qlat

#ifndef USE_NAMESPACE
    using namespace qlat;
#endif
