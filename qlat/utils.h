// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <qlat/config.h>

#include <endian.h>

#include <vector>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>
#include <fstream>

QLAT_START_NAMESPACE

const double PI = 3.141592653589793;

const Complex ii(0, 1);

template <class T>
T sqr(const T& x)
{
  return x * x;
}

inline void set_zero(double& x)
{
  x = 0;
}

inline void set_zero(Complex& x)
{
  x = 0;
}

template <class M, unsigned long N>
void set_zero(std::array<M,N>& arr)
{
  long size = N * sizeof(M);
  std::memset(arr.data(), 0, size);
}

inline void set_unit(double& x, const double& coef = 1.0)
{
  x = coef;
}

inline void set_unit(Complex& x, const Complex& coef = 1.0)
{
  x = coef;
}

template <class M>
void set_zero(std::vector<M>& vec)
{
  long size = vec.size() * sizeof(M);
  std::memset(vec.data(), 0, size);
}

template <class M>
void clear(std::vector<M>& vec)
{
  std::vector<M> empty;
  swap(empty, vec);
}

inline double norm(const double& x)
{
  return x*x;
}

inline double norm(const Complex& x)
{
  return std::norm(x);
}

template <class T, size_t N>
inline double norm(const std::array<T,N>& mm)
{
  double sum = 0.0;
  for (size_t i = 0; i < N; ++i) {
    sum += norm(mm[i]);
  }
  return sum;
}

template <class M, unsigned long N>
bool operator==(const std::array<M,N>& x, const std::array<M,N>& y)
{
  return 0 == memcmp(x.data(), y.data(), N * sizeof(M));
}

template <class M>
bool operator==(const std::vector<M>& x, const std::vector<M>& y)
{
  return x.size() == y.size() && 0 == memcmp(x.data(), y.data(), x.size() * sizeof(M));
}

template <class M>
struct Handle
{
  M* p;
  //
  Handle<M>()
  {
    init();
  }
  Handle<M>(M& obj)
  {
    init(obj);
  }
  //
  void init()
  {
    p = NULL;
  }
  void init(M& obj)
  {
    p = (M*)&obj;
  }
  //
  bool null() const
  {
    return p == NULL;
  }
  //
  M& operator()() const
  {
    qassert(NULL != p);
    return *p;
  }
};

template <class M>
struct ConstHandle
{
  const M* p;
  //
  ConstHandle<M>()
  {
    init();
  }
  ConstHandle<M>(const M& obj)
  {
    init(obj);
  }
  ConstHandle<M>(const Handle<M>& h)
  {
    init(h());
  }
  //
  void init()
  {
    p = NULL;
  }
  void init(const M& obj)
  {
    p = (M*)&obj;
  }
  //
  bool null() const
  {
    return p == NULL;
  }
  //
  const M& operator()() const
  {
    qassert(NULL != p);
    return *p;
  }
};

template <class M> struct Vector;

template <class M, int N>
struct Array
{
  M* p;
  //
  Array<M,N>()
  {
    p = NULL;
  }
  Array<M,N>(const Array<M,N>& arr)
  {
    p = arr.p;
  }
  Array<M,N>(const Vector<M>& vec)
  {
    qassert(N == vec.size());
    p = vec.p;
  }
  Array<M,N>(const std::array<M,N>& arr)
  {
    p = (M*)arr.data();
  }
  Array<M,N>(const M* p_)
  {
    p = (M*)p_;
  }
  Array<M,N>(const M& x)
  {
    qassert(N == 1);
    p = (M*)&x;
  }
  //
  const M& operator[](int i) const
  {
    qassert(0 <= i && i < N);
    return p[i];
  }
  M& operator[](int i)
  {
    qassert(0 <= i && i < N);
    return p[i];
  }
  //
  M* data()
  {
    return p;
  }
  const M* data() const
  {
    return p;
  }
  //
  int size() const
  {
    return N;
  }
  //
  long data_size() const
  {
    return N * sizeof(M);
  }
  //
  const Array<M,N>& operator=(const Array<M,N>& v)
  {
    p = v.p;
    return *this;
  }
  const Array<M,N>& operator=(const Vector<M>& v)
  {
    qassert(N == v.size());
    p = v.p;
    return *this;
  }
};

template <class M, int N>
void set_zero(Array<M,N> arr)
{
  long size = N * sizeof(M);
  std::memset(arr.data(), 0, size);
}

template <class M>
struct Vector
{
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
  Vector<M>(const Array<M,N>& arr)
  {
    p = arr.p;
    n = N;
  }
  template <int N>
  Vector<M>(const std::array<M,N>& arr)
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
  M* data()
  {
    return p;
  }
  const M* data() const
  {
    return p;
  }
  //
  long size() const
  {
    return n;
  }
  //
  long data_size() const
  {
    return n * sizeof(M);
  }
  //
  const Vector<M>& operator=(const Vector<M>& v)
  {
    n = v.n;
    p = v.p;
    return *this;
  }
  template <int N>
  const Vector<M>& operator=(const Array<M,N>& v)
  {
    n = N;
    p = v.p;
    return *this;
  }
};

template <class M>
void set_zero(Vector<M> vec)
{
  std::memset(vec.data(), 0, vec.data_size());
}

template <class M, int N>
Vector<M> get_data(Array<M,N> arr)
{
  return Vector<M>(arr);
}

template <class M>
Vector<M> get_data(Vector<M> vec)
{
  return vec;
}

template <class M, unsigned long N>
Vector<M> get_data(const std::array<M, N>& vec)
{
  return Vector<M>((M*)vec.data(), vec.size());
}

template <class M>
Vector<M> get_data(const std::vector<M>& vec)
{
  return Vector<M>((M*)vec.data(), vec.size());
}

inline Vector<char> get_data(const std::string& str)
{
  return Vector<char>(&str[0], str.length());
}

template <class M>
Vector<M> get_data(const Handle<M>& h)
{
  return Vector<M>(h.p, 1);
}

template <class M>
Vector<M> get_data(const ConstHandle<M>& h)
{
  return Vector<M>(h.p, 1);
}

inline Vector<long> get_data(const long& x)
{
  return Vector<long>(&x, 1);
}

inline Vector<double> get_data(const double& x)
{
  return Vector<double>(&x, 1);
}

inline Vector<int> get_data(const int& x)
{
  return Vector<int>(&x, 1);
}

inline Vector<float> get_data(const float& x)
{
  return Vector<float>(&x, 1);
}

template <class M>
Vector<double> get_data_double(const M& v)
{
  return Vector<double>(&v, sizeof(M) / sizeof(double));
}

template <class M>
Vector<long> get_data_long(const M& v)
{
  return Vector<long>(&v, sizeof(M) / sizeof(long));
}

template <class M>
long get_data_size(const M& x)
{
  return get_data(x).data_size();
}

template <class M, int N>
void assign(std::array<M,N>& vec, const Array<M,N>& src)
{
  memcpy(vec.data(), src.data(), src.data_size());
}

template <class M, int N>
void assign(std::array<M,N>& vec, const Vector<M>& src)
{
  qassert(N == src.size());
  memcpy(vec.data(), src.data(), src.data_size());
}

template <class M, int N>
void assign(std::vector<M>& vec, const Array<M,N>& src)
{
  vec.resize(src.size());
  memcpy(vec.data(), src.data(), src.data_size());
}

template <class M>
void assign(std::vector<M>& vec, const Vector<M>& src)
{
  vec.resize(src.size());
  memcpy(vec.data(), src.data(), src.data_size());
}

template <class M>
void assign(Vector<M> vec, const Vector<M>& src)
{
  qassert(vec.size() == src.size());
  memcpy(vec.data(), src.data(), src.data_size());
}

template <class M, int N>
void assign(Vector<M> vec, const Array<M,N>& src)
{
  qassert(vec.size() == N);
  memcpy(vec.data(), src.data(), src.data_size());
}

template <class M, int N>
void assign(Array<M,N> vec, const Array<M,N>& src)
{
  memcpy(vec.data(), src.data(), src.data_size());
}

template <class M, int N>
void assign(Array<M,N> vec, const Vector<M>& src)
{
  qassert(src.size() == N);
  memcpy(vec.data(), src.data(), src.data_size());
}

template <class M, class N>
void assign_truncate(M& x, const N& y)
{
  if (sizeof(M) <= sizeof(N)) {
    memcpy(&x, &y, sizeof(M));
  } else {
    // if M has a larger size, than leave the extra space untouched
    memcpy(&x, &y, sizeof(N));
  }
}

inline bool is_integer(const double& x)
{
  const double diff = x - (long)x;
  return 1e-6 > diff || diff > 1-1e-6;
}

template <class M>
inline bool is_integer(const std::vector<M>& v)
{
  for (int i = 0; i < (int)v.size(); ++i) {
    if (!is_integer(v[i])) {
      return false;
    }
  }
  return true;
}

template <class M, unsigned long N>
inline bool is_integer(const std::array<M,N>& v)
{
  for (int i = 0; i < N; ++i) {
    if (!is_integer(v[i])) {
      return false;
    }
  }
  return true;
}

inline uint16_t flip_endian_16(uint16_t x)
{
  return
    ((x >> 8)) |
    ((x << 8));
}

inline uint32_t flip_endian_32(uint32_t x)
{
  return
    ((x >> 24)) |
    ((x >>  8) & 0x0000FF00) |
    ((x <<  8) & 0x00FF0000) |
    ((x << 24));
}

inline uint64_t flip_endian_64(uint64_t x)
{
  return
    ((x >> 56)) |
    ((x >> 40) & 0xFF00) |
    ((x >> 24) & 0xFF0000) |
    ((x >>  8) & 0xFF000000) |
    ((x <<  8) & 0xFF00000000) |
    ((x << 24) & 0xFF0000000000) |
    ((x << 40) & 0xFF000000000000) |
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
#if defined(__BYTE_ORDER) && (__BYTE_ORDER != 0) && (__BYTE_ORDER == __BIG_ENDIAN)
  return true;
#else
  return false;
#endif
}

inline bool is_little_endian()
{
  return not is_big_endian();
}

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

QLAT_END_NAMESPACE

namespace qshow {

inline std::string show(const qlat::Complex& x) {
  return ssprintf("(%24.17E + %24.17E j)", x.real(), x.imag());
}

}

#ifndef USE_NAMESPACE
using namespace qshow;
#endif
