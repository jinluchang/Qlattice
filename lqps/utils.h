#pragma once

#include <lqps/config.h>

#include <array>
#include <vector>
#include <iostream>
#include <cassert>
#include <cstdio>
#include <cstdlib>

LQPS_START_NAMESPACE

const double PI = 3.141592653589793;

const Complex ii(0, 1);

template <class T>
T sqr(const T& x)
{
  return x * x;
}

template <class M>
void setZero(std::vector<M>& vec)
{
  long size = vec.size() * sizeof(M);
  memset(vec.data(), 0, size);
}

template <class M>
void swap(std::vector<M>& v1, std::vector<M>& v2)
{
  std::vector<M> tmp(v2);
  v2 = v1;
  v1 = tmp;
}

inline double norm(const double& x)
{
  return x*x;
}

inline double norm(const Complex& x)
{
  return std::norm(x);
}

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
    assert(N == vec.size());
    p = vec.p;
  }
  Array<M,N>(std::array<M,N>& arr)
  {
    p = arr.data();
  }
  Array<M,N>(M* p_)
  {
    p = p_;
  }
  //
  const M& operator[](int i) const
  {
    assert(0 <= i && i < N);
    return p[i];
  }
  M& operator[](int i)
  {
    assert(0 <= i && i < N);
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
  const Array<M,N>& operator=(const Array<M,N>& v)
  {
    memcpy(this, v.data(), N * sizeof(M));
    return *this;
  }
  const Array<M,N>& operator=(const Vector<M>& v)
  {
    assert(N == v.size());
    memcpy(this, v.data(), N * sizeof(M));
    return *this;
  }
};

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
  Vector<M>(std::array<M,N>& arr)
  {
    p = arr.data();
    n = arr.size();
  }
  Vector<M>(std::vector<M>& vec)
  {
    p = vec.data();
    n = vec.size();
  }
  Vector<M>(M* p_, const long n_)
  {
    p = p_;
    n = n_;
  }
  //
  const M& operator[](long i) const
  {
    assert(0 <= i && i < n);
    return p[i];
  }
  M& operator[](long i)
  {
    assert(0 <= i && i < n);
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
  long dataSize() const
  {
    return n * sizeof(M);
  }
  //
  const Vector<M>& operator=(const Vector<M>& v)
  {
    assert(v.size() == n);
    memcpy(this, v.data(), v.size() * sizeof(M));
    return *this;
  }
  template <int N>
  const Vector<M>& operator=(const Array<M,N>& v)
  {
    assert(v.size() == n);
    memcpy(this, v.data(), v.size() * sizeof(M));
    return *this;
  }
};

template <class M>
Vector<M> getData(const std::vector<M>& vec)
{
  return Vector<M>((M*)vec.data(), vec.size());
}

template <class T>
long getDataSize(const T& x)
{
  return getData(x).dataSize();
}

template <class M, int N>
void assign(std::array<M,N>& vec, const Array<M,N>& src)
{
  memcpy(vec.data(), src.data(), src.size() * sizeof(M));
}

template <class M, int N>
void assign(std::array<M,N>& vec, const Vector<M>& src)
{
  assert(N == src.size());
  memcpy(vec.data(), src.data(), src.size() * sizeof(M));
}

template <class M, int N>
void assign(std::vector<M>& vec, const Array<M,N>& src)
{
  vec.resize(src.size());
  memcpy(vec.data(), src.data(), src.size() * sizeof(M));
}

template <class M>
void assign(std::vector<M>& vec, const Vector<M>& src)
{
  vec.resize(src.size());
  memcpy(vec.data(), src.data(), src.size() * sizeof(M));
}

inline int mod(const int x, const int len) {
  assert(0 < len);
  const int m = x % len;
  if (0 <= m) {
    return m;
  } else {
    return m + len;
  }
}

inline int signMod(const int x, const int len) {
  assert(0 < len);
  const int m = mod(x, len);
  if (m * 2 < len) {
    return m;
  } else {
    return m - len;
  }
}

inline int middleMod(const int x, const int y, const int len) {
  assert(0 < len);
  const int xm = mod(x, len);
  const int ym = mod(y, len);
  if (xm <= ym) {
    const int r = signMod(ym - xm, len);
    return mod(xm + r/2, len);
  } else {
    const int r = signMod(xm - ym, len);
    return mod(ym + r/2, len);
  }
}

inline void regularizeCoordinate(Coordinate& x, const Coordinate& size) {
  x[0] = mod(x[0], size[0]);
  x[1] = mod(x[1], size[1]);
  x[2] = mod(x[2], size[2]);
  x[3] = mod(x[3], size[3]);
}

inline long distance2RelativeCoordinateG(const Coordinate& xg) {
  return sqr((long)xg[0]) + sqr((long)xg[1]) + sqr((long)xg[2]) + sqr((long)xg[3]);
}

inline double distanceRelativeCoordinateG(const Coordinate& xg) {
  return sqrt(distance2RelativeCoordinateG(xg));
}

inline void coordinateFromIndex(Coordinate& x, long index, const Coordinate& size) {
  x[0] = index % size[0];
  index /= size[0];
  x[1] = index % size[1];
  index /= size[1];
  x[2] = index % size[2];
  index /= size[2];
  x[3] = index % size[3];
}

inline long indexFromCoordinate(const Coordinate& x, const Coordinate& size) {
  return (((x[3] * size[2]) + x[2]) * size[1] + x[1]) * size[0] + x[0];
}

inline void shiftCoordinateAdd(Coordinate& x, const Coordinate& shift)
{
  x[0] += shift[0];
  x[1] += shift[1];
  x[2] += shift[2];
  x[3] += shift[3];
}

inline void shiftCoordinateSub(Coordinate& x, const Coordinate& shift)
{
  x[0] -= shift[0];
  x[1] -= shift[1];
  x[2] -= shift[2];
  x[3] -= shift[3];
}

inline int vssprintf(std::string& str, const char* fmt, va_list args)
{
  char* cstr;
  int ret = vasprintf(&cstr, fmt, args);
  str += std::string(cstr);
  std::free(cstr);
  return ret;
}

inline int ssprintf(std::string& str, const char* fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  return vssprintf(str, fmt, args);
}

inline std::string vssprintf(const char* fmt, va_list args)
{
  std::string str;
  vssprintf(str, fmt, args);
  return str;
}

inline std::string ssprintf(const char* fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  return vssprintf(fmt, args);
}

inline std::string show() {
  return "";
}

inline std::string show(const long& x) {
  return ssprintf("%ld", x);
}

inline std::string show(const double& x) {
  std::ostringstream out;
  out.precision(16);
  out << std::scientific;
  out << x;
  return out.str();
}

inline std::string show(const bool& x) {
  return x ? "true" : "false";
}

inline std::string show(const std::string& x) {
  std::ostringstream out;
  out << x;
  return out.str();
}

inline std::string show(const Coordinate& x) {
  return ssprintf("%dx%dx%dx%d", x[0], x[1], x[2], x[3]);
}

LQPS_END_NAMESPACE
