#pragma once

#include <lqps/config.h>

#include <array>
#include <vector>
#include <cassert>

LQPS_START_NAMESPACE

const double PI = 3.141592653589793;

template <class T>
T sqr(T x) {
  return x * x;
}

template <class M, int N>
struct Array
{
  M* p;
  //
  Array<M,N>()
  {
    p = NULL;
  }
  Array<M,N>(std::array<M,N>& arr)
  {
    p = arr.data();
  }
  Array<M,N>(const M* p_)
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
  template <int N>
  Vector<M>(std::array<M,N>& arr)
  {
    p = arr.data();
    n = arr.szie();
  }
  Vector<M>(std::vector<M>& vec)
  {
    p = vec.data();
    n = vec.size();
  }
  template <int N>
  Vector<M>(Array<M,N>& arr)
  {
    p = arr.data;
    n = N;
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
    return n;
  }
  //
  const Vector<M>& operator=(const Vector<M>& v)
  {
    assert(v.size() == n);
    memcpy(this, v.data(), n * sizeof(M));
    return *this;
  }
  template <int N>
  const Vector<M>& operator=(const Array<M,N>& v)
  {
    assert(v.size() == n);
    memcpy(this, v.data(), n * sizeof(M));
    return *this;
  }
};

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

LQPS_END_NAMESPACE
