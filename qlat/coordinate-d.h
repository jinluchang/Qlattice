#pragma once

#include <qlat/utils.h>

#include <cmath>

namespace std {

  template <class M, unsigned long N>
  inline bool isnan(const array<M,N>& arr)
  {
    for (int i = 0; i < (int)N; ++i) {
      if (isnan(arr[i])) {
        return true;
      }
    }
    return false;
  }

  template <class M>
  inline bool isnan(const vector<M>& arr)
  {
    for (size_t i = 0; i < arr.size(); ++i) {
      if (isnan(arr[i])) {
        return true;
      }
    }
    return false;
  }

}

QLAT_START_NAMESPACE

struct CoordinateD : public std::array<double,DIM>
{
  CoordinateD()
  {
    memset(this, 0, sizeof(CoordinateD));
  }
  CoordinateD(const std::array<double,DIM>& arr)
  {
    CoordinateD& c = *this;
    c = arr;
    qassert(false == isnan(c));
  }
  CoordinateD(const double x0, const double x1, const double x2, const double x3)
  {
    qassert(DIM == 4);
    CoordinateD& c = *this;
    c[0] = x0;
    c[1] = x1;
    c[2] = x2;
    c[3] = x3;
    qassert(false == isnan(c));
  }
  CoordinateD(const Coordinate& x)
  {
    CoordinateD& c = *this;
    for (int i = 0; i < DIM; ++i) {
      c[i] = x[i];
    }
  }
};

inline double coordinate_len(const CoordinateD& c)
{
  const double ans = std::sqrt(sqr(c[0]) + sqr(c[1]) + sqr(c[2]) + sqr(c[3]));
  qassert(false == std::isnan(ans));
  return ans;
}

inline CoordinateD operator+(const CoordinateD& c1, const CoordinateD& c2)
{
  qassert(false == isnan(c1));
  qassert(false == isnan(c2));
  return CoordinateD(c1[0] + c2[0], c1[1] + c2[1], c1[2] + c2[2], c1[3] + c2[3]);
}

inline CoordinateD operator-(const CoordinateD& c1, const CoordinateD& c2)
{
  qassert(false == isnan(c1));
  qassert(false == isnan(c2));
  return CoordinateD(c1[0] - c2[0], c1[1] - c2[1], c1[2] - c2[2], c1[3] - c2[3]);
}

inline CoordinateD operator-(const CoordinateD& c)
{
  qassert(false == isnan(c));
  return CoordinateD(-c[0], -c[1], -c[2], -c[3]);
}

inline CoordinateD operator*(const double a, const CoordinateD& c)
{
  qassert(false == isnan(c));
  return CoordinateD(c[0]*a, c[1]*a, c[2]*a, c[3]*a);
}

inline CoordinateD operator*(const CoordinateD& c, const double a)
{
  return a * c;
}

inline CoordinateD operator/(const CoordinateD& c, const double a)
{
  qassert(false == isnan(c));
  return (1.0 / a) * c;
}

inline double dot_product(const CoordinateD& c1, const CoordinateD& c2)
{
  qassert(false == isnan(c1));
  qassert(false == isnan(c2));
  return c1[0] * c2[0] + c1[1] * c2[1] + c1[2] * c2[2] + c1[3] * c2[3];
}

QLAT_END_NAMESPACE

namespace qshow {

inline std::string show(const qlat::CoordinateD& c)
{
  return ssprintf("(%23.16e,%23.16e,%23.16e,%23.16e)", c[0], c[1], c[2], c[3]);
}

}

#ifndef USE_NAMESPACE qshow;
using namespace qshow;
#endif
