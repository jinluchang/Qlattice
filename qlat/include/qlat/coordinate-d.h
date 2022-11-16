#pragma once

#include <cmath>

namespace qlat
{  //

qacc CoordinateD operator+(const CoordinateD& c1, const CoordinateD& c2)
{
  qassert(false == qisnan(c1));
  qassert(false == qisnan(c2));
  return CoordinateD(c1[0] + c2[0], c1[1] + c2[1], c1[2] + c2[2],
                     c1[3] + c2[3]);
}

qacc CoordinateD operator-(const CoordinateD& c1, const CoordinateD& c2)
{
  qassert(false == qisnan(c1));
  qassert(false == qisnan(c2));
  return CoordinateD(c1[0] - c2[0], c1[1] - c2[1], c1[2] - c2[2],
                     c1[3] - c2[3]);
}

qacc CoordinateD operator-(const CoordinateD& c)
{
  qassert(false == qisnan(c));
  return CoordinateD(-c[0], -c[1], -c[2], -c[3]);
}

qacc CoordinateD operator*(const double a, const CoordinateD& c)
{
  qassert(false == qisnan(c));
  qassert(false == std::isnan(a));
  return CoordinateD(c[0] * a, c[1] * a, c[2] * a, c[3] * a);
}

qacc CoordinateD operator*(const CoordinateD& c, const double a)
{
  return a * c;
}

qacc CoordinateD operator/(const double a, const CoordinateD& c)
{
  qassert(false == qisnan(c));
  qassert(false == std::isnan(a));
  return CoordinateD(a / c[0], a / c[1], a / c[2], a / c[3]);
}

qacc CoordinateD operator/(const CoordinateD& c, const double a)
{
  qassert(false == qisnan(c));
  qassert(false == std::isnan(a));
  return CoordinateD(c[0] / a, c[1] / a, c[2] / a, c[3] / a);
}

qacc CoordinateD operator*(const CoordinateD& c1, const CoordinateD& c2)
{
  qassert(false == qisnan(c1));
  qassert(false == qisnan(c2));
  return CoordinateD(c1[0] * c2[0], c1[1] * c2[1], c1[2] * c2[2],
                     c1[3] * c2[3]);
}

qacc CoordinateD operator/(const CoordinateD& c1, const CoordinateD& c2)
{
  qassert(false == qisnan(c1));
  qassert(false == qisnan(c2));
  return CoordinateD(c1[0] / c2[0], c1[1] / c2[1], c1[2] / c2[2],
                     c1[3] / c2[3]);
}

inline double coordinate_len(const CoordinateD& c)
{
  const double ans = std::sqrt(sqr(c[0]) + sqr(c[1]) + sqr(c[2]) + sqr(c[3]));
  double cmax = 0.0;
  for (int i = 0; i < DIMN; ++i) {
    if (std::abs(c[i]) > cmax) {
      cmax = std::abs(c[i]);
    }
  }
  qassert(false == std::isnan(ans));
  if (0.0 == ans) {
    if (0.0 == cmax) {
      return 0.0;
    } else {
      const double ans = cmax * coordinate_len(c / cmax);
      return std::max(ans, cmax);
    }
  } else {
    return std::max(ans, cmax);
  }
}

qacc double dot_product(const CoordinateD& c1, const CoordinateD& c2)
{
  qassert(false == qisnan(c1));
  qassert(false == qisnan(c2));
  return c1[0] * c2[0] + c1[1] * c2[1] + c1[2] * c2[2] + c1[3] * c2[3];
}

}  // namespace qlat
