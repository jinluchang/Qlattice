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

qacc CoordinateD mod(const CoordinateD& x, const CoordinateD& size)
{
  CoordinateD ret;
  ret[0] = mod(x[0], size[0]);
  ret[1] = mod(x[1], size[1]);
  ret[2] = mod(x[2], size[2]);
  ret[3] = mod(x[3], size[3]);
  return ret;
}

qacc CoordinateD smod(const CoordinateD& x, const CoordinateD& size)
{
  CoordinateD ret;
  ret[0] = smod(x[0], size[0]);
  ret[1] = smod(x[1], size[1]);
  ret[2] = smod(x[2], size[2]);
  ret[3] = smod(x[3], size[3]);
  return ret;
}

qacc CoordinateD smod_sym(const CoordinateD& x, const CoordinateD& size)
{
  CoordinateD ret;
  ret[0] = smod_sym(x[0], size[0]);
  ret[1] = smod_sym(x[1], size[1]);
  ret[2] = smod_sym(x[2], size[2]);
  ret[3] = smod_sym(x[3], size[3]);
  return ret;
}

qacc CoordinateD middle_mod(const CoordinateD& x, const CoordinateD& y,
                            const CoordinateD& size)
{
  CoordinateD ret;
  ret[0] = middle_mod(x[0], y[0], size[0]);
  ret[1] = middle_mod(x[1], y[1], size[1]);
  ret[2] = middle_mod(x[2], y[2], size[2]);
  ret[3] = middle_mod(x[3], y[3], size[3]);
  return ret;
}

}  // namespace qlat
