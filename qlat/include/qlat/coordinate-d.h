#pragma once

#include <cmath>

namespace qlat
{  //

qacc bool qisnan(const float& arg) { return std::isnan(arg); }

qacc bool qisnan(const double& arg) { return std::isnan(arg); }

template <class T>
bool qisnan(const std::complex<T>& arg)
{
  return qisnan(arg.real()) or qisnan(arg.imag());
}

template <class M, unsigned long N>
qacc bool qisnan(const array<M, N>& arr)
{
  for (int i = 0; i < (int)N; ++i) {
    if (qisnan(arr[i])) {
      return true;
    }
  }
  return false;
}

template <class M>
bool qisnan(const std::vector<M>& arr)
{
  for (size_t i = 0; i < arr.size(); ++i) {
    if (qisnan(arr[i])) {
      return true;
    }
  }
  return false;
}

struct CoordinateD : public array<double, DIMN> {
  qacc CoordinateD() { memset(this, 0, sizeof(CoordinateD)); }
  qacc CoordinateD(const array<double, DIMN>& arr)
  {
    CoordinateD& c = *this;
    c = arr;
    qassert(false == qisnan(c));
  }
  qacc CoordinateD(const double x0, const double x1, const double x2,
              const double x3)
  {
    qassert(DIMN == 4);
    CoordinateD& c = *this;
    c[0] = x0;
    c[1] = x1;
    c[2] = x2;
    c[3] = x3;
    qassert(false == qisnan(c));
  }
  qacc CoordinateD(const Coordinate& x)
  {
    CoordinateD& c = *this;
    for (int i = 0; i < DIMN; ++i) {
      c[i] = x[i];
    }
  }
};

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

inline std::string show(const CoordinateD& c)
{
  return ssprintf("(%23.16e,%23.16e,%23.16e,%23.16e)", c[0], c[1], c[2], c[3]);
}

}  // namespace qlat
