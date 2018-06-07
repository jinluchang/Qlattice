#pragma once

#include <qlat/config.h>
#include <eigen3/Eigen/Eigen>

#include <cmath>

QLAT_START_NAMESPACE

template <int DIMN>
struct Mvector 
{
  Complex p[DIMN];
  //
  // convert to double array
  double* d()
  {
    return (double*)p;
  }
  const double* d() const
  {
    return (const double*)p;
  }
  //
  // convert to Eigen Matrix 
  Eigen::Matrix<Complex,DIMN,1>& em()
  {
    return *((Eigen::Matrix<Complex,DIMN,1>*)this);
  }
  const Eigen::Matrix<Complex,DIMN,1>& em() const
  {
    return *((Eigen::Matrix<Complex,DIMN,1>*)this);
  }
  //
  Complex& operator()(int i)
  {
    qassert(0 <= i && i < DIMN);
    return p[i];
  }
  const Complex& operator()(int i) const
  {
    qassert(0 <= i && i < DIMN);
    return p[i];
  }
  //
  const Mvector& operator+=(const Mvector& x)
  {
    *this = *this + x;
    return *this;
  }
  //
  const Mvector& operator-=(const Mvector& x)
  {
    *this = *this - x;
    return *this;
  }
  //
  const Mvector& operator*=(const Complex& x)
  {
    *this = *this * x;
    return *this;
  }
  //
  const Mvector& operator/=(const Complex& x)
  {
    *this = *this / x;
    return *this;
  }
};

template <int DIMN>
Mvector<DIMN> operator+(const Mvector<DIMN>& x, const Mvector<DIMN>& y)
{
  Mvector<DIMN> ret;
  ret.em() = x.em() + y.em();
  return ret;
}

template <int DIMN>
Mvector<DIMN> operator-(const Mvector<DIMN>& x, const Mvector<DIMN>& y)
{
  Mvector<DIMN> ret;
  ret.em() = x.em() - y.em();
  return ret;
}

template <int DIMN>
Mvector<DIMN> operator*(const Complex& x, const Mvector<DIMN>& y)
{
  Mvector<DIMN> ret;
  ret.em() = x * y.em();
  return ret;
}

template <int DIMN>
Mvector<DIMN> operator*(const Mvector<DIMN>& x, const Complex& y)
{
  Mvector<DIMN> ret;
  ret.em() = x.em() * y;
  return ret;
}

template <int DIMN>
Mvector<DIMN> operator/(const Mvector<DIMN>& x, const Complex& y)
{
  Mvector<DIMN> ret;
  ret.em() = x.em() / y;
  return ret;
}

template <int DIMN>
void set_zero(Mvector<DIMN>& m)
{
  memset(&m, 0, sizeof(Mvector<DIMN>));
}

template <int DIMN>
double norm(const Mvector<DIMN>& m)
{
  return m.em().squaredNorm();
}

template <int DIMN>
Mvector<DIMN> vector_conjugate(const Mvector<DIMN>& x)
{
  Mvector<DIMN> ret;
  ret.em() = x.em().conjugate();
  return ret;
}

struct WilsonVector: Mvector<4*NUM_COLOR>
{
  WilsonVector()
  {
  }
  WilsonVector(const Mvector<4*NUM_COLOR>& m)
  {
    *this = m;
  }
  //
  const WilsonVector& operator=(const Mvector<NUM_COLOR>& m)
  {
    *this = (const WilsonVector&)m;
    return *this;
  }
};

inline WilsonVector operator*(const Complex& x, const WilsonVector& y)
{
  WilsonVector ret;
  ret.em() = x * y.em();
  return ret;
}

inline WilsonVector operator+(const WilsonVector& x, const WilsonVector& y)
{
  WilsonVector ret;
  ret.em() = x.em() + y.em();
  return ret;
}

QLAT_END_NAMESPACE

namespace qshow {

template <int DIMN>
std::string show(const qlat::Mvector<DIMN>& m)
{
  std::ostringstream out;
  out << m.em();
  return out.str();
}

}

#ifndef USE_NAMESPACE
using namespace qshow;
#endif
