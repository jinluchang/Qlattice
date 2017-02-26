#pragma once

#include <qlat/config.h>
#include <eigen3/Eigen/Eigen>

#include <cmath>

QLAT_START_NAMESPACE

template <int DIM>
struct Mvector 
{
  static const int dim = DIM;
  Complex p[DIM];
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
  Eigen::Matrix<Complex,DIM,1>& em()
  {
    return *((Eigen::Matrix<Complex,DIM,1>*)this);
  }
  const Eigen::Matrix<Complex,DIM,1>& em() const
  {
    return *((Eigen::Matrix<Complex,DIM,1>*)this);
  }
  //
  Complex& operator()(int i)
  {
    qassert(0 <= i && i < DIM);
    return p[i];
  }
  const Complex& operator()(int i, int j) const
  {
    qassert(0 <= i && i < DIM);
    qassert(0 <= j && j < DIM);
    return p[i * DIM + j];
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
  const Mvector& operator*=(const Mvector& x)
  {
    *this = *this * x;
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

template <int DIM>
Mvector<DIM> operator+(const Mvector<DIM>& x, const Mvector<DIM>& y)
{
  Mvector<DIM> ret;
  ret.em() = x.em() + y.em();
  return ret;
}

template <int DIM>
Mvector<DIM> operator-(const Mvector<DIM>& x, const Mvector<DIM>& y)
{
  Mvector<DIM> ret;
  ret.em() = x.em() - y.em();
  return ret;
}

template <int DIM>
Mvector<DIM> operator*(const Mvector<DIM>& x, const Mvector<DIM>& y)
{
  Mvector<DIM> ret;
  ret.em() = x.em() * y.em();
  return ret;
}

template <int DIM>
Mvector<DIM> operator*(const Complex& x, const Mvector<DIM>& y)
{
  Mvector<DIM> ret;
  ret.em() = x * y.em();
  return ret;
}

template <int DIM>
Mvector<DIM> operator*(const Mvector<DIM>& x, const Complex& y)
{
  Mvector<DIM> ret;
  ret.em() = x.em() * y;
  return ret;
}

template <int DIM>
Mvector<DIM> operator/(const Mvector<DIM>& x, const Complex& y)
{
  Mvector<DIM> ret;
  ret.em() = x.em() / y;
  return ret;
}

template <int DIM>
void set_zero(Mvector<DIM>& m)
{
  memset(&m, 0, sizeof(Mvector<DIM>));
}

template <int DIM>
double norm(const Mvector<DIM>& m)
{
  return m.em().squaredNorm();
}

template <int DIM>
Mvector<DIM> vector_conjugate(const Mvector<DIM>& x)
{
  Mvector<DIM> ret;
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

QLAT_END_NAMESPACE

namespace qshow {

template <int DIM>
std::string show(const qlat::Mvector<DIM>& m)
{
  std::ostringstream out;
  out << m.em();
  return out.str();
}

}

#ifndef USE_NAMESPACE
using namespace qshow;
#endif
