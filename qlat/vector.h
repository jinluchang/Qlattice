#pragma once

#include <eigen3/Eigen/Eigen>

#include <cmath>

QLAT_START_NAMESPACE

template <int DIM>
struct Vector 
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
  // convert to Eigen Vector 
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
  const Vector& operator+=(const Vector& x)
  {
    *this = *this + x;
    return *this;
  }
  //
  const Vector& operator-=(const Vector& x)
  {
    *this = *this - x;
    return *this;
  }
  //
  const Vector& operator*=(const Vector& x)
  {
    *this = *this * x;
    return *this;
  }
  //
  const Vector& operator*=(const Complex& x)
  {
    *this = *this * x;
    return *this;
  }
  //
  const Vector& operator/=(const Complex& x)
  {
    *this = *this / x;
    return *this;
  }
};

template <int DIM>
Vector<DIM> operator+(const Vector<DIM>& x, const Vector<DIM>& y)
{
  Vector<DIM> ret;
  ret.em() = x.em() + y.em();
  return ret;
}

template <int DIM>
Vector<DIM> operator-(const Vector<DIM>& x, const Vector<DIM>& y)
{
  Vector<DIM> ret;
  ret.em() = x.em() - y.em();
  return ret;
}

template <int DIM>
Vector<DIM> operator*(const Vector<DIM>& x, const Vector<DIM>& y)
{
  Vector<DIM> ret;
  ret.em() = x.em() * y.em();
  return ret;
}

template <int DIM>
Vector<DIM> operator*(const Complex& x, const Vector<DIM>& y)
{
  Vector<DIM> ret;
  ret.em() = x * y.em();
  return ret;
}

template <int DIM>
Vector<DIM> operator*(const Vector<DIM>& x, const Complex& y)
{
  Vector<DIM> ret;
  ret.em() = x.em() * y;
  return ret;
}

template <int DIM>
Vector<DIM> operator/(const Vector<DIM>& x, const Complex& y)
{
  Vector<DIM> ret;
  ret.em() = x.em() / y;
  return ret;
}

template <int DIM>
void set_zero(Vector<DIM>& m)
{
  memset(&m, 0, sizeof(Vector<DIM>));
}

template <int DIM>
double norm(const Vector<DIM>& m)
{
  return m.em().squaredNorm();
}

template <int DIM>
Vector<DIM> vector_conjugate(const Vector<DIM>& x)
{
  Vector<DIM> ret;
  ret.em() = x.em().conjugate();
  return ret;
}

struct WilsonVector: Vector<4*NUM_COLOR>
{
  WilsonVector()
  {
  }
  WilsonVector(const Vector<4*NUM_COLOR>& m)
  {
    *this = m;
  }
  //
  const WilsonVector& operator=(const Vector<NUM_COLOR>& m)
  {
    *this = (const WilsonVector&)m;
    return *this;
  }
};

QLAT_END_NAMESPACE

namespace qshow {

template <int DIM>
std::string show(const qlat::Vector<DIM>& m)
{
  std::ostringstream out;
  out << m.em();
  return out.str();
}

}

#ifndef USE_NAMESPACE
using namespace qshow;
#endif
