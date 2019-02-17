#pragma once

#include <qlat/config.h>

#include <cmath>

QLAT_START_NAMESPACE

template <int DIMN, class T = ComplexT>
struct MvectorT {
  T p[DIMN];
  //
  // convert to double array
  double* d() { return (double*)p; }
  const double* d() const { return (const double*)p; }
  //
  // convert to Eigen Matrix
  Eigen::Matrix<T, DIMN, 1>& em()
  {
    return *((Eigen::Matrix<T, DIMN, 1>*)this);
  }
  const Eigen::Matrix<T, DIMN, 1>& em() const
  {
    return *((Eigen::Matrix<T, DIMN, 1>*)this);
  }
  //
  T& operator()(int i)
  {
    qassert(0 <= i && i < DIMN);
    return p[i];
  }
  const T& operator()(int i) const
  {
    qassert(0 <= i && i < DIMN);
    return p[i];
  }
  //
  const MvectorT& operator+=(const MvectorT& x)
  {
    *this = *this + x;
    return *this;
  }
  //
  const MvectorT& operator-=(const MvectorT& x)
  {
    *this = *this - x;
    return *this;
  }
  //
  const MvectorT& operator*=(const T& x)
  {
    *this = *this * x;
    return *this;
  }
  //
  const MvectorT& operator/=(const T& x)
  {
    *this = *this / x;
    return *this;
  }
};

template <int DIMN, class T>
MvectorT<DIMN, T> operator+(const MvectorT<DIMN, T>& x,
                            const MvectorT<DIMN, T>& y)
{
  MvectorT<DIMN, T> ret;
  ret.em() = x.em() + y.em();
  return ret;
}

template <int DIMN, class T>
MvectorT<DIMN, T> operator-(const MvectorT<DIMN, T>& x,
                            const MvectorT<DIMN, T>& y)
{
  MvectorT<DIMN, T> ret;
  ret.em() = x.em() - y.em();
  return ret;
}

template <int DIMN, class T>
MvectorT<DIMN, T> operator-(const MvectorT<DIMN, T>& x)
{
  MvectorT<DIMN, T> ret;
  ret.em() = -x.em();
  return ret;
}

template <int DIMN, class T>
MvectorT<DIMN, T> operator*(const T& x, const MvectorT<DIMN, T>& y)
{
  MvectorT<DIMN, T> ret;
  ret.em() = x * y.em();
  return ret;
}

template <int DIMN, class T>
MvectorT<DIMN, T> operator*(const MvectorT<DIMN, T>& x, const T& y)
{
  MvectorT<DIMN, T> ret;
  ret.em() = x.em() * y;
  return ret;
}

template <int DIMN, class T>
MvectorT<DIMN, T> operator/(const MvectorT<DIMN, T>& x, const T& y)
{
  MvectorT<DIMN, T> ret;
  ret.em() = x.em() / y;
  return ret;
}

template <int DIMN, class T>
Complex dot_product(const MvectorT<DIMN, T>& x, const MvectorT<DIMN, T>& y)
{
  const Complex ret = x.em().adjoint() * y.em();
  return ret;
}

template <int DIMN, class T>
void set_zero(MvectorT<DIMN, T>& m)
{
  memset(&m, 0, sizeof(MvectorT<DIMN, T>));
}

template <int DIMN, class T>
double norm(const MvectorT<DIMN, T>& m)
{
  return m.em().squaredNorm();
}

template <int DIMN, class T>
MvectorT<DIMN, T> vector_conjugate(const MvectorT<DIMN, T>& x)
{
  MvectorT<DIMN, T> ret;
  ret.em() = x.em().conjugate();
  return ret;
}

template <class T = ComplexT>
struct WilsonVectorT : MvectorT<4 * NUM_COLOR, T> {
  WilsonVectorT() {}
  WilsonVectorT(const MvectorT<4 * NUM_COLOR, T>& m) { *this = m; }
  //
  const WilsonVectorT& operator=(const MvectorT<4 * NUM_COLOR, T>& m)
  {
    *this = (const WilsonVectorT&)m;
    return *this;
  }
};

template <class T = ComplexT>
struct SpinVectorT : MvectorT<4, T> {
  SpinVectorT() {}
  SpinVectorT(const MvectorT<4, T>& m) { *this = m; }
  //
  const SpinVectorT& operator=(const MvectorT<4, T>& m)
  {
    *this = (const SpinVectorT&)m;
    return *this;
  }
};

#ifndef QLAT_NO_DEFAULT_TYPE

typedef WilsonVectorT<> WilsonVector;

typedef SpinVectorT<> SpinVector;

#endif

QLAT_END_NAMESPACE

namespace qshow
{
template <int DIMN, class T>
std::string show(const qlat::MvectorT<DIMN, T>& m)
{
  std::ostringstream out;
  out << m.em();
  return out.str();
}

}  // namespace qshow

#ifndef USE_NAMESPACE
using namespace qshow;
#endif
