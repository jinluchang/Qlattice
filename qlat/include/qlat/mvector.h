#pragma once

#include <qlat/config.h>

#include <cmath>

namespace qlat
{  //

template <int DIMN, class T = ComplexT>
struct alignas(QLAT_ALIGNED_BYTES) MvectorT
{
  T p[DIMN];
  //
  // convert to double array
  qacc double* d() { return (double*)p; }
  qacc const double* d() const { return (const double*)p; }
  //
  // convert to Eigen Matrix
  qacc Eigen::Matrix<T, DIMN, 1>& em()
  {
    return *((Eigen::Matrix<T, DIMN, 1>*)this);
  }
  qacc const Eigen::Matrix<T, DIMN, 1>& em() const
  {
    return *((Eigen::Matrix<T, DIMN, 1>*)this);
  }
  //
  qacc T& operator()(int i)
  {
    qassert(0 <= i && i < DIMN);
    return p[i];
  }
  qacc const T& operator()(int i) const
  {
    qassert(0 <= i && i < DIMN);
    return p[i];
  }
  //
  qacc const MvectorT& operator+=(const MvectorT& x)
  {
    *this = *this + x;
    return *this;
  }
  //
  qacc const MvectorT& operator-=(const MvectorT& x)
  {
    *this = *this - x;
    return *this;
  }
  //
  qacc const MvectorT& operator*=(const T& x)
  {
    *this = *this * x;
    return *this;
  }
  //
  qacc const MvectorT& operator/=(const T& x)
  {
    *this = *this / x;
    return *this;
  }
};

template <int DIMN, class T>
qacc MvectorT<DIMN, T> operator+(const MvectorT<DIMN, T>& x,
                                 const MvectorT<DIMN, T>& y)
{
  MvectorT<DIMN, T> ret;
  ret.em() = x.em() + y.em();
  return ret;
}

template <int DIMN, class T>
qacc MvectorT<DIMN, T> operator-(const MvectorT<DIMN, T>& x,
                                 const MvectorT<DIMN, T>& y)
{
  MvectorT<DIMN, T> ret;
  ret.em() = x.em() - y.em();
  return ret;
}

template <int DIMN, class T>
qacc MvectorT<DIMN, T> operator-(const MvectorT<DIMN, T>& x)
{
  MvectorT<DIMN, T> ret;
  ret.em() = -x.em();
  return ret;
}

template <int DIMN, class T>
qacc MvectorT<DIMN, T> operator*(const T& x, const MvectorT<DIMN, T>& y)
{
  MvectorT<DIMN, T> ret;
  ret.em() = x * y.em();
  return ret;
}

template <int DIMN, class T>
qacc MvectorT<DIMN, T> operator*(const MvectorT<DIMN, T>& x, const T& y)
{
  MvectorT<DIMN, T> ret;
  ret.em() = x.em() * y;
  return ret;
}

template <int DIMN, class T>
qacc MvectorT<DIMN, T> operator/(const MvectorT<DIMN, T>& x, const T& y)
{
  MvectorT<DIMN, T> ret;
  ret.em() = x.em() / y;
  return ret;
}

template <int DIMN, class T>
qacc Complex dot_product(const MvectorT<DIMN, T>& x, const MvectorT<DIMN, T>& y)
{
  const Complex ret = x.em().adjoint() * y.em();
  return ret;
}

template <int DIMN, class T>
qacc void set_zero(MvectorT<DIMN, T>& m)
{
  std::memset((void*)&m, 0, sizeof(MvectorT<DIMN, T>));
}

template <int DIMN, class T>
qacc void set_unit(MvectorT<DIMN, T>& m, const Complex& coef = 1.0)
{
  for (int i = 0; i < DIMN; ++i) {
    m(i) = coef;
  }
}

template <int DIMN, class T>
qacc double qnorm(const MvectorT<DIMN, T>& m)
{
  return m.em().squaredNorm();
}

template <int DIMN, class T>
qacc MvectorT<DIMN, T> vector_conjugate(const MvectorT<DIMN, T>& x)
{
  MvectorT<DIMN, T> ret;
  ret.em() = x.em().conjugate();
  return ret;
}

template <class T = ComplexT>
struct WilsonVectorT : MvectorT<4 * NUM_COLOR, T> {
  qacc WilsonVectorT() {}
  qacc WilsonVectorT(const MvectorT<4 * NUM_COLOR, T>& m) { *this = m; }
  //
  qacc const WilsonVectorT& operator=(const MvectorT<4 * NUM_COLOR, T>& m)
  {
    *this = (const WilsonVectorT&)m;
    return *this;
  }
};

template <class T = ComplexT>
struct SpinVectorT : MvectorT<4, T> {
  qacc SpinVectorT() {}
  qacc SpinVectorT(const MvectorT<4, T>& m) { *this = m; }
  //
  qacc const SpinVectorT& operator=(const MvectorT<4, T>& m)
  {
    *this = (const SpinVectorT&)m;
    return *this;
  }
};

#ifndef QLAT_NO_DEFAULT_TYPE

typedef WilsonVectorT<> WilsonVector;

typedef SpinVectorT<> SpinVector;

#endif

template <int DIMN, class T>
std::string show(const MvectorT<DIMN, T>& m)
{
  std::ostringstream out;
  out << m.em();
  return out.str();
}

}  // namespace qlat
