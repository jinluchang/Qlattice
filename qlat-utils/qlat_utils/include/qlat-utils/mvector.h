#pragma once

#include <qlat-utils/mat-vec.h>

#include <cmath>

namespace qlat
{  //

template <int DIMN, class T>
qacc Vector<T> get_data(const MvectorT<DIMN, T>& m)
{
  return Vector<T>(m.p, DIMN);
}

template <int DIMN, class T>
qacc MvectorT<DIMN, T> operator+(const MvectorT<DIMN, T>& x,
                                 const MvectorT<DIMN, T>& y)
{
  MvectorT<DIMN, T> ret;
  for (int i = 0; i < DIMN; ++i) {
    ret.p[i] = x.p[i] + y.p[i];
  }
  return ret;
}

template <int DIMN, class T>
qacc MvectorT<DIMN, T> operator-(const MvectorT<DIMN, T>& x,
                                 const MvectorT<DIMN, T>& y)
{
  MvectorT<DIMN, T> ret;
  for (int i = 0; i < DIMN; ++i) {
    ret.p[i] = x.p[i] - y.p[i];
  }
  return ret;
}

template <int DIMN, class T>
qacc MvectorT<DIMN, T> operator-(const MvectorT<DIMN, T>& x)
{
  MvectorT<DIMN, T> ret;
  for (int i = 0; i < DIMN; ++i) {
    ret.p[i] = -x.p[i];
  }
  return ret;
}

template <int DIMN, class T>
qacc MvectorT<DIMN, T> operator*(const T& x, const MvectorT<DIMN, T>& y)
{
  MvectorT<DIMN, T> ret;
  for (int i = 0; i < DIMN; ++i) {
    ret.p[i] = x * y.p[i];
  }
  return ret;
}

template <int DIMN, class T>
qacc MvectorT<DIMN, T> operator*(const MvectorT<DIMN, T>& x, const T& y)
{
  MvectorT<DIMN, T> ret;
  for (int i = 0; i < DIMN; ++i) {
    ret.p[i] = x.p[i] * y;
  }
  return ret;
}

template <int DIMN, class T>
qacc MvectorT<DIMN, T> operator/(const MvectorT<DIMN, T>& x, const T& y)
{
  MvectorT<DIMN, T> ret;
  for (int i = 0; i < DIMN; ++i) {
    ret.p[i] = x.p[i] / y;
  }
  return ret;
}

template <int DIMN, class T>
qacc ComplexD dot_product(const MvectorT<DIMN, T>& x, const MvectorT<DIMN, T>& y)
{
  ComplexD ret = 0.0;
  for (int i = 0; i < DIMN; ++i) {
    ret += qconj(x.p[i]) * y.p[i];
  }
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
    m.p[i] = coef;
  }
}

template <int DIMN, class T>
qacc RealD qnorm(const MvectorT<DIMN, T>& m)
{
  RealD ret = 0.0;
  for (int i = 0; i < DIMN; ++i) {
    ret += qnorm(m.p[i]);
  }
  return ret;
}

template <int DIMN, class T>
qacc MvectorT<DIMN, T> vector_conjugate(const MvectorT<DIMN, T>& x)
{
  MvectorT<DIMN, T> ret;
  for (int i = 0; i < DIMN; ++i) {
    ret.p[i] = qconj(x.p[i]);
  }
  return ret;
}

template <int DIMN, class T>
std::string show(const MvectorT<DIMN, T>& m)
{
  std::ostringstream out;
  out << "[ ";
  for (int i = 0; i < DIMN; ++i) {
    out << m.p[i] << ", ";
  }
  out << "]";
  return out.str();
}

}  // namespace qlat
