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

template <int DIMN, class T>
std::string show(const MvectorT<DIMN, T>& m)
{
  std::ostringstream out;
  out << m.em();
  return out.str();
}

}  // namespace qlat
