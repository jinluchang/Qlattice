#pragma once

#include <qlat-utils/config.h>
#include <qlat-utils/qacc.h>
#include <qlat-utils/show.h>
#include <qlat-utils/utils-ldouble.h>

#include <complex>
#include <type_traits>

#ifdef QLAT_USE_ACC
#include <thrust/complex.h>
#endif

#if defined QLAT_NO_ALIGNED_ALLOC
#define QLAT_ALIGNED_BYTES 1
constexpr Int qlat_aligned_bytes(Int size) { (void)size; return 1; }
#define QLAT_ALIGN(SIZE) __attribute__((aligned(1)))
// #define QLAT_ALIGN(SIZE) alignas(1)
#else
#define QLAT_ALIGNED_BYTES 16 // should divide all matrix sizes (which can convert with GPT).
constexpr int qlat_aligned_bytes(int size)
{
  int ret = 0;
  if (size % 2 != 0) {
    ret = 1;
  } else if (size % 4 != 0) {
    ret = 2;
  } else if (size % 8 != 0) {
    ret = 4;
  } else if (size % 16 != 0) {
    ret = 8;
  } else if (size % 32 != 0) {
    ret = 16;
  } else if (size % 64 != 0) {
    ret = 32;
  } else if (size % 128 != 0) {
    ret = 64;
  } else if (size % 256 != 0) {
    ret = 128;
  } else {
    ret = 256;
  }
  return ret;
}
#define QLAT_ALIGN(SIZE) __attribute__((aligned(qlat_aligned_bytes(SIZE))))
// #define QLAT_ALIGN(SIZE) alignas(SIZE)
#endif

#define QLAT_ENABLE_IF(X) std::enable_if_t<(X), bool> = true

namespace qlat
{  //

template <class M, class N>
qacc constexpr bool is_same()
{
  return std::is_same<M, N>::value;
}

template <class M>
qacc constexpr bool is_real()
{
  bool ret = false;
  if (is_same<M, RealD>()) {
    ret = true;
  } else if (is_same<M, RealF>()) {
    ret = true;
  } else if (is_same<M, RealDD>()) {
    ret = true;
  }
  return ret;
}

template <class T, QLAT_ENABLE_IF(is_real<T>())>
qacc bool qisnan(const T &x)
{
  return std::isnan(x);
}

template <class T, QLAT_ENABLE_IF(is_real<T>())>
qacc RealD qnorm(const T &x)
{
  return x * x;
}

template <class M, class N, QLAT_ENABLE_IF(is_real<M>() and is_real<N>())>
qacc RealD qnorm(const M &x, const N &y)
{
  return x * y;
}

template <class T, QLAT_ENABLE_IF(is_real<T>())>
qacc T qconj(const T &x)
{
  return x;
}

#ifdef QLAT_USE_ACC

template <class T, QLAT_ENABLE_IF(is_real<T>())>
using ComplexT = thrust::complex<T>;

template <class T>
qacc ComplexT<T> qconj(const ComplexT<T> &x)
{
  return thrust::conj(x);
}

template <class T>
qacc RealD qnorm(const ComplexT<T> &x)
{
  return thrust::norm(x);
}

template <class T>
qacc ComplexT<T> qpolar(const T &r, const T &theta = T())
{
  return thrust::polar(r, theta);
}

#else

template <class T, QLAT_ENABLE_IF(is_real<T>())>
using ComplexT = std::complex<T>;

template <class T>
ComplexT<T> qconj(const ComplexT<T> &x)
{
  return std::conj(x);
}

template <class T>
RealD qnorm(const ComplexT<T> &x)
{
  return std::norm(x);
}

template <class T>
ComplexT<T> qpolar(const T &r, const T &theta = T())
{
  return std::polar(r, theta);
}

#endif

template <class T1, class T2>
qacc RealD qnorm(const ComplexT<T1> &x, const ComplexT<T2> &y)
{
  return qnorm(x.real(), y.real()) + qnorm(x.imag(), y.imag());
}

template <class T>
qacc bool qisnan(const ComplexT<T> &arg)
{
  return qisnan(arg.real()) or qisnan(arg.imag());
}

using ComplexD = ComplexT<RealD>;

using ComplexF = ComplexT<RealF>;

using Complex = ComplexT<Real>;

template <class M>
qacc constexpr bool is_complex()
{
  bool ret = false;
  if (is_same<M, ComplexD>()) {
    ret = true;
  } else if (is_same<M, ComplexF>()) {
    ret = true;
  }
  return ret;
}

const ComplexD ii(0, 1);

template <class T>
inline std::string show(const ComplexT<T> &x)
{
  return ssprintf("(%24.17E + %24.17E j)", x.real(), x.imag());
}

qacc ComplexT<RealD> operator*(const RealD &a, const ComplexT<RealD> &b)
{
  return ComplexT<RealD>(b.real() * a, b.imag() * a);
}

qacc ComplexT<RealD> operator*(const ComplexT<RealD> &b, const RealD &a)
{
  return ComplexT<RealD>(b.real() * a, b.imag() * a);
}

qacc ComplexT<RealF> operator*(const RealD &a, const ComplexT<RealF> &b)
{
  return ComplexT<RealF>(b.real() * a, b.imag() * a);
}

qacc ComplexT<RealF> operator*(const ComplexT<RealF> &b, const RealD &a)
{
  return ComplexT<RealF>(b.real() * a, b.imag() * a);
}

qacc ComplexT<RealF> operator*(const RealF &a, const ComplexT<RealF> &b)
{
  return ComplexT<RealF>(b.real() * a, b.imag() * a);
}

qacc ComplexT<RealF> operator*(const ComplexT<RealF> &b, const RealF &a)
{
  return ComplexT<RealF>(b.real() * a, b.imag() * a);
}

qacc ComplexT<RealDD> operator*(const RealD &a, const ComplexT<RealDD> &b)
{
  return ComplexT<RealDD>(b.real() * RealDD(a), b.imag() * RealDD(a));
}

qacc ComplexT<RealDD> operator*(const ComplexT<RealDD> &a,
                                const ComplexT<RealDD> &b)
{
  return ComplexT<RealDD>(a.real() * b.real() - a.imag() * b.imag(),
                          a.imag() * b.real() + a.real() * b.imag());
}

qacc ComplexT<RealDD> operator*(const RealDD &a, const ComplexT<RealDD> &b)
{
  return ComplexT<RealDD>(b.real() * a, b.imag() * a);
}

qacc ComplexT<RealDD> operator*(const ComplexT<RealDD> &b, const RealD &a)
{
  return ComplexT<RealDD>(b.real() * RealDD(a), b.imag() * RealDD(a));
}

qacc ComplexT<RealDD> operator*(const ComplexT<RealDD> &b, const RealDD &a)
{
  return ComplexT<RealDD>(b.real() * a, b.imag() * a);
}

template <class Ta, class Tb>
qacc void copy_complex(ComplexT<Ta> &r, const ComplexT<Tb> &a)
{
  r = a;
  // r = ComplexT<Ta>(a.real(), a.imag());
}

template <class Ta>
qacc void copy_complex(ComplexT<Ta> &r, const ComplexT<RealDD> &a)
{
  r = a;
}

template <class Ta>
qacc Ta qfabs(Ta a)
{
  return std::fabs(a);
}

template <>
qacc RealDD qfabs(RealDD a)
{
  return fabsT(a);
}

template <class Ta>
qacc Ta qsqrt(Ta a)
{
  return std::sqrt(a);
}

template <>
qacc RealDD qsqrt(RealDD a)
{
  return sqrtT(a);
}

template <class Ta>
qacc Ta qfmin(Ta a, Ta b)
{
  return fmin(a, b);
}

template <>
qacc RealDD qfmin(RealDD a, RealDD b)
{
  if (a.Y() < b.Y()) {
    return a;
  }
  if (a.Y() > b.Y()) {
    return b;
  }
  if (a.X() < b.X()) {
    return a;
  }
  if (a.X() > b.X()) {
    return b;
  }
  return b;
}

template <class Ta>
qacc Ta qsin(Ta a)
{
  return sin(a);
}

template <class Ta>
qacc Ta qcos(Ta a)
{
  return cos(a);
}

template <class Ta>
qacc Ta qacos(Ta a)
{
  return acos(a);
}

template <>
qacc RealDD qsin(RealDD a)
{
  RealD y = a.Y();
  RealD x = a.X();
  RealDD z, e;
  RealDD cosy = std::cos(y);
  RealDD cosx = std::cos(x);
  RealDD siny = std::sin(y);
  RealDD sinx = std::sin(x);
  e = siny * cosx + cosy * sinx;
  z.Y() = e.X() + e.Y();
  z.X() = 0.0;
  z.X() = e - z;
  return z;
  // return std::sin(y);
}

template <>
qacc RealDD qcos(RealDD a)
{
  RealD y = a.Y();
  RealD x = a.X();
  RealDD z, e;
  RealDD cosy = std::cos(y);
  RealDD cosx = std::cos(x);
  RealDD siny = std::sin(y);
  RealDD sinx = std::sin(x);
  e = cosy * cosx - siny * sinx;
  z.Y() = e.X() + e.Y();
  z.X() = 0.0;
  z.X() = e - z;
  return z;
  // return std::cos(y);
}

template <>
qacc RealDD qacos(RealDD a)
{
  RealD y = a.Y();
  RealD x = a.X();
  RealDD z,e;
  e.Y() = std::acos(y);
  e.X() = -1.0 * (1.0 / std::sqrt(1.0 - y * y)) * x;
  z.Y() = e.X() + e.Y();
  z.X() = 0.0;
  z.X() = e - z;
  return z;
  // return std::acos(y);
}

// qacc std::complex<double>& std::complex<double>::operator=(const
// std::complex<qlat::RealDD>& a) {
//   this->real() = a.real().Y();
//   this->imag() = a.imag().Y();
//   return *this;
// }

// qacc ComplexT<float>& ComplexT<float>::operator=(const ComplexT<RealDD> &a) {
//   return ComplexT<float>(a.real().Y(), a.imag().Y());
// }
//
// qacc ComplexT<RealDD>& ComplexT<RealDD>::operator=(const ComplexT<double> &a)
// {
//   return ComplexT<RealDD>(a.real(), a.imag());
// }
//
// qacc ComplexT<RealDD>& ComplexT<RealDD>::operator=(const ComplexT<float> &a)
// {
//   return ComplexT<RealDD>(a.real(), a.imag());
// }

qacc ComplexT<RealDD> operator/(const ComplexT<RealDD> &a,
                                const ComplexT<RealDD> &b)
{
  RealDD sq = b.real() * b.real() + b.imag() * b.imag();
  RealDD r = a.real() * b.real() + a.imag() * b.imag();
  RealDD i = a.imag() * b.real() - a.real() * b.imag();
  return ComplexT<RealDD>(r / sq, i / sq);
}

qacc ComplexT<RealDD> operator+(const ComplexT<RealDD> &a,
                                const ComplexT<RealDD> &b)
{
  return ComplexT<RealDD>(a.real() + b.real(), a.imag() + b.imag());
}

qacc ComplexT<RealDD> operator-(const ComplexT<RealDD> &a,
                                const ComplexT<RealDD> &b)
{
  return ComplexT<RealDD>(a.real() - b.real(), a.imag() - b.imag());
}

qacc RealDD qnorm(const ComplexT<RealDD> &a)
{
  RealDD sq = a.real() * a.real() + a.imag() * a.imag();
  return sq;
}

template <class T>
qacc ComplexT<RealDD> qconj(const ComplexT<RealDD> &x)
{
  // RealDD tmp = minus(x.imag());
  return ComplexT<RealDD>(x.real(), -x.imag());
}

}  // namespace qlat

#ifndef __HIPCC__

namespace std
{

template <>
qacc complex<qlat::RealD> &complex<qlat::RealD>::operator=(
    const complex<qlat::RealDD> &__z)
{
  *this = complex<qlat::RealD>(__z.real(), __z.imag());
  return *this;
}

}  // namespace std

#endif
