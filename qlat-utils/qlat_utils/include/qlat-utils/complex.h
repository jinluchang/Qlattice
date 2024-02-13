#pragma once

#include <qlat-utils/config.h>
#include <qlat-utils/qacc.h>
#include <qlat-utils/show.h>

#include <complex>
#include <type_traits>

#ifdef QLAT_USE_ACC
#include <thrust/complex.h>
#endif

#if defined QLAT_NO_ALIGNED_ALLOC
#define QLAT_ALIGNED_BYTES 1
constexpr int qlat_aligned_bytes(int size) { return 1; }
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
  }
  return ret;
}

template <class T, QLAT_ENABLE_IF(is_real<T>())>
qacc bool qisnan(const T& x)
{
  return std::isnan(x);
}

template <class T, QLAT_ENABLE_IF(is_real<T>())>
qacc RealD qnorm(const T& x)
{
  return x * x;
}

template <class T, QLAT_ENABLE_IF(is_real<T>())>
qacc T qconj(const T& x)
{
  return x;
}

#ifdef QLAT_USE_ACC

template <class T, QLAT_ENABLE_IF(is_real<T>())>
using ComplexT = thrust::complex<T>;

template <class T>
qacc ComplexT<T> qconj(const ComplexT<T>& x)
{
  return thrust::conj(x);
}

template <class T>
qacc RealD qnorm(const ComplexT<T>& x)
{
  return thrust::norm(x);
}

template <class T>
qacc ComplexT<T> qpolar(const T& r, const T& theta = T())
{
  return thrust::polar(r, theta);
}

#else

template <class T, QLAT_ENABLE_IF(is_real<T>())>
using ComplexT = std::complex<T>;

template <class T>
ComplexT<T> qconj(const ComplexT<T>& x)
{
  return std::conj(x);
}

template <class T>
RealD qnorm(const ComplexT<T>& x)
{
  return std::norm(x);
}

template <class T>
ComplexT<T> qpolar(const T& r, const T& theta = T())
{
  return std::polar(r, theta);
}

#endif

template <class T>
bool qisnan(const ComplexT<T>& arg)
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
inline std::string show(const ComplexT<T>& x)
{
  return ssprintf("(%24.17E + %24.17E j)", x.real(), x.imag());
}

}  // namespace qlat
