#pragma once

#include <qlat-utils/config.h>
#include <qlat-utils/qacc.h>
#include <qlat-utils/show.h>

#include <complex>

#ifdef QLAT_USE_ACC
#include <thrust/complex.h>
#endif

#if defined QLAT_NO_ALIGNED_ALLOC
#define QLAT_ALIGNED_BYTES 1
constexpr int qlat_aligned_bytes(int size) { return 1; }
#define ALIGN(SIZE) __attribute__((aligned(1)))
// #define ALIGN(SIZE) alignas(1)
#else
#define QLAT_ALIGNED_BYTES 16 // should divide all matrix sizes (which can convert with GPT).
constexpr int qlat_aligned_bytes(int size)
{
  return (size % 2 != 0)     ? 1
         : (size % 4 != 0)   ? 2
         : (size % 8 != 0)   ? 4
         : (size % 16 != 0)  ? 8
         : (size % 32 != 0)  ? 16
         : (size % 64 != 0)  ? 32
         : (size % 128 != 0) ? 64
         : (size % 256 != 0) ? 128
                             : 256;
}
#define ALIGN(SIZE) __attribute__((aligned(qlat_aligned_bytes(SIZE))))
// #define ALIGN(SIZE) alignas(SIZE)
#endif

namespace qlat
{  //

using RealD = double;

using RealF = float;

using Real = RealD;  // default Real type

qacc RealD qnorm(const RealD& x) { return x * x; }

qacc RealF qnorm(const RealF& x) { return x * x; }

qacc RealD qconj(const RealD& x) { return x; }

qacc RealF qconj(const RealF& x) { return x; }

#ifdef QLAT_USE_ACC

template <class T = Real>
using ComplexT = thrust::complex<T>;

template <class T>
qacc ComplexT<T> qconj(const ComplexT<T>& x)
{
  return thrust::conj(x);
}

template <class T>
qacc T qnorm(const ComplexT<T>& x)
{
  return thrust::norm(x);
}

template <class T>
qacc ComplexT<T> qpolar(const T& r, const T& theta = T())
{
  return thrust::polar(r, theta);
}

#else

template <class T = Real>
using ComplexT = std::complex<T>;

template <class T>
ComplexT<T> qconj(const ComplexT<T>& x)
{
  return std::conj(x);
}

template <class T>
T qnorm(const ComplexT<T>& x)
{
  return std::norm(x);
}

template <class T>
qacc ComplexT<T> qpolar(const T& r, const T& theta = T())
{
  return std::polar(r, theta);
}

#endif

using ComplexD = ComplexT<RealD>;

using ComplexF = ComplexT<RealF>;

using Complex = ComplexT<Real>;

const ComplexD ii(0, 1);

template <class T>
inline std::string show(const ComplexT<T>& x)
{
  return ssprintf("(%24.17E + %24.17E j)", x.real(), x.imag());
}

}  // namespace qlat
