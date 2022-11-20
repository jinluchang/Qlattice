#pragma once

#include <qlat-utils/qacc.h>
#include <complex>

#ifdef QLAT_USE_ACC
#include <thrust/complex.h>
#endif

namespace qlat
{  //

using RealD = double;

using RealF = float;

using Real = RealD;  // default Real type

#ifdef QLAT_USE_ACC

template <class T = Real>
using ComplexT = thrust::complex<T>;

template <class T>
qacc thrust::complex<T> qconj(const thrust::complex<T>& x)
{
  return thrust::conj(x);
}

template <class T>
qacc double qnorm(const thrust::complex<T>& x)
{
  return thrust::norm(x);
}

#else

template <class T = Real>
using ComplexT = std::complex<T>;

template <class T>
std::complex<T> qconj(const std::complex<T>& x)
{
  return std::conj(x);
}

template <class T>
double qnorm(const std::complex<T>& x)
{
  return std::norm(x);
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
