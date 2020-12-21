#pragma once

#include <complex>

#ifdef __NVCC__
#define QLAT_USE_ACC
#endif

#ifdef QLAT_USE_ACC
#include <thrust/complex.h>
#endif

namespace qlat
{  //

#ifdef QLAT_USE_ACC

typedef thrust::complex<double> Complex;

typedef thrust::complex<float> ComplexF;

template <class T>
thrust::complex<T> qconj(const thrust::complex<T>& x)
{
  return thrust::conj(x);
}

template <class T>
double qnorm(const thrust::complex<T>& x)
{
  return thrust::norm(x);
}

#else

typedef std::complex<double> Complex;

typedef std::complex<float> ComplexF;

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

const Complex ii(0, 1);

}  // namespace qlat
