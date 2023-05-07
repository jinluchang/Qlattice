#pragma once

#include <qlat-utils/config.h>
#include <qlat-utils/complex.h>

// -------------------------------------------------------------------------------------

#if defined QLAT_NO_ALIGNED_ALLOC
// #define EIGEN_MALLOC_ALREADY_ALIGNED 0
#define EIGEN_MAX_ALIGN_BYTES 0
#define EIGEN_MAX_STATIC_ALIGN_BYTES 0
#else
// #ifndef EIGEN_MALLOC_ALREADY_ALIGNED
// #define EIGEN_MALLOC_ALREADY_ALIGNED 1
// #endif
#ifndef EIGEN_MAX_ALIGN_BYTES
#define EIGEN_MAX_ALIGN_BYTES QLAT_ALIGNED_BYTES
#endif
#ifndef EIGEN_MAX_STATIC_ALIGN_BYTES
#define EIGEN_MAX_STATIC_ALIGN_BYTES QLAT_ALIGNED_BYTES
#endif
#endif

// -------------------------------------------------------------------------------------

#ifdef QLAT_USE_GRID_EIGEN

#include <Grid/Eigen/Eigen>

#else

#include <Eigen/Eigen>

#endif

// -------------------------------------------------------------------------------------

#ifdef QLAT_USE_ACC

namespace Eigen
{  //

#ifndef QLAT_GRID

template <>
struct NumTraits<thrust::complex<float>> : NumTraits<std::complex<float>> {
  typedef float Real;
  typedef typename NumTraits<Real>::Literal Literal;
  enum {
    IsComplex = 1,
    RequireInitialization = NumTraits<Real>::RequireInitialization,
    ReadCost = 2 * NumTraits<Real>::ReadCost,
    AddCost = 2 * NumTraits<Real>::AddCost,
    MulCost = 4 * NumTraits<Real>::MulCost + 2 * NumTraits<Real>::AddCost
  };
};

template <>
struct NumTraits<thrust::complex<double>> : NumTraits<std::complex<double>> {
  typedef double Real;
  typedef typename NumTraits<Real>::Literal Literal;
  enum {
    IsComplex = 1,
    RequireInitialization = NumTraits<Real>::RequireInitialization,
    ReadCost = 2 * NumTraits<Real>::ReadCost,
    AddCost = 2 * NumTraits<Real>::AddCost,
    MulCost = 4 * NumTraits<Real>::MulCost + 2 * NumTraits<Real>::AddCost
  };
};

#endif

}  // namespace Eigen

namespace thrust
{  //

template <typename T>
EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE T real(complex<T> const& x)
{
  return x.real();
}

template <typename T>
EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE T imag(complex<T> const& x)
{
  return x.imag();
}

}  // namespace thrust

#endif

// -------------------------------------------------------------------------------------
