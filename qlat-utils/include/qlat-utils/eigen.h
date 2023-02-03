#pragma once

#include <qlat-utils/config.h>
#include <qlat-utils/complex.h>

// -------------------------------------------------------------------------------------

#if defined QLAT_NO_ALIGNED_ALLOC

#define QLAT_ALIGNED_BYTES 1
// #define EIGEN_MALLOC_ALREADY_ALIGNED 0
#define EIGEN_MAX_ALIGN_BYTES 0
#define EIGEN_MAX_STATIC_ALIGN_BYTES 0

#else

#define QLAT_ALIGNED_BYTES 16 // should divide all Eigen matrix sizes (which can convert with GPT).
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

#define ALIGN __attribute__((aligned(QLAT_ALIGNED_BYTES)))
// #define ALIGN alignas(QLAT_ALIGNED_BYTES)

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

template <typename RealType>
struct NumTraits<thrust::complex<RealType> >
    : NumTraits<std::complex<RealType> > {
  typedef RealType Real;
  typedef typename NumTraits<Real>::Literal Literal;
  enum {
    IsComplex = 1,
    RequireInitialization = NumTraits<Real>::RequireInitialization,
    ReadCost = 2 * NumTraits<Real>::ReadCost,
    AddCost = 2 * NumTraits<Real>::AddCost,
    MulCost = 4 * NumTraits<Real>::MulCost + 2 * NumTraits<Real>::AddCost
  };
};

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
