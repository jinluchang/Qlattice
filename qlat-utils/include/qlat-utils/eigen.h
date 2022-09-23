#pragma once

#include <qlat-utils/complex.h>

#ifdef QLAT_GRID

#include <Grid/Eigen/Eigen>

#else

#include <Eigen/Eigen>

#endif

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
