#pragma once

#include <complex>

#ifdef __NVCC__
#define QLAT_USE_GPU
#endif

#ifdef QLAT_USE_GPU
#include <thrust/complex.h>
#endif

namespace qlat
{  //

#ifdef QLAT_USE_GPU

typedef std::thrust<double> Complex;

typedef std::thrust<float> ComplexF;

#else

typedef std::complex<double> Complex;

typedef std::complex<float> ComplexF;

#endif

const Complex ii(0, 1);

}  // namespace qlat
