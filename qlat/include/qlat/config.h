// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <unistd.h>
#include <array>
#include <cassert>
#include <complex>

#ifdef NO_OMP
#include <omp-compatible.h>
#else
#include <omp.h>
#endif

#if defined NO_ALIGNED_ALLOC
#define QLAT_ALIGNED_BYTES 1
#define EIGEN_MALLOC_ALREADY_ALIGNED 0
#define EIGEN_MAX_ALIGN_BYTES 0
#define EIGEN_MAX_STATIC_ALIGN_BYTES 0
#else
#define QLAT_ALIGNED_BYTES 16 // should divide all Eigen matrix sizes (which can convert with GPT).
#ifndef EIGEN_MALLOC_ALREADY_ALIGNED
#define EIGEN_MALLOC_ALREADY_ALIGNED 1
#endif
#ifndef EIGEN_MAX_ALIGN_BYTES
#define EIGEN_MAX_ALIGN_BYTES QLAT_ALIGNED_BYTES
#endif
#ifndef EIGEN_MAX_STATIC_ALIGN_BYTES
#define EIGEN_MAX_STATIC_ALIGN_BYTES QLAT_ALIGNED_BYTES
#endif
#endif

#ifndef USE_SINGLE_NODE
#define USE_MULTI_NODE
#endif

#define USE_NAMESPACE

#include <qlat-utils/crc32.h>
#include <qlat-utils/eigen.h>
#include <qlat-utils/lat-io.h>
#include <qlat-utils/qacc-func.h>
#include <qlat-utils/qacc.h>
#include <qlat-utils/qutils-io.h>
#include <qlat-utils/qutils-vec.h>
#include <qlat-utils/qutils.h>
#include <qlat-utils/rng-state.h>
#include <qlat-utils/show.h>
#include <qlat-utils/timer.h>
#include <qlat-utils/vector.h>

namespace qlat
{  //

const int DIMN = 4;

const int NUM_COLOR = 3;

typedef Complex ComplexT;  // default Complex type

inline void warn(const std::string& str = "")
{
  if (str != "") {
    displayln_info(ssprintf("WARNING: %s", str.c_str()));
  }
}

}  // namespace qlat

#include <qlat/coordinate.h>
