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

#define EIGEN_MALLOC_ALREADY_ALIGNED 1
#define EIGEN_MAX_ALIGN_BYTES 256
#define EIGEN_MAX_STATIC_ALIGN_BYTES 256

#if defined NO_ALIGNED_ALLOC
#define EIGEN_MALLOC_ALREADY_ALIGNED 0
#define EIGEN_MAX_ALIGN_BYTES 0
#define EIGEN_MAX_STATIC_ALIGN_BYTES 0
#endif

#ifndef USE_SINGLE_NODE
#define USE_MULTI_NODE
#endif

#define USE_NAMESPACE

#include <qutils/crc32.h>
#include <qutils/eigen.h>
#include <qutils/lat-io.h>
#include <qutils/qacc-func.h>
#include <qutils/qacc.h>
#include <qutils/qutils-io.h>
#include <qutils/qutils-vec.h>
#include <qutils/qutils.h>
#include <qutils/rng-state.h>
#include <qutils/show.h>
#include <qutils/timer.h>
#include <qutils/vector.h>

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
