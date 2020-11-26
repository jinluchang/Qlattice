// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <unistd.h>
#include <cassert>
#include <complex>

#ifdef OLD_CPP
#include <array-compatible.h>
#else
#include <array>
#endif

#ifdef NO_OMP
#include <omp-compatible.h>
#else
#include <omp.h>
#endif

// #define EIGEN_DONT_ALIGN_STATICALLY

#include <Eigen/Eigen>

#define QLAT_START_NAMESPACE \
  namespace qlat             \
  {
#define QLAT_END_NAMESPACE }

#ifndef USE_SINGLE_NODE
#define USE_MULTI_NODE
#endif

#define USE_NAMESPACE

#include <qutils/crc32.h>
#include <qutils/lat-io.h>
#include <qutils/qutils-io.h>
#include <qutils/qutils-vec.h>
#include <qutils/qutils.h>
#include <qutils/rng-state.h>
#include <qutils/show.h>
#include <qutils/timer.h>

namespace qlat
{  //

const int DIMN = 4;

const int NUM_COLOR = 3;

typedef Complex ComplexT;  // default Complex type

inline const std::string& cname()
{
  static const std::string s = "Qlat";
  return s;
}

inline void warn(const std::string& str = "")
{
  if (str != "") {
    displayln_info(ssprintf("WARNING: %s", str.c_str()));
  }
}

}  // namespace qlat

#include <qlat/coordinate.h>
