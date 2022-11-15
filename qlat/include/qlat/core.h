#pragma once

#include <qlat-utils/timer.h>
#include <qlat-utils/complex.h>

#include <qlat/config.h>

#if defined NO_ALIGNED_ALLOC
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

// #define ALIGN alignas(QLAT_ALIGNED_BYTES)
#define ALIGN __attribute__((aligned(QLAT_ALIGNED_BYTES)))

#ifndef USE_SINGLE_NODE
#define USE_MULTI_NODE
#endif

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

struct API Coordinate : public array<int, DIMN> {
  //
  qacc Coordinate() { array<int, DIMN>::fill(0); }
  //
  qacc Coordinate(int first, int second, int third, int fourth)
  {
    int* p = data();
    p[0] = first;
    p[1] = second;
    p[2] = third;
    p[3] = fourth;
  }
  //
  qacc long product() const
  {
    long ret = 1;
    int size_ = size();
    for (int i = 0; i < size_; i++) {
      ret *= operator[](i);
    }
    return ret;
  }
};

}  // namespace qlat
