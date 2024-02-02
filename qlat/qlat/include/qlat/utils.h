#pragma once

#include <gsl/gsl_math.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <string>
#include <qlat/qlat.h>

inline bool is_very_close(const double x, const double y, const double epsabs = 1.0e-10, const double epsrel = 1.0e-10)
{
  const double diff = std::abs(x - y);
  if (diff <= epsabs) {
    return true;
  } else {
    return diff * 2.0 <= (std::abs(x) + std::abs(y)) * epsrel;
  }
}
