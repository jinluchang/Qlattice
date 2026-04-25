#pragma once

#include <gsl/gsl_math.h>
#include <qlat/qlat.h>

#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <string>

inline bool is_very_close(const double x, const double y,
                          const double epsabs = 1.0e-10,
                          const double epsrel = 1.0e-10)
{
  const double diff = std::abs(x - y);
  if (diff <= epsabs) {
    return true;
  } else {
    return diff * 2.0 <= (std::abs(x) + std::abs(y)) * epsrel;
  }
}
