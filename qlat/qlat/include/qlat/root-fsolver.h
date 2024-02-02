#pragma once

#include "utils.h"
#include "function.h"

#include <qlat/qlat.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include <cassert>

namespace qlat
{  //

template <class F>
void rootFSolver(double& result, double& abserr, const F& f, const double x0,
                 const double a, const double b, const double epsabs,
                 const double epsrel)
{
  gsl_function gf = make_gsl_function(f);
  gsl_root_fsolver* s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
  gsl_root_fsolver_set(s, &gf, a, b);
  assert(a <= x0 && x0 <= b);
  double x = x0, low = a, high = b;
  while (GSL_SUCCESS != gsl_root_test_interval(low, high, epsabs, epsrel)) {
    gsl_root_fsolver_iterate(s);
    x = gsl_root_fsolver_root(s);
    low = gsl_root_fsolver_x_lower(s);
    high = gsl_root_fsolver_x_upper(s);
  }
  gsl_root_fsolver_free(s);
  result = x;
  abserr = std::abs(high - low);
}

template <class F>
double rootFSolver(const F& f, const double x0, const double a, const double b,
                   const double epsabs = 0.0, const double epsrel = 1.0e-8)
{
  double result = 0.0, abserr = 0.0;
  rootFSolver(result, abserr, f, x0, a, b, epsabs, epsrel);
  return result;
}

inline void test_rootFSolver()
{
  TIMER_VERBOSE("test_rootFSolver");
  test_Functor f;
  f.scale = std::sqrt(3.0);
  DisplayInfo("", fname.c_str(), "%23.16e = %23.16e\n", 1.3467736870885985,
              rootFSolver(f, 1.0, 0.1, 1.8));
}

}  // namespace qlat
