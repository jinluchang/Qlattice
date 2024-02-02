#pragma once

#include "utils.h"
#include "function.h"

#include <qlat/qlat.h>

#include <gsl/gsl_integration.h>

#include <iostream>
#include <cmath>

namespace qlat
{  //

template <class F>
int integrate(double& result, double& abserr, const F& f, const double a,
              const double b, const double epsabs, const double epsrel)
{
  const int limit = 1000;
  const int key = 2;
  gsl_function gf = make_gsl_function(f);
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
  int ret = gsl_integration_qag(&gf, a, b, epsabs, epsrel, limit, key, w,
                                &result, &abserr);
  gsl_integration_workspace_free(w);
  return ret;
}

template <class F>
double integrate(const F& f, const double a, const double b,
                 const double epsabs = 0.0, const double epsrel = 1.0e-8)
{
  double result = 0.0, abserr = 0.0;
  integrate(result, abserr, f, a, b, epsabs, epsrel);
  return result;
}

inline void test_integrate()
{
  TIMER_VERBOSE("test_integrate");
  test_Functor f;
  f.scale = std::sqrt(3.0);
  DisplayInfo("", fname.c_str(), "2 = %23.16e\n",
              integrate(f, 0, std::sqrt(PI / f.scale)));
}

}  // namespace qlat
