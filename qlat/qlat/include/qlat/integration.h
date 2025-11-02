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
Int integrate(RealD& result, RealD& abserr, const F& f, const RealD a,
              const RealD b, const RealD epsabs, const RealD epsrel)
{
  const Int limit = 1000;
  const Int key = 2;
  gsl_function gf = make_gsl_function(f);
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
  Int ret = gsl_integration_qag(&gf, a, b, epsabs, epsrel, limit, key, w,
                                &result, &abserr);
  gsl_integration_workspace_free(w);
  return ret;
}

template <class F>
RealD integrate(const F& f, const RealD a, const RealD b,
                 const RealD epsabs = 0.0, const RealD epsrel = 1.0e-8)
{
  RealD result = 0.0, abserr = 0.0;
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
