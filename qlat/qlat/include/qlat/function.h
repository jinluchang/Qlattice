#pragma once

#include <gsl/gsl_math.h>

#include <cmath>

namespace qlat
{  //

template <class F>
inline double gsl_util_function(double x, void* params)
{
  return (*((const F*)params))(x);
}

template <class F>
inline gsl_function make_gsl_function(const F& f)
{
  gsl_function gf;
  gf.function = &gsl_util_function<F>;
  gf.params = (void*)&f;
  return gf;
}

struct test_Functor {
  double scale;
  double operator()(const double x) const
  {
    return scale * std::sin(scale * sqr(x)) * 2.0 * x;
  }
};

}  // namespace qlat
