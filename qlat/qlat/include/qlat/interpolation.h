#pragma once

#include "utils.h"
#include "function.h"
#include "root-fsolver.h"

#include <qlat/qlat.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>

#include <vector>
#include <cassert>

namespace qlat
{  //

struct Interpolation {
  gsl_spline* spline;
  //
  Interpolation() { init(); }
  template <class F>
  Interpolation(const F& f, const RealD high, const RealD low, const Int size)
  {
    init();
    init(f, high, low, size);
  }
  ~Interpolation() { free(); }
  //
  void init() { spline = NULL; }
  template <class F>
  void init(const F& f, const RealD high, const RealD low, const Int size)
  {
    free();
    spline = gsl_spline_alloc(gsl_interp_cspline, size);
    const RealD sep = (high - low) / (size - 1);
    std::vector<RealD> xs(size);
    std::vector<RealD> ys(size);
    for (Int i = 0; i < size; ++i) {
      const RealD x = low + sep * i;
      xs[i] = x;
    }
    xs[size - 1] = high;
    for (Int i = 0; i < size; ++i) {
      const RealD x = xs[i];
      ys[i] = f(x);
    }
    gsl_spline_init(spline, xs.data(), ys.data(), size);
  }
  void free()
  {
    if (NULL != spline) {
      gsl_spline_free(spline);
    }
  }
  //
  RealD operator()(const RealD x) const
  {
    return gsl_spline_eval(spline, x, NULL);
  }
};

struct Interpolation2d {
  gsl_spline2d* spline;
  //
  Interpolation2d() { init(); }
  template <class F>
  Interpolation2d(const F& f, const RealD xhigh, const RealD xlow,
                  const Long xsize, const RealD yhigh, const RealD ylow,
                  const Long ysize)
  {
    init();
    init(f, xhigh, xlow, xsize, yhigh, ylow, ysize);
  }
  ~Interpolation2d() { free(); }
  //
  void init() { spline = NULL; }
  template <class F>
  void init(const F& f, const RealD xhigh, const RealD xlow, const Long xsize,
            const RealD yhigh, const RealD ylow, const Long ysize)
  {
    free();
    spline = gsl_spline2d_alloc(gsl_interp2d_bicubic, xsize, ysize);
    const RealD xsep = (xhigh - xlow) / (xsize - 1);
    const RealD ysep = (yhigh - ylow) / (ysize - 1);
    const Long size = xsize * ysize;
    std::vector<RealD> xs(xsize, 0.0);
    std::vector<RealD> ys(ysize, 0.0);
    std::vector<RealD> zs(size, 0.0);
    for (Long i = 0; i < xsize; ++i) {
      const RealD x = xlow + xsep * i;
      xs[i] = x;
    }
    for (Long j = 0; j < ysize; ++j) {
      const RealD y = ylow + ysep * j;
      ys[j] = y;
    }
    xs[xsize - 1] = xhigh;
    ys[ysize - 1] = yhigh;
    for (Long i = 0; i < xsize; ++i) {
      for (Long j = 0; j < ysize; ++j) {
        const RealD x = xs[i];
        const RealD y = ys[j];
        gsl_spline2d_set(spline, zs.data(), i, j, f(x, y));
      }
    }
    gsl_spline2d_init(spline, xs.data(), ys.data(), zs.data(), xsize, ysize);
  }
  void free()
  {
    if (NULL != spline) {
      gsl_spline2d_free(spline);
    }
  }
  //
  RealD operator()(const RealD x, const RealD y) const
  {
    return gsl_spline2d_eval(spline, x, y, NULL, NULL);
  }
};

inline void test_interpolation()
{
  TIMER_VERBOSE("test_interpolation");
  test_Functor f;
  f.scale = std::sqrt(3.0);
  Interpolation fi(f, 10.0, -10.0, 16 * 1024);
  DisplayInfo("", fname.c_str(), "%23.16e = %23.16e\n", 1.3467736870885985,
              rootFSolver(f, 1.0, 0.1, 1.8));
  DisplayInfo("", fname.c_str(), "%23.16e = %23.16e\n", 1.3467736870885985,
              rootFSolver(fi, 1.0, 0.1, 1.8));
}

}  // namespace qlat
