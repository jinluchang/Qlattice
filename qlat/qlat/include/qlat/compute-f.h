#pragma once

#include "muon-line-config.h"
#include "utils.h"
#include "integration.h"
#include "interpolation.h"

#include <qlat/qlat.h>

#include <gsl/gsl_specfunc.h>

#include <iostream>
#include <cmath>
#include <cassert>

namespace qlat
{  //

inline RealD specialK0(const RealD x) { return gsl_sf_bessel_K0(x); }

inline RealD specialK0x(const RealD x)
{
  if (0.0 == x) {
    return 0.0;
  } else {
    return x * specialK0(x);
  }
}

inline RealD specialK1(const RealD x) { return gsl_sf_bessel_K1(x); }

inline RealD specialK1x(const RealD x)
{
  if (0.0 == x) {
    return 1.0;
  } else {
    return x * specialK1(x);
  }
}

struct FgxxIntegrand {
  RealD len, ratio;
  //
  FgxxIntegrand() { init(); }
  //
  void init()
  {
    len = 0.0;
    ratio = 0.0;
  }
  void init(const RealD len_, const RealD ratio_)
  {
    len = len_;
    ratio = ratio_;
    assert(len >= 0.0);
    assert(-1.0 <= ratio && ratio <= 1.0);
  }
  //
  RealD operator()(const RealD y) const
  {
    if (0.0 == y || 0.0 == len || 0.0 == ratio) {
      return 0.0;
    } else {
      return y * (std::exp(-sqr(y) * ratio * len) - 1.0) *
             specialK0(sqr(y) * len);
    }
  }
};

inline RealD gxxCalcIntegrate(const RealD len, const RealD ratio)
{
  FgxxIntegrand gxx;
  gxx.init(len, ratio);
  RealD ret = 1.0 / (4.0 * sqr(PI)) * integrate(gxx, 0.0, 1.0, 0.0, 1.0e-12);
  if (std::isnan(ret)) {
    displayln(ssprintf("%e %e %e %d", len, ratio, ret, std::isnan(ret)));
    assert(false);
  }
  return ret;
}

inline void recordGxx(const RealD len, const RealD lenp, const RealD ratio,
                      const RealD gxx)
{
  displayln(ssprintf("gxx: %s %s %s %s", show(len).c_str(), show(lenp).c_str(),
                     show(ratio).c_str(), show(gxx).c_str()));
}

inline RealD gxxInterpPrime(const RealD lenp, const RealD psi)
{
  if (1.0 == lenp || PI / 2.0 == psi) {
    return 0.0;
  } else {
    const RealD len = lenp / (1.0 - lenp);
    const RealD ratio = std::cos(psi);
    const RealD gxx = gxxCalcIntegrate(len, ratio);
    // recordGxx(len, lenp, ratio, gxx);
    return gxx;
  }
}

inline RealD gxxCalc(const RealD len, const RealD ratio)
{
  if (false == IS_USING_INTERPOLATION) {
    return gxxCalcIntegrate(len, ratio);
  } else {
    // ADJUST ME
    const Int nlen = 512 + 1;
    const Int nratio = 128 + 1;
    static Interpolation2d fi(gxxInterpPrime, 1.0, 0.0, nlen, PI, 0.0, nratio);
    const RealD lenp = len / (1.0 + len);
    const RealD psi = std::acos(ratio);
    return fi(lenp, psi);
  }
}

inline RealD gxxCalc(const qlat::CoordinateD& c)
{
  const RealD len = coordinate_len(c);
  if (0.0 == len) {
    return 0.0;
  } else {
    RealD ratio = c[3] / len;
    if (ratio > 1.0) {
      ratio = 1.0;
    } else if (ratio < -1.0) {
      ratio = -1.0;
    }
    const RealD ans = gxxCalc(len, ratio);
    if (std::isnan(ans)) {
      displayln(ssprintf("%e %e %e %e %d", len, ratio, std::acos(ratio), ans,
                         std::isnan(ans)));
      assert(false);
    }
    return ans;
  }
}

struct FhxIntegrand {
  RealD len;
  //
  FhxIntegrand() { init(); }
  //
  void init() { len = 0.0; }
  void init(const RealD len_)
  {
    len = len_;
    assert(len >= 0.0);
    if (0.0 == len) {
      len = 1.0e-5;  // ADJUST ME
    }
  }
  //
  RealD operator()(const RealD y) const
  {
    RealD y2len = sqr(y) * len;
    if (0.0 == y || y2len > 128.0) {
      return 0.0;
    } else {
      return y * specialK0(y2len);
    }
  }
};

inline RealD hxCalcIntegrate(const RealD len)
{
  FhxIntegrand hx;
  hx.init(len);
  const RealD ans =
      1.0 / (4.0 * sqr(PI)) * integrate(hx, 0.0, 1.0, 0.0, 1.0e-12);
  assert(false == isnan(ans));
  return ans;
}

inline void recordHx(const RealD len, const RealD lenp, const RealD hx)
{
  displayln(ssprintf("hx: %s %s %s", show(len).c_str(), show(lenp).c_str(),
                     show(hx).c_str()));
}

inline RealD hxInterpPrime(const RealD lenp)
{
  if (1.0 == lenp) {
    return 0.0;
  } else {
    const RealD len = lenp / (1.0 - lenp);
    const RealD hx = hxCalcIntegrate(len);
    // recordHx(len, lenp, hx);
    return hx;
  }
}

inline RealD hxCalc(const RealD len)
{
  if (false == IS_USING_INTERPOLATION) {
    return hxCalcIntegrate(len);
  } else {
    // ADJUST ME
    const Int nlen = 16 * 1024 + 1;
    static Interpolation fi(hxInterpPrime, 1.0, 0.0, nlen);
    const RealD lenp = len / (1.0 + len);
    return fi(lenp);
  }
}

inline RealD hxCalc(const qlat::CoordinateD& c)
{
  const RealD ans = hxCalc(coordinate_len(c));
  assert(false == std::isnan(ans));
  return ans;
}

inline RealD fxxCalc(const qlat::CoordinateD& c)
{
  // ADJUST ME
  return hxCalc(c) + gxxCalc(c);
  // return hxCalc(c);
  // return gxxCalc(c);
}

struct Ff1xxIntegrand {
  RealD len, ratio;
  //
  Ff1xxIntegrand() { init(); }
  //
  void init()
  {
    len = 0.0;
    ratio = 0.0;
  }
  void init(const RealD len_, const RealD ratio_)
  {
    len = len_;
    ratio = ratio_;
    assert(len >= 0.0);
    assert(-1.0 <= ratio && ratio <= 1.0);
  }
  //
  RealD operator()(const RealD y) const
  {
    return std::exp(-y * ratio * len) * specialK1x(y * len);
  }
};

inline RealD f1xxCalcIntegrate(const RealD len, const RealD ratio)
{
  Ff1xxIntegrand f1xx;
  f1xx.init(len, ratio);
  return -1.0 / (8.0 * sqr(PI)) * integrate(f1xx, 0.0, 1.0, 0.0, 1.0e-12);
}

inline void recordF1xx(const RealD len, const RealD lenp, const RealD ratio,
                       const RealD f1xx)
{
  displayln(ssprintf("f1xx: %s %s %s %s", show(len).c_str(), show(lenp).c_str(),
                     show(ratio).c_str(), show(f1xx).c_str()));
}

inline RealD f1xxInterpPrime(const RealD lenp, const RealD psi)
{
  if (1.0 == lenp) {
    return 0.0;
  } else {
    const RealD len = lenp / (1.0 - lenp);
    const RealD ratio = std::cos(psi);
    const RealD f1xx = f1xxCalcIntegrate(len, ratio);
    // recordF1xx(len, lenp, ratio, f1xx);
    return f1xx;
  }
}

inline RealD f1xxCalc(const RealD len, const RealD ratio)
{
  if (0.0 == len) {
    return 1000.0;  // ADJUST ME
  } else if (false == IS_USING_INTERPOLATION) {
    return 1.0 / len * f1xxCalcIntegrate(len, ratio);
  } else {
    // ADJUST ME
    const Int nlen = 512 + 1;
    const Int nratio = 128 + 1;
    static Interpolation2d fi(f1xxInterpPrime, 1.0, 0.0, nlen, PI, 0.0, nratio);
    const RealD lenp = len / (1.0 + len);
    const RealD psi = std::acos(ratio);
    return 1.0 / len * fi(lenp, psi);
  }
}

inline RealD f1xxCalc(const qlat::CoordinateD& c)
{
  const RealD len = coordinate_len(c);
  RealD ratio = c[3] / len;
  if (ratio > 1.0) {
    ratio = 1.0;
  } else if (ratio < -1.0) {
    ratio = -1.0;
  }
  const RealD ans = f1xxCalc(len, ratio);
  assert(false == std::isnan(ans));
  return ans;
}

struct Fg0xxIntegrand {
  RealD len, ratio;
  //
  Fg0xxIntegrand() { init(); }
  //
  void init()
  {
    len = 0.0;
    ratio = 0.0;
  }
  void init(const RealD len_, const RealD ratio_)
  {
    len = len_;
    ratio = ratio_;
    assert(len >= 0.0);
    assert(-1.0 <= ratio && ratio <= 1.0);
  }
  //
  RealD operator()(const RealD y) const
  {
    return std::exp(-y * ratio * len) * specialK0x(y * len);
  }
};

inline RealD g0xxCalcIntegrate(const RealD len, const RealD ratio)
{
  Fg0xxIntegrand g0xx;
  g0xx.init(len, ratio);
  return -1.0 / (8.0 * sqr(PI)) * integrate(g0xx, 0.0, 1.0, 0.0, 1.0e-12);
}

inline void recordG0xx(const RealD len, const RealD lenp, const RealD ratio,
                       const RealD g0xx)
{
  displayln(ssprintf("g0xx: %s %s %s %s", show(len).c_str(), show(lenp).c_str(),
                     show(ratio).c_str(), show(g0xx).c_str()));
}

inline RealD g0xxInterpPrime(const RealD lenp, const RealD psi)
{
  if (1.0 == lenp || 0.0 == lenp) {
    return 0.0;
  } else {
    const RealD len = lenp / (1.0 - lenp);
    const RealD ratio = std::cos(psi);
    const RealD g0xx = g0xxCalcIntegrate(len, ratio);
    // recordG0xx(len, lenp, ratio, g0xx);
    return g0xx;
  }
}

inline RealD g0xxCalc(const RealD len, const RealD ratio)
{
  if (0.0 == len) {
    return 10.0;  // ADJUST ME
  } else if (false == IS_USING_INTERPOLATION) {
    return 1.0 / len * g0xxCalcIntegrate(len, ratio);
  } else {
    // ADJUST ME
    const Int nlen = 512 + 1;
    const Int nratio = 128 + 1;
    static Interpolation2d fi(g0xxInterpPrime, 1.0, 0.0, nlen, PI, 0.0, nratio);
    const RealD lenp = len / (1.0 + len);
    const RealD psi = std::acos(ratio);
    return 1.0 / len * fi(lenp, psi);
  }
}

inline RealD g0xxCalc(const qlat::CoordinateD& c)
{
  const RealD len = coordinate_len(c);
  RealD ratio = c[3] / len;
  if (ratio > 1.0) {
    ratio = 1.0;
  } else if (ratio < -1.0) {
    ratio = -1.0;
  }
  const RealD ans = g0xxCalc(len, ratio);
  assert(false == std::isnan(ans));
  return ans;
}

inline void test_fderiv()
{
  TIMER_VERBOSE("test_fderiv");
  RngState rs("test_fderiv");
  const RealD low = -3.0;
  const RealD high = 3.0;
  const Long size = 4;
  for (Long i = 0; i < size; ++i) {
    CoordinateD c;
    for (Int mu = 0; mu < 4; ++mu) {
      c[mu] = u_rand_gen(rs, high, low);
    }
    const RealD len = coordinate_len(c);
    const RealD f1 = f1xxCalc(c);
    const RealD g0 = g0xxCalc(c);
    array<RealD, 4> fds;
    const RealD eps = 1.0e-5;
    for (Int mu = 0; mu < 4; ++mu) {
      CoordinateD cdp = c;
      CoordinateD cdm = c;
      cdp[mu] += eps;
      cdm[mu] -= eps;
      fds[mu] = (fxxCalc(cdp) - fxxCalc(cdm)) / (2 * eps);
    }
    array<RealD, 4> reffds;
    for (Int mu = 0; mu < 4; ++mu) {
      reffds[mu] = c[mu] / len * f1;
    }
    reffds[3] += g0;
    if (0 == qlat::get_id_node()) {
      for (Int mu = 0; mu < 4; ++mu) {
        displayln(ssprintf("%s : %d : %20e %20e %20e", show(c).c_str(), mu,
                           (fds[mu] - reffds[mu]) / reffds[mu], reffds[mu],
                           fds[mu]));
      }
    }
  }
}

inline void profile_gxxCalc()
{
  TIMER_VERBOSE_FLOPS("profile_gxxCalc");
  const Long size = 128;
  timer.flops += size * size;
  RealD sum = 0.0;
  for (Long i = 0; i < size; ++i) {
    for (Long j = 0; j < size; ++j) {
      const RealD len = (RealD)i / (RealD)size;
      const RealD ratio = (RealD)(2 * j) / (RealD)size - 1.0;
      sum += gxxCalc(len, ratio);
    }
  }
  DisplayInfo("", fname.c_str(), "%23.16e\n", sum);
}

inline void profile_f1xxCalc()
{
  TIMER_VERBOSE_FLOPS("profile_f1xxCalc");
  const Long size = 128;
  timer.flops += size * size;
  RealD sum = 0.0;
  for (Long i = 0; i < size; ++i) {
    for (Long j = 0; j < size; ++j) {
      const RealD len = (RealD)i / (RealD)size;
      const RealD ratio = (RealD)(2 * j) / (RealD)size - 1.0;
      sum += f1xxCalc(len, ratio);
      assert(false == std::isnan(sum));
    }
  }
  DisplayInfo("", fname.c_str(), "%23.16e\n", sum);
}

inline void profile_g0xxCalc()
{
  TIMER_VERBOSE_FLOPS("profile_g0xxCalc");
  const Long size = 128;
  timer.flops += size * size;
  RealD sum = 0.0;
  for (Long i = 0; i < size; ++i) {
    for (Long j = 0; j < size; ++j) {
      const RealD len = (RealD)i / (RealD)size;
      const RealD ratio = (RealD)(2 * j) / (RealD)size - 1.0;
      sum += g0xxCalc(len, ratio);
      assert(false == std::isnan(sum));
    }
  }
  DisplayInfo("", fname.c_str(), "%23.16e\n", sum);
}

inline void test_gxxCalc()
{
  TIMER_VERBOSE("test_gxxCalc");
  DisplayInfo("", fname.c_str(), "%23.16e\n", gxxCalc(1.0, -0.5));
  DisplayInfo("", fname.c_str(), "%23.16e\n", gxxCalc(0.0001, 1.0));
  DisplayInfo("", fname.c_str(), "%23.16e\n", gxxCalc(1.0, -1.0));
  DisplayInfo("", fname.c_str(), "%23.16e\n", gxxCalc(10.0, 1.0));
  DisplayInfo("", fname.c_str(), "%23.16e\n", gxxCalc(10.0, -1.0));
  profile_gxxCalc();
}

inline void test_hxCalc()
{
  TIMER_VERBOSE("test_hxCalc");
  DisplayInfo("", fname.c_str(), "%23.16e\n", hxCalc(0.0));
  DisplayInfo("", fname.c_str(), "%23.16e\n", hxCalc(0.0001));
  DisplayInfo("", fname.c_str(), "%23.16e\n", hxCalc(1.0));
  DisplayInfo("", fname.c_str(), "%23.16e\n", hxCalc(10.0));
}

inline void test_g0xxCalc()
{
  TIMER_VERBOSE("test_g0xxCalc");
  DisplayInfo("", fname.c_str(), "%23.16e\n", g0xxCalc(1.0, -0.5));
  DisplayInfo("", fname.c_str(), "%23.16e\n", g0xxCalc(0.0001, 1.0));
  DisplayInfo("", fname.c_str(), "%23.16e\n", g0xxCalc(1.0, -1.0));
  DisplayInfo("", fname.c_str(), "%23.16e\n", g0xxCalc(10.0, 1.0));
  DisplayInfo("", fname.c_str(), "%23.16e\n", g0xxCalc(10.0, -1.0));
  profile_g0xxCalc();
}

inline void test_f1xxCalc()
{
  TIMER_VERBOSE("test_f1xxCalc");
  DisplayInfo("", fname.c_str(), "%23.16e\n", f1xxCalc(1.0, -0.5));
  DisplayInfo("", fname.c_str(), "%23.16e\n", f1xxCalc(0.0001, 1.0));
  DisplayInfo("", fname.c_str(), "%23.16e\n", f1xxCalc(1.0, -1.0));
  DisplayInfo("", fname.c_str(), "%23.16e\n", f1xxCalc(10.0, 1.0));
  DisplayInfo("", fname.c_str(), "%23.16e\n", f1xxCalc(10.0, -1.0));
  profile_f1xxCalc();
}

inline void test_fCalc()
{
  TIMER_VERBOSE("test_fCalc");
  test_gxxCalc();
  test_hxCalc();
  test_f1xxCalc();
  test_g0xxCalc();
  test_fderiv();
}

}  // namespace qlat
