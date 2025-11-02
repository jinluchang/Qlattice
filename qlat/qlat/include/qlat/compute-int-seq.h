#pragma once

#include "muon-line-config.h"
#include "utils.h"
#include "root-fsolver.h"
#include "interpolation.h"
#include "compute-f.h"

#include <qlat/qlat.h>

#include <cmath>
#include <cstdlib>
#include <cassert>

namespace qlat
{  //

inline RealD vFromPsi(const RealD psi)
{
  return psi - std::sin(psi) * std::cos(psi);
}

struct EquationPsiFromV {
  RealD v;
  //
  RealD operator()(const RealD psi) const { return vFromPsi(psi) - v; }
};

inline RealD psiFromVSolve(const RealD v)
{
  if (v == 0 || v == PI) {
    return v;
  }
  EquationPsiFromV f;
  f.v = v;
  return rootFSolver(f, v, 0.0, PI, 1.0e-12, 0.0);
}

inline RealD psiFromV(const RealD v)
{
  if (false == IS_USING_INTERPOLATION) {
    return psiFromVSolve(v);
  } else {
    // ADJUST ME
    const Int nv = 16 * 1024;
    static Interpolation fi(psiFromVSolve, PI, 0.0, nv);
    return fi(v);
  }
}

inline RealD uFromTheta(const RealD theta) { return std::cos(theta); }

inline RealD thetaFromU(const RealD u) { return std::acos(u); }

inline RealD rpFromR(const RealD r, const RealD r0)
{
  const RealD ra = r / r0;
  return ra / (1.0 + ra);
}

inline RealD regularizeRp(const RealD rp)
{
  assert(0.0 <= rp && rp <= 1.0);
  // ADJUST ME
  const RealD eps = 1.0e-5;
  if (1.0 == rp) {
    return 1.0 - eps;
  } else {
    return rp;
  }
}

inline RealD rFromRp(const RealD rp, const RealD r0)
{
  const RealD rrp = regularizeRp(rp);
  return rrp / (1.0 - rrp) * r0;
}

inline qlat::CoordinateD etaFromParam(const qlat::CoordinateD& param,
                                      const RealD r0)
{
  // TIMER_VERBOSE("etaFromParam");
  const RealD rp = param[0];
  const RealD v = param[1];
  const RealD u = param[2];
  const RealD phi = param[3];
  assert(0.0 <= rp && rp <= 1.0);
  assert(0.0 <= v && v <= PI);
  assert(-1.0 <= u && u <= 1.0);
  assert(-PI <= phi && phi <= PI);
  const RealD r = rFromRp(rp, r0);
  const RealD psi = psiFromV(v);
  const RealD theta = thetaFromU(u);
  const RealD t = r * std::cos(psi);
  const RealD rSinPsi = r * std::sin(psi);
  const RealD z = rSinPsi * std::cos(theta);
  const RealD rSinPsiSinTheta = rSinPsi * std::sin(theta);
  const RealD x = rSinPsiSinTheta * std::cos(phi);
  const RealD y = rSinPsiSinTheta * std::sin(phi);
  CoordinateD eta(x, y, z, t);
  return eta;
}

inline qlat::CoordinateD paramFromEta(const qlat::CoordinateD& eta,
                                      const RealD r0)
{
  // TIMER_VERBOSE("paramFromEta");
  const RealD x = eta[0];
  const RealD y = eta[1];
  const RealD z = eta[2];
  const RealD t = eta[3];
  const RealD r = std::sqrt(sqr(x) + sqr(y) + sqr(z) + sqr(t));
  const RealD rp = rpFromR(r, r0);
  if (0.0 == r) {
    return CoordinateD(rp, PI / 2, 0.0, 0.0);
  } else {
    const RealD psi = std::acos(t / r);
    const RealD v = vFromPsi(psi);
    const RealD rspatial = std::sqrt(sqr(x) + sqr(y) + sqr(z));
    if (0.0 == rspatial) {
      return CoordinateD(rp, v, 0.0, 0.0);
    } else {
      const RealD u = z / rspatial;
      if (0.0 == x && 0.0 == y) {
        return CoordinateD(rp, v, u, 0.0);
      } else {
        const RealD phi = std::atan2(y, x);
        CoordinateD param(rp, v, u, phi);
        return param;
      }
    }
  }
}

struct MuonLineIntegrand {
  qlat::CoordinateD x, y;
  RealD r0;
  //
  MuonLineIntegrand() { init(); }
  //
  void init()
  {
    for (Int i = 0; i < 4; ++i) {
      x[i] = 0.0;
      y[i] = 0.0;
    }
    r0 = (coordinate_len(x) + coordinate_len(y)) * 0.5 + 0.001;
  }
  void init(const qlat::CoordinateD& x_, const qlat::CoordinateD& y_)
  {
    x = x_;
    y = y_;
    r0 = (coordinate_len(x) + coordinate_len(y)) * 0.5 + 0.001;
  }
  //
  std::vector<RealD> integrand(const RealD rp,
                                const qlat::CoordinateD& eta) const
  {
    // TIMER_VERBOSE("MuonLineIntegrand");
    std::vector<RealD> ans(25, 0.0);
    const RealD leta = coordinate_len(eta);
    const CoordinateD c1 = eta - y;
    const CoordinateD c2 = x - eta;
    const RealD lc1 = coordinate_len(c1);
    const RealD lc2 = coordinate_len(c2);
    const CoordinateD s1 = -eta + y;
    const CoordinateD s2 = -x + eta;
    const RealD ls1 = coordinate_len(s1);
    const RealD ls2 = coordinate_len(s2);
    if (leta < 1e-10 || lc1 < 1e-10 || lc2 < 1e-10 || ls1 < 1e-10 ||
        ls2 < 1e-10) {
      return ans;
    }
    const RealD fxxC1 = fxxCalc(c1);
    const RealD fxxC2 = fxxCalc(c2);
    const RealD f1xxC1 = f1xxCalc(c1);
    const RealD f1xxC2 = f1xxCalc(c2);
    const RealD f0xxC1 = c1[3] / lc1 * f1xxC1 + g0xxCalc(c1);
    const RealD f0xxC2 = c2[3] / lc2 * f1xxC2 + g0xxCalc(c2);
    const RealD fxxS1 = fxxCalc(s1);
    const RealD fxxS2 = fxxCalc(s2);
    const RealD f1xxS1 = f1xxCalc(s1);
    const RealD f1xxS2 = f1xxCalc(s2);
    const RealD f0xxS1 = s1[3] / ls1 * f1xxS1 + g0xxCalc(s1);
    const RealD f0xxS2 = s2[3] / ls2 * f1xxS2 + g0xxCalc(s2);
    const RealD fxx1sub = fxxCalc(eta);
    const RealD fxx2sub = fxxCalc(-eta);
    const RealD f1xx1sub = f1xxCalc(eta);
    const RealD f1xx2sub = f1xxCalc(-eta);
    const RealD f0xx1sub = eta[3] / leta * f1xx1sub + g0xxCalc(eta);
    const RealD f0xx2sub = -eta[3] / leta * f1xx2sub + g0xxCalc(-eta);
    array<RealD, 5> fc1xs;
    fc1xs[0] = c1[0] / lc1 * f1xxC1;
    fc1xs[1] = c1[1] / lc1 * f1xxC1;
    fc1xs[2] = c1[2] / lc1 * f1xxC1;
    fc1xs[3] = f0xxC1;
    fc1xs[4] = fxxC1;
    array<RealD, 5> fc2xs;
    fc2xs[0] = c2[0] / lc2 * f1xxC2;
    fc2xs[1] = c2[1] / lc2 * f1xxC2;
    fc2xs[2] = c2[2] / lc2 * f1xxC2;
    fc2xs[3] = f0xxC2;
    fc2xs[4] = fxxC2;
    array<RealD, 5> fs1xs;
    fs1xs[0] = s1[0] / ls1 * f1xxS1;
    fs1xs[1] = s1[1] / ls1 * f1xxS1;
    fs1xs[2] = s1[2] / ls1 * f1xxS1;
    fs1xs[3] = f0xxS1;
    fs1xs[4] = fxxS1;
    array<RealD, 5> fs2xs;
    fs2xs[0] = s2[0] / ls2 * f1xxS2;
    fs2xs[1] = s2[1] / ls2 * f1xxS2;
    fs2xs[2] = s2[2] / ls2 * f1xxS2;
    fs2xs[3] = f0xxS2;
    fs2xs[4] = fxxS2;
    array<RealD, 5> f1subxs;
    f1subxs[0] = eta[0] / leta * f1xx1sub;
    f1subxs[1] = eta[1] / leta * f1xx1sub;
    f1subxs[2] = eta[2] / leta * f1xx1sub;
    f1subxs[3] = f0xx1sub;
    f1subxs[4] = fxx1sub;
    array<RealD, 5> f2subxs;
    f2subxs[0] = -eta[0] / leta * f1xx2sub;
    f2subxs[1] = -eta[1] / leta * f1xx2sub;
    f2subxs[2] = -eta[2] / leta * f1xx2sub;
    f2subxs[3] = f0xx2sub;
    f2subxs[4] = fxx2sub;
    if (IS_USING_AGGRESSIVE_SUBTRACTION) {
      for (Int i = 0; i < 5; ++i) {
        if (IS_USING_MORE_AGGRESSIVE_SUBTRACTION || 4 == i) {
          if (0 == IS_USING_ONLY_SUBTRACTION) {
            fc1xs[i] -= f1subxs[i];
            fc2xs[i] -= f2subxs[i];
            fs1xs[i] -= f2subxs[i];
            fs2xs[i] -= f1subxs[i];
          } else if (1 == IS_USING_ONLY_SUBTRACTION) {
            fc1xs[i] = -f1subxs[i];
            fs2xs[i] = -f1subxs[i];
          } else if (2 == IS_USING_ONLY_SUBTRACTION) {
            fc2xs[i] = -f2subxs[i];
            fs1xs[i] = -f2subxs[i];
          }
        }
      }
    }
    const RealD factor = sqr(r0) / sqr(1.0 - rp) * rp / (1.0 - rp);
    for (Int i = 0; i < 5; ++i) {
      for (Int j = 0; j < 5; ++j) {
        ans[i * 5 + j] =
            0.5 * factor * (fc1xs[i] * fc2xs[j] - fs1xs[i] * fs2xs[j]);
        // the first 0.5 because we need to address the RealD counting in the
        // subtraction term
      }
    }
    if (0.0 == x[2] && 0.0 == y[2]) {
      if (0.0 == x[1] || 0.0 == y[1] || x[1] == y[1]) {
        // We are computing muon-line for the x,y generated from the params
        ans[2] = 0.0;
        ans[7] = 0.0;
        ans[10] = 0.0;
        ans[11] = 0.0;
        ans[13] = 0.0;
        ans[14] = 0.0;
        ans[17] = 0.0;
        ans[22] = 0.0;
      }
    }
    // the last component is not relevant
    ans[24] = 0.0;
    for (size_t i = 0; i < ans.size(); ++i) {
      if (std::isnan(ans[i])) {
        displayln(ssprintf("MuonLineIntegrand: x   =%s", show(x).c_str()));
        displayln(ssprintf("MuonLineIntegrand: y   =%s", show(y).c_str()));
        displayln(ssprintf("MuonLineIntegrand: eta =%s", show(eta).c_str()));
        for (size_t k = 0; k < ans.size(); ++k) {
          displayln(ssprintf("MuonLineIntegrand: ans[%d]=%s", k,
                             show(ans[k]).c_str()));
        }
        assert(false);
      }
    }
    return ans;
  }
  //
  std::vector<RealD> integrand_eta_reflection_avg(
      const RealD rp, const qlat::CoordinateD& eta) const
  {
    const std::vector<RealD> pos = integrand(rp, eta);
    const std::vector<RealD> neg = integrand(rp, -eta);
    qassert(25 == pos.size());
    qassert(25 == neg.size());
    std::vector<RealD> ret(25);
    for (size_t i = 0; i < ret.size(); ++i) {
      ret[i] = 1.0 / 2.0 * (pos[i] + neg[i]);
    }
    return ret;
  }
  //
  std::vector<RealD> operator()(const qlat::CoordinateD& param) const
  {
    const RealD rp = regularizeRp(param[0]);
    const CoordinateD eta = etaFromParam(param, r0);
    // ADJUST ME
    // return integrand(rp, eta);
    return integrand_eta_reflection_avg(rp, eta);
  }
};

struct MuonLineIntegrandPhi {
  MuonLineIntegrand mli;
  RealD rp, v, u;
  //
  RealD operator()(const RealD phi) const
  {
    TIMER("MuonLineIntegrandPhi");
    const CoordinateD param(rp, v, u, phi);
    return mli(param)[0];
  }
};

struct MuonLineIntegrandU {
  MuonLineIntegrand mli;
  RealD rp, v;
  //
  RealD operator()(const RealD u) const
  {
    TIMER("MuonLineIntegrandU");
    MuonLineIntegrandPhi mliPhi;
    mliPhi.mli = mli;
    mliPhi.rp = rp;
    mliPhi.v = v;
    mliPhi.u = u;
    return integrate(mliPhi, -PI, PI, 1.0e-10, 0.0);
  }
};

struct MuonLineIntegrandV {
  MuonLineIntegrand mli;
  RealD rp;
  //
  RealD operator()(const RealD v) const
  {
    TIMER("MuonLineIntegrandV");
    MuonLineIntegrandU mliU;
    mliU.mli = mli;
    mliU.rp = rp;
    mliU.v = v;
    return integrate(mliU, -1.0, 1.0, 1.0e-9, 0.0);
  }
};

struct MuonLineIntegrandRp {
  MuonLineIntegrand mli;
  //
  RealD operator()(const RealD rp) const
  {
    Timer::display();
    TIMER_VERBOSE("MuonLineIntegrandRp");
    MuonLineIntegrandV mliV;
    mliV.mli = mli;
    mliV.rp = rp;
    return integrate(mliV, 0.0, PI, 1.0e-8, 0.0);
  }
};

inline RealD integrateMuonLineTest(const MuonLineIntegrand& mli)
{
  TIMER_VERBOSE("integrateMuonLineTest");
  MuonLineIntegrandRp mliRp;
  mliRp.mli = mli;
  return integrate(mliRp, 0.0, 1.0, 1e9, 0.0);
}

inline void profile_computeIntSeq_integrand()
{
  TIMER_VERBOSE_FLOPS("profile_computeIntSeq_integrand");
  RngState rs("profile_computeIntSeq_integrand");
  const Int size = 128;
  const RealD low = -3.0;
  const RealD high = 3.0;
  RealD sum = 0.0;
  for (Int i = 0; i < size; ++i) {
    CoordinateD x, y, eta;
    for (Int m = 0; m < 4; ++m) {
      x[m] = u_rand_gen(rs, high, low);
      y[m] = u_rand_gen(rs, high, low);
      eta[m] = u_rand_gen(rs, high, low);
    }
    MuonLineIntegrand mli;
    mli.init(x, y);
    sum += mli(paramFromEta(eta, mli.r0))[0];
  }
  DisplayInfo("", fname.c_str(), "sum=%23.16e\n", sum);
  timer.flops += size;
}

inline void profile_computeIntSeq()
{
  TIMER_VERBOSE("profile_computeIntSeq");
  profile_computeIntSeq_integrand();
}

inline void test_computeIntSeq()
{
  TIMER_VERBOSE("test_computeIntSeq");
  const RealD v = 2.0;
  const RealD psi = psiFromV(v);
  DisplayInfo("", fname.c_str(), "v=%23.16e ; psi=%23.16e -> v=%23.16e\n", v, psi,
              vFromPsi(psi));
  const CoordinateD eta(0.2, 0.2, 0.4, 1.4);
  const CoordinateD param(paramFromEta(eta, 1.1));
  DisplayInfo("", fname.c_str(), "eta   =%s\n", show(eta).c_str());
  DisplayInfo("", fname.c_str(), "param =%s\n", show(param).c_str());
  DisplayInfo("", fname.c_str(), "eta   =%s\n", show(etaFromParam(param, 1.1)).c_str());
  const CoordinateD x(0.1, 0.2, 0.0, 0.5);
  const CoordinateD y(0.3, 0.0, -0.2, 0.1);
  MuonLineIntegrand mli;
  mli.init(x, y);
  DisplayInfo("", fname.c_str(), "x =%s\n", show(x).c_str());
  DisplayInfo("", fname.c_str(), "y =%s\n", show(y).c_str());
  DisplayInfo("", fname.c_str(), "MuonLineIntegrand=%23.16e\n", mli(param)[0]);
  {
    MuonLineIntegrandPhi mliPhi;
    mliPhi.mli.init(x, y);
    const CoordinateD param(paramFromEta(eta, mli.r0));
    mliPhi.rp = param[0];
    mliPhi.v = param[1];
    mliPhi.u = param[2];
    const RealD mliPhiInt = integrate(mliPhi, -PI, PI);
    DisplayInfo("", fname.c_str(), "mliPhiInt=%23.16e\n", mliPhiInt);
  }
  // DisplayInfo("", fname.c_str(), "mliInt=%23.16e\n", integrateMuonLineTest(mli));
  profile_computeIntSeq();
}

}  // namespace qlat
