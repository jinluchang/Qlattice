#pragma once

#include <qlat/qlat.h>

#ifndef QLAT_NO_CUBA
#include <cuba.h>
#else
using peakfinder_t = void*;
#endif

#include <vector>
#include <cstring>
#include <cassert>

namespace qlat
{

inline bool has_cuba()
{
#ifndef QLAT_NO_CUBA
  return true;
#else
  return false;
#endif
}

template <class F>
inline Int cubaFunction(const Int* ndim, const RealD x[], const Int* ncomp,
                        RealD f[], void* userdata)
{
  // TIMER("cubaFunction");
  std::vector<RealD> vx(*ndim, 0.0);
  std::memcpy(vx.data(), x, *ndim * sizeof(RealD));
  std::vector<RealD> vf((*((const F*)userdata))(vx));
  assert(vf.size() == (size_t)*ncomp);
  std::memcpy(f, vf.data(), *ncomp * sizeof(RealD));
  return 0;
}

template <class F>
void integrateDivonne(
    std::vector<RealD>& integral, std::vector<RealD>& error,
    std::vector<RealD>& prob, Int& nregions, long long& neval, Int& fail,
    const Int ndim, const Int ncomp, const F& f, const RealD epsabs = 0.0,
    const RealD epsrel = 1.0e-5, const Int flags = 0, const Int seed = 23,
    const long mineval = 128, const long maxeval = 16 * 1024 * 1024 * 4,
    const Int key1 = 7, const Int key2 = 7, const Int key3 = 1,
    const Int maxpass = 5, const RealD border = 0.0,
    const RealD maxchisq = 10.0, const RealD mindeviation = 0.25,
    const long long ngiven = 0, Int ldxgiven = 0, RealD xgiven[] = NULL,
    const long long nextra = 0, peakfinder_t peakfinder = NULL,
    const char* statefile = NULL, void* spin = NULL)
{
  TIMER("integrateDivonne");
#ifndef QLAT_NO_CUBA
  if (0 == ldxgiven) {
    ldxgiven = ndim;
  }
  integral.resize(ncomp);
  error.resize(ncomp);
  prob.resize(ncomp);
  // cubacores(0, 0);
  llDivonne(ndim, ncomp, cubaFunction<F>, (void*)&f, 1, epsrel, epsabs, flags,
          seed, mineval, maxeval, key1, key2, key3, maxpass, border, maxchisq,
          mindeviation, ngiven, ldxgiven, xgiven, nextra, peakfinder, statefile,
          spin, &nregions, &neval, &fail, integral.data(), error.data(),
          prob.data());
  // DisplayInfo("", fname, "nregions=%d ; neval=%d ; fail=%d\n", nregions,
  // neval, fail); for (Int i = 0; i < ncomp; ++i) {
  //   DisplayInfo("", fname, "i=%d integral=%23.16e ; error=%23.16e ;
  //   fail=%f\n", i, integral[i], error[i], prob[i]);
  // }
#else
  (void)integral;
  (void)error;
  (void)prob;
  (void)nregions;
  (void)neval;
  (void)fail;
  (void)ndim;
  (void)ncomp;
  (void)f;
  (void)epsabs;
  (void)epsrel;
  (void)flags;
  (void)seed;
  (void)mineval;
  (void)maxeval;
  (void)key1;
  (void)key2;
  (void)key3;
  (void)maxpass;
  (void)border;
  (void)maxchisq;
  (void)mindeviation;
  (void)ngiven;
  (void)ldxgiven;
  (void)xgiven;
  (void)nextra;
  (void)peakfinder;
  (void)statefile;
  (void)spin;
  qerr("integrateDivonne: QLAT_NO_CUBA defined. (Cuba library not found)");
#endif
}

template <class F>
void integrateCuhre(std::vector<RealD>& integral, std::vector<RealD>& error,
                    std::vector<RealD>& prob, Int& nregions, long long& neval,
                    Int& fail, const Int ndim, const Int ncomp, const F& f,
                    const RealD epsabs = 0.0, const RealD epsrel = 1.0e-5,
                    const Int flags = 0, const long mineval = 128,
                    const long maxeval = 16 * 1024 * 1024 * 4, const Int key = 7,
                    const char* statefile = NULL, void* spin = NULL)
{
  TIMER("integrateCuhre");
#ifndef QLAT_NO_CUBA
  integral.resize(ncomp);
  error.resize(ncomp);
  prob.resize(ncomp);
  // cubacores(0, 0);
  llCuhre(ndim, ncomp, cubaFunction<F>, (void*)&f, 1, epsrel, epsabs, flags,
        mineval, maxeval, key, statefile, spin, &nregions, &neval, &fail,
        integral.data(), error.data(), prob.data());
  // DisplayInfo("", fname, "nregions=%d ; neval=%d ; fail=%d\n", nregions,
  // neval, fail); for (Int i = 0; i < ncomp; ++i) {
  //   DisplayInfo("", fname, "i=%d integral=%23.16e ; error=%23.16e ;
  //   fail=%f\n", i, integral[i], error[i], prob[i]);
  // }
#else
  (void)integral;
  (void)error;
  (void)prob;
  (void)nregions;
  (void)neval;
  (void)fail;
  (void)ndim;
  (void)ncomp;
  (void)f;
  (void)epsabs;
  (void)epsrel;
  (void)flags;
  (void)mineval;
  (void)maxeval;
  (void)key;
  (void)statefile;
  (void)spin;
  qerr("integrateDivonne: QLAT_NO_CUBA defined. (Cuba library not found)");
#endif
}

inline std::vector<RealD> test_integrand4d(const std::vector<RealD>& vx)
{
  // TIMER_VERBOSE("test_integrand4d");
  assert(4 == vx.size());
  std::vector<RealD> ans(1);
  ans[0] = vx[0] * vx[1] * vx[2] * vx[3] + sin(vx[0] * PI) * sin(vx[1] * PI) +
           sqrt(vx[3]) * sqrt(vx[2]) * sqrt(vx[1]);
  // DisplayInfo("", fname, "%f %f %f %f %f\n", ans[0], vx[0], vx[1], vx[2],
  // vx[3]);
  return ans;
}

inline void test_integrationMultidimensional()
{
  TIMER_VERBOSE("test_integrationMultidimensional");
  std::vector<RealD> integral, error, prob;
  Int nregions, fail;
  long long neval;
  integrateCuhre(integral, error, prob, nregions, neval, fail, 4, 1,
                 test_integrand4d);
  fdisplayln(stdout, ssprintf("%f %f %f", integral[0], error[0], prob[0]));
  integrateDivonne(integral, error, prob, nregions, neval, fail, 4, 1,
                   test_integrand4d);
  fdisplayln(stdout, ssprintf("%f %f %f", integral[0], error[0], prob[0]));
}

}  // namespace qlat
