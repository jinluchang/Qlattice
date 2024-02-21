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
inline int cubaFunction(const int* ndim, const double x[], const int* ncomp,
                        double f[], void* userdata)
{
  // TIMER("cubaFunction");
  std::vector<double> vx(*ndim, 0.0);
  std::memcpy(vx.data(), x, *ndim * sizeof(double));
  std::vector<double> vf((*((const F*)userdata))(vx));
  assert(vf.size() == (size_t)*ncomp);
  std::memcpy(f, vf.data(), *ncomp * sizeof(double));
  return 0;
}

template <class F>
void integrateDivonne(std::vector<double>& integral, std::vector<double>& error,
                      std::vector<double>& prob, int& nregions, long long int& neval,
                      int& fail, const int ndim, const int ncomp, const F& f,
                      const double epsabs = 0.0, const double epsrel = 1.0e-5,
                      const int flags = 0, const int seed = 23,
                      const long mineval = 128,
                      const long maxeval = 16 * 1024 * 1024 * 4,
                      const int key1 = 7, const int key2 = 7,
                      const int key3 = 1, const int maxpass = 5,
                      const double border = 0.0, const double maxchisq = 10.0,
                      const double mindeviation = 0.25, const long long int ngiven = 0,
                      int ldxgiven = 0, double xgiven[] = NULL,
                      const long long int nextra = 0, peakfinder_t peakfinder = NULL,
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
  // neval, fail); for (int i = 0; i < ncomp; ++i) {
  //   DisplayInfo("", fname, "i=%d integral=%23.16e ; error=%23.16e ;
  //   fail=%f\n", i, integral[i], error[i], prob[i]);
  // }
#else
  qerr("integrateDivonne: QLAT_NO_CUBA defined. (Cuba library not found)");
#endif
}

template <class F>
void integrateCuhre(std::vector<double>& integral, std::vector<double>& error,
                    std::vector<double>& prob, int& nregions, long long int& neval,
                    int& fail, const int ndim, const int ncomp, const F& f,
                    const double epsabs = 0.0, const double epsrel = 1.0e-5,
                    const int flags = 0, const long mineval = 128,
                    const long maxeval = 16 * 1024 * 1024 * 4, const int key = 7,
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
  // neval, fail); for (int i = 0; i < ncomp; ++i) {
  //   DisplayInfo("", fname, "i=%d integral=%23.16e ; error=%23.16e ;
  //   fail=%f\n", i, integral[i], error[i], prob[i]);
  // }
#else
  qerr("integrateDivonne: QLAT_NO_CUBA defined. (Cuba library not found)");
#endif
}

inline std::vector<double> test_integrand4d(const std::vector<double>& vx)
{
  // TIMER_VERBOSE("test_integrand4d");
  assert(4 == vx.size());
  std::vector<double> ans(1);
  ans[0] = vx[0] * vx[1] * vx[2] * vx[3] + sin(vx[0] * PI) * sin(vx[1] * PI) +
           sqrt(vx[3]) * sqrt(vx[2]) * sqrt(vx[1]);
  // DisplayInfo("", fname, "%f %f %f %f %f\n", ans[0], vx[0], vx[1], vx[2],
  // vx[3]);
  return ans;
}

inline void test_integrationMultidimensional()
{
  TIMER_VERBOSE("test_integrationMultidimensional");
  std::vector<double> integral, error, prob;
  int nregions, fail;
  long long int neval;
  integrateCuhre(integral, error, prob, nregions, neval, fail, 4, 1,
                 test_integrand4d);
  fdisplayln(stdout, ssprintf("%f %f %f", integral[0], error[0], prob[0]));
  integrateDivonne(integral, error, prob, nregions, neval, fail, 4, 1,
                   test_integrand4d);
  fdisplayln(stdout, ssprintf("%f %f %f", integral[0], error[0], prob[0]));
}

}  // namespace qlat
