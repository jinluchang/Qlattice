#pragma once

#include <qlat/field.h>
#include <qlat/field-fft.h>

#include <Eigen/Dense>

#include <cmath>

QLAT_START_NAMESPACE

struct QedGaugeField : FieldM<Complex,DIM>
{
  virtual const char* cname()
  {
    return "QedGaugeField";
  }
};

struct ComplexScalerField : FieldM<Complex,1>
{
  virtual const char* cname()
  {
    return "ComplexScalerField";
  }
};

using SpinMatrix = Eigen::Matrix<Complex,4,4>;

inline void setZero(SpinMatrix& sm)
{
  sm.setZero();
}

inline void setUnit(SpinMatrix& sm, const Complex& coef = 1.0)
{
  setZero(sm);
  for (int i = 0; i < sm.rows() && i < sm.cols(); ++i) {
    sm(i,i) = coef;
  }
}

inline double norm(const SpinMatrix& sm)
{
  return sm.squaredNorm();
}

inline std::string show(const SpinMatrix& sm)
{
  std::ostringstream out;
  out << sm;
  return out.str();
}

struct SpinMatrixConstants
{
  SpinMatrix unit;
  std::array<SpinMatrix,4> gammas;
  // Not using CPS's convention, but a more standard one.
  SpinMatrix gamma5;
  // Same as CPS's gamma5
  std::array<SpinMatrix,3> capSigmas;
  //
  SpinMatrixConstants()
  {
    init();
  }
  //
  void init()
  {
    unit <<
        1,   0,   0,   0,
        0,   1,   0,   0,
        0,   0,   1,   0,
        0,   0,   0,   1;
    // gamma_x
    gammas[0] <<
        0,   0,   0,   1,
        0,   0,   1,   0,
        0,  -1,   0,   0,
       -1,   0,   0,   0;
    // gamma_y
    gammas[1] <<
        0,   0,   0, -ii,
        0,   0,  ii,   0,
        0,  ii,   0,   0,
      -ii,   0,   0,   0;
    // gamma_z
    gammas[2] <<
        0,   0,   1,   0,
        0,   0,   0,  -1,
       -1,   0,   0,   0,
        0,   1,   0,   0;
    gammas[0] *= -ii;
    gammas[1] *= -ii;
    gammas[2] *= -ii;
    // gamma_t
    gammas[3] <<
        0,   0,   1,   0,
        0,   0,   0,   1,
        1,   0,   0,   0,
        0,   1,   0,   0;
    // gamma_5
    gamma5 <<
        1,   0,   0,   0,
        0,   1,   0,   0,
        0,   0,  -1,   0,
        0,   0,   0,  -1;
    // Sigma_x
    capSigmas[0] <<
        0,   1,   0,   0,
        1,   0,   0,   0,
        0,   0,   0,   1,
        0,   0,   1,   0;
    // Sigma_y
    capSigmas[1] <<
        0, -ii,   0,   0,
       ii,   0,   0,   0,
        0,   0,   0, -ii,
        0,   0,  ii,   0;
    // Sigma_z
    capSigmas[2] <<
        1,   0,   0,   0,
        0,  -1,   0,   0,
        0,   0,   1,   0,
        0,   0,   0,  -1;
  }
  //
  static const SpinMatrixConstants& getInstance()
  {
    static SpinMatrixConstants smcs;
    return smcs;
  }
  //
  static const SpinMatrix& getUnit()
  {
    return getInstance().unit;
  }
  static const SpinMatrix& getGamma(int mu)
  {
    assert(0 <= mu && mu < 4);
    return getInstance().gammas[mu];
  }
  static const SpinMatrix& getGamma5()
  {
    return getInstance().gamma5;
  }
  static const SpinMatrix& getCapSigma(int i)
  {
    return getInstance().capSigmas[i];
  }
};

struct SpinPropagator4d : FieldM<SpinMatrix,1>
{
  virtual const char* cname()
  {
    return "SpinPropagator4d";
  }
};

inline void propMomPhotonInvert(QedGaugeField& egf, const std::array<double,DIM>& momtwist)
  // Feynman Gauge
  // All spatial zero mode removed.
  // egf in momentum space.
{
  TIMER("propMomPhotonInvert");
  const Geometry& geo = egf.geo;
  for (long index = 0; index < geo.localVolume(); ++index) {
    Coordinate kl; geo.coordinateFromIndex(kl, index);
    Coordinate kg; geo.coordinateGfL(kg, kl);
    std::array<double,DIM> kk;
    std::array<double,DIM> ks;
    double s2 = 0.0;
    for (int i = 0; i < DIM; i++) {
      kg[i] = signMod(kg[i], geo.totalSite(i));
      kk[i] = 2.0 * PI * (kg[i] + momtwist[i]) / (double)geo.totalSite(i);
      s2 += 4.0 * sqr(std::sin(kk[i] / 2.0));
    }
    if (0.0 == kk[0] && 0.0 == kk[1] && 0.0 == kk[2]) {
      for (int mu = 0; mu < geo.multiplicity; ++mu) {
        egf.getElem(kl, mu) *= 0.0;
      }
    } else {
      double s2inv = 1.0 / s2;
      for (int mu = 0; mu < geo.multiplicity; ++mu) {
        egf.getElem(kl, mu) *= s2inv;
      }
    }
  }
}

inline void propPhotonInvert(QedGaugeField& egf, const std::array<double,DIM>& momtwist)
  // Feynman Gauge
  // All spatial zero mode removed.
  // egf in coordinate space.
{
  TIMER_VERBOSE("propPhotonInvert");
  const Geometry& geo = egf.geo;
  fftComplexField(egf, true);
  propMomPhotonInvert(egf, momtwist);
  fftComplexField(egf, false);
  egf *= 1.0 / geo.totalVolume();
}

inline void propMomComplexScalerInvert(ComplexScalerField& csf, const double mass, const std::array<double,DIM>& momtwist)
{
  TIMER("propMomComplexScalerInvert");
  const Geometry& geo = csf.geo;
  for (long index = 0; index < geo.localVolume(); ++index) {
    Coordinate kl; geo.coordinateFromIndex(kl, index);
    Coordinate kg; geo.coordinateGfL(kg, kl);
    std::array<double,DIM> kk;
    std::array<double,DIM> ks;
    double s2 = 0.0;
    for (int i = 0; i < DIM; i++) {
      kg[i] = signMod(kg[i], geo.totalSite(i));
      kk[i] = 2.0 * PI * (kg[i] + momtwist[i]) / (double)geo.totalSite(i);
      s2 += 4.0 * sqr(std::sin(kk[i] / 2.0));
    }
    double s2inv = 1.0 / s2;
    for (int mu = 0; mu < geo.multiplicity; ++mu) {
      csf.getElem(kl, mu) *= s2inv;
    }
  }
}

inline void propComplexScalerInvert(ComplexScalerField& csf, const double mass, const std::array<double,DIM>& momtwist)
{
  TIMER_VERBOSE("propComplexScalerInvert");
  const Geometry& geo = csf.geo;
  fftComplexField(csf, true);
  propMomComplexScalerInvert(csf, mass, momtwist);
  fftComplexField(csf, false);
  csf *= 1.0 / geo.totalVolume();
}

inline void propMomSpinPropagator4d(SpinPropagator4d& sp4d, const double mass, const std::array<double,DIM>& momtwist)
  // DWF infinite L_s
  // M_5 = 1.0
{
  TIMER("propMomSpinPropagator4d");
  const Geometry& geo = sp4d.geo;
  const double m5 = 1.0;
#pragma omp parallel for
  for (long index = 0; index < geo.localVolume(); ++index) {
    Coordinate kl; geo.coordinateFromIndex(kl, index);
    Coordinate kg; geo.coordinateGfL(kg, kl);
    std::array<double,DIM> kk, ks;
    double p2 = 0.0;
    double wp = 1.0 - m5;
    SpinMatrix pg; setZero(pg);
    for (int i = 0; i < DIM; ++i) {
      kg[i] = signMod(kg[i], geo.totalSite(i));
      kk[i] = 2.0 * PI * (kg[i] + momtwist[i]) / (double)geo.totalSite(i);
      ks[i] = sin(kk[i]);
      pg += SpinMatrixConstants::getGamma(i) * (Complex)ks[i];
      p2 += sqr(ks[i]);
      wp += 2.0 * sqr(sin(kk[i]/2.0));
    }
    const double calpha = (1.0 + sqr(wp) + p2) / 2.0 / wp;
    const double alpha = std::acosh(calpha);
    const double lwa = 1.0 - wp * exp(-alpha);
    SpinMatrix m; setUnit(m, mass * lwa);
    SpinMatrix ipgm = pg;
    ipgm *= -ii;
    ipgm += m;
    ipgm *= lwa / (p2 + sqr(mass * lwa));
    if (1.0e-10 > p2 && 1.0e-10 > lwa) {
      // if (0.0 != norm(ipgm)) {
      //   Display(cname, fname, "kg = %s\n", show(kg).c_str());
      //   Display(cname, fname, "p2         = %13.5E\n", p2);
      //   Display(cname, fname, "wp         = %13.5E\n", wp);
      //   Display(cname, fname, "alpha      = %13.5E\n", alpha);
      //   Display(cname, fname, "lwa        = %13.5E\n", lwa);
      //   Display(cname, fname, "norm(ipgm) = %13.5E\n", norm(ipgm));
      // }
      setZero(sp4d.getElem(kl));
      continue;
    }
    ipgm *= sp4d.getElem(kl);
    sp4d.getElem(kl) = ipgm;
  }
}

inline void propSpinPropagator4d(SpinPropagator4d& sp4d, const double mass, const std::array<double,DIM>& momtwist)
{
  TIMER_VERBOSE("propSpinPropagator4d");
  const Geometry& geo = sp4d.geo;
  fftComplexField(sp4d, true);
  propMomSpinPropagator4d(sp4d, mass, momtwist);
  fftComplexField(sp4d, false);
  sp4d *= 1.0 / geo.totalVolume();
}

inline void setPointSourcePlusM(QedGaugeField& f, const Complex& coef, const Coordinate& xg, const int mu)
{
  TIMER("setPointSourcePlusM");
  const Geometry& geo = f.geo;
  Coordinate xl; geo.coordinateLfG(xl, xg);
  if (geo.isLocal(xl)) {
    f.getElem(xl,mu) += coef;
  }
}

inline void setBoxSourcePlusM(SpinPropagator4d& f, const Complex& coef, const Coordinate& xg1, const Coordinate& xg2)
  // 0 <= xg1 <= xg2 <= totalSite
  // FIXME: Do not handle the cross boundary case very well.
{
  TIMER("setBoxSourcePlusM");
  const Geometry& geo = f.geo;
  SpinMatrix sm; setUnit(sm, coef);
#pragma omp parallel for
  for (long index = 0; index < geo.localVolume(); ++index) {
    Coordinate xl; geo.coordinateFromIndex(xl, index);
    Coordinate xg; geo.coordinateGfL(xg, xl);
    if (xg1 <= xg && xg < xg2) {
      f.getElem(xl) += sm;
    }
  }
}

inline void setWallSourcePlusM(SpinPropagator4d& f, const Complex& coef, const int t)
{
  TIMER("setWallSourcePlusM");
  const Geometry& geo = f.geo;
  assert(0 <= t && t < geo.totalSite(3));
  const Coordinate xg1({ 0, 0, 0, t });
  const Coordinate xg2({ geo.totalSite(0), geo.totalSite(1), geo.totalSite(2), t+1 });
  setBoxSourcePlusM(f, coef, xg1, xg2);
}

inline void sequentialPhotonSpinPropagatorPlusM(
    SpinPropagator4d& src, const Complex coef,
    const QedGaugeField& egf, const SpinPropagator4d& sol)
{
  TIMER("sequentialPhotonSpinPropagatorPlusM");
  const Geometry& geo = sol.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.localVolume(); ++index) {
    Coordinate xl; geo.coordinateFromIndex(xl, index);
    for (int mu = 0; mu < DIM; mu++) {
      // a = A_\mu(x)
      Complex a = egf.getElem(xl, mu);
      // tmp = \gamma_\mu \psi(x)
      SpinMatrix tmp = SpinMatrixConstants::getGamma(mu) * sol.getElem(xl);
      // tmp = coef * \gamma_\mu A_\mu(x) \psi(x)
      tmp *= a * coef;
      src.getElem(xl) += tmp;
    }
  }
}

inline SpinMatrix contractSpinPropagator4d(const SpinPropagator4d& snk, const SpinPropagator4d& src)
{
  TIMER("contractSpinPropagator");
  const Geometry& geo = src.geo;
  SpinMatrix sum; setZero(sum);
#pragma omp parallel
  {
    SpinMatrix psum; setZero(psum);
#pragma omp for nowait
    for (long index = 0; index < geo.localVolume(); ++index) {
      Coordinate xl; geo.coordinateFromIndex(xl, index);
      psum += snk.getElem(xl).adjoint() * src.getElem(xl);
    }
    for (int i = 0; i < omp_get_num_threads(); ++i) {
#pragma omp barrier
      if (omp_get_thread_num() == i) {
          sum += psum;
      }
    }
  }
  sumVector(Vector<double>((double*)&sum, sizeof(SpinMatrix)/sizeof(double)));
  return sum;
}

QLAT_END_NAMESPACE
