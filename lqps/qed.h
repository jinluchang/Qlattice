#pragma once

#include <lqps/field.h>
#include <lqps/field-fft.h>

#include <Eigen/Dense>

LQPS_START_NAMESPACE

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

struct SpinMatrix
{
  Eigen::Matrix<Complex,4,4> matrix;
  //
  Complex& operator()(int i, int j)
  {
    return matrix(i,j);
  }
  const Complex& operator()(int i, int j) const
  {
    return matrix(i,j);
  }
};

inline void setZero(SpinMatrix& sm)
{
  sm.matrix.setZero();
}

struct SpinPropagator4d : FieldM<SpinMatrix,1>
{
  virtual const char* cname()
  {
    return "SpinPropagator4d";
  }
};

inline void propMomPhotonInvert(QedGaugeField& qgf, const std::array<double,4>& momtwist)
  // Feynman Gauge
  // All spatial zero mode removed.
  // qgf in momentum space.
{
  TIMER("propMomPhotonInvert");
  const Geometry& geo = qgf.geo;
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
      for (int mu = 0; mu < DIM; ++mu) {
        qgf.getElem(kl, mu) *= 0.0;
      }
    } else {
      double s2inv = 1.0 / s2;
      for (int mu = 0; mu < geo.multiplicity; mu++) {
        qgf.getElem(kl, mu) *= s2inv;
      }
    }
  }
}

inline void propPhotonInvert(QedGaugeField& qgf, const std::array<double,4>& momtwist)
  // Feynman Gauge
  // All spatial zero mode removed.
  // qgf in coordinate space.
{
  TIMER("propPhotonInvert");
  const Geometry& geo = qgf.geo;
  fftComplexField(qgf, true);
  propMomPhotonInvert(qgf, momtwist);
  fftComplexField(qgf, false);
  qgf *= 1.0 / geo.totalVolume();
}

inline void propMomComplexScalerInvert(ComplexScalerField& csf, const double mass, const std::array<double,4>& momtwist)
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
    for (int mu = 0; mu < geo.multiplicity; mu++) {
      csf.getElem(kl, mu) *= s2inv;
    }
  }
}

inline void propComplexScalerInvert(ComplexScalerField& csf, const double mass, const std::array<double,4>& momtwist)
{
  TIMER("propComplexScalerInvert");
  const Geometry& geo = csf.geo;
  fftComplexField(csf, true);
  propMomComplexScalerInvert(csf, mass, momtwist);
  fftComplexField(csf, false);
  csf *= 1.0 / geo.totalVolume();
}

LQPS_END_NAMESPACE
