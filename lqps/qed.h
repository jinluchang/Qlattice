#pragma once

#include <lqps/field.h>
#include <lqps/field-fft.h>

LQPS_START_NAMESPACE

struct QedGaugeField : FieldM<Complex,DIM>
{
  virtual const char* cname()
  {
    return "QedGaugeField";
  }
};

void propMomPhotonInvert(QedGaugeField& qgf, const std::array<double,4>& momtwist)
  // Feynman Gauge
  // All spatial zero mode removed.
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
      double s2inv = 1.0 / s2 / geo.localVolume() / geo.geon.numNode;
      for (int mu = 0; mu < 4; mu++) {
        qgf.getElem(kl, mu) *= s2inv;
      }
    }
  }
}

LQPS_END_NAMESPACE
