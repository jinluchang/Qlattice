#include <lqps/lqps.h>

#include <iostream>

const char* cname = "Main";

void setField(lqps::Field<lqps::Complex>& f)
{
  TIMER("setField");
  const lqps::Geometry& geo = f.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.localVolume(); ++index) {
    lqps::Coordinate x; geo.coordinateFromIndex(x, index);
    lqps::Vector<lqps::Complex> fx = f.getElems(x);
    for (int m = 0; m < geo.multiplicity; ++m) {
      fx[m] = geo.geon.idNode * sqrt(2) + index * sqrt(3) * lqps::ii + m * sqrt(5);
    }
  }
}

lqps::SpinMatrix projPositiveState(const lqps::SpinMatrix& x)
{
  const lqps::SpinMatrix psm = (lqps::SpinMatrixConstants::getUnit() + lqps::SpinMatrixConstants::getGamma(3)) / 2.0;
  return psm * x * psm;
}

lqps::SpinMatrix lblMuonLine(const int tsnk, const int tsrc,
    const lqps::QedGaugeField& egf1, const lqps::QedGaugeField& egf2, const lqps::QedGaugeField& egf3,
    const double mass, const std::array<double,4>& momtwist)
{
  TIMER("lblMuonLine");
  const lqps::Geometry& geo = egf1.geo;
  lqps::SpinPropagator4d sol; sol.init(geo);
  lqps::SpinPropagator4d src; src.init(geo);
  lqps::setZero(src);
  lqps::setWallSourcePlusM(src, 1.0, tsrc);
  lqps::propSpinPropagator4d(src, mass, momtwist);
  //
  sol = src;
  setZero(src);
  lqps::sequentialPhotonSpinPropagatorPlusM(src, lqps::ii, egf1, sol);
  lqps::propSpinPropagator4d(src, mass, momtwist);
  //
  sol = src;
  setZero(src);
  lqps::sequentialPhotonSpinPropagatorPlusM(src, lqps::ii, egf2, sol);
  lqps::propSpinPropagator4d(src, mass, momtwist);
  //
  sol = src;
  setZero(src);
  lqps::sequentialPhotonSpinPropagatorPlusM(src, lqps::ii, egf3, sol);
  lqps::propSpinPropagator4d(src, mass, momtwist);
  //
  lqps::SpinPropagator4d snk; snk.init(geo);
  lqps::setZero(snk);
  lqps::setWallSourcePlusM(snk, 1.0, tsnk);
  return lqps::contractSpinPropagator4d(snk, src);
}

lqps::SpinMatrix lblMuonLineC(const int tsnk, const int tsrc,
    const lqps::QedGaugeField& egf1, const lqps::QedGaugeField& egf2, const lqps::QedGaugeField& egf3,
    const double mass, const std::array<double,4>& momtwist)
{
  TIMER("lblMuonLineC");
  lqps::SpinMatrix sm; lqps::setZero(sm);
  sm += lblMuonLine(tsnk, tsrc, egf1, egf2, egf3, mass, momtwist);
  sm += lblMuonLine(tsnk, tsrc, egf2, egf3, egf1, mass, momtwist);
  sm += lblMuonLine(tsnk, tsrc, egf3, egf1, egf2, mass, momtwist);
  sm += lblMuonLine(tsnk, tsrc, egf3, egf2, egf1, mass, momtwist);
  sm += lblMuonLine(tsnk, tsrc, egf2, egf1, egf3, mass, momtwist);
  sm += lblMuonLine(tsnk, tsrc, egf1, egf3, egf2, mass, momtwist);
  return sm;
}

lqps::SpinMatrix lblMuonPartPointSrc(const lqps::Geometry& geo, const int tsnk, const int tsrc,
    const lqps::Coordinate& xg1, const int mu1,
    const lqps::Coordinate& xg2, const int mu2,
    const lqps::Coordinate& xg3, const int mu3,
    const double mass, const std::array<double,4>& momtwist)
{
  TIMER("lblMuonPartPointSrc");
  lqps::QedGaugeField egf1; egf1.init(geo);
  lqps::QedGaugeField egf2; egf2.init(geo);
  lqps::QedGaugeField egf3; egf3.init(geo);
  lqps::setZero(egf1);
  lqps::setZero(egf2);
  lqps::setZero(egf3);
  lqps::setPointSourcePlusM(egf1, 1.0, xg1, mu1);
  lqps::setPointSourcePlusM(egf2, 1.0, xg2, mu2);
  lqps::setPointSourcePlusM(egf3, 1.0, xg3, mu3);
  lqps::propPhotonInvert(egf1, momtwist);
  lqps::propPhotonInvert(egf2, momtwist);
  lqps::propPhotonInvert(egf3, momtwist);
  return lblMuonLineC(tsnk, tsrc, egf1, egf2, egf3, mass, momtwist);
}

void lblMagneticMomentSpinMatrix(lqps::Array<lqps::SpinMatrix,3> bs, const lqps::Geometry& geo, const int tsnk, const int tsrc,
    const double mass, const std::array<double,4>& momtwist)
  // pretend to the operator to be
  // \Sigma_i * mass / 2
{
  TIMER("lblMagneticMomentSpinMatrix");
  lqps::SpinPropagator4d snk; snk.init(geo);
  lqps::SpinPropagator4d src; src.init(geo);
  lqps::setZero(snk);
  lqps::setZero(src);
  lqps::setWallSourcePlusM(snk, 1.0, tsnk);
  lqps::setWallSourcePlusM(src, 1.0, tsrc);
  lqps::propSpinPropagator4d(snk, mass, momtwist);
  lqps::propSpinPropagator4d(src, mass, momtwist);
  const int top = lqps::mod(tsrc + lqps::mod(tsnk - tsrc, geo.totalSite(3)) / 2, geo.totalSite(3));
  lqps::Coordinate xgop({ 0, 0, 0, top });
  lqps::Coordinate xlop; geo.coordinateLfG(xlop, xgop);
  lqps::setZero(bs);
  if (geo.isLocal(xlop)) {
    for (int i = 0; i < 3; ++i) {
      bs[i] = lqps::SpinMatrixConstants::getGamma5() * snk.getElem(xlop).adjoint() * lqps::SpinMatrixConstants::getGamma5();
      bs[i] *= lqps::SpinMatrixConstants::getCapSigma(i) * mass / 2.0;
      bs[i] *= src.getElem(xlop);
      bs[i] = projPositiveState(bs[i]);
    }
  }
  lqps::sumVector(lqps::Vector<double>((double*)bs.data(), getDataSize(bs)/sizeof(double)));
}

lqps::Complex linearFit(const lqps::SpinMatrix& x, const lqps::SpinMatrix& base)
{
  return (conj(base.array()) * x.array()).sum() / (conj(base.array()) * base.array()).sum();
}

void lblShowMuonPartPointSrc(const lqps::Geometry& geo, const int tsnk, const int tsrc,
    const lqps::Array<lqps::SpinMatrix,3>& bs,
    const lqps::Coordinate& xg1, const int mu1,
    const lqps::Coordinate& xg2, const int mu2,
    const lqps::Coordinate& xg3, const int mu3,
    const double mass, const std::array<double,4>& momtwist)
{
  TIMER("lblShowMuonPartPointSrc");
  DisplayInfo(cname, fname, "mu1 = %d ; mu2 = %d ; mu3 = %d .\n", mu1, mu2, mu3);
  lqps::SpinMatrix muonline = lblMuonPartPointSrc(geo, tsnk, tsrc, xg1, mu1, xg2, mu2, xg3, mu3, mass, momtwist);
  muonline = projPositiveState(muonline);
  DisplayInfo(cname, fname, "norm(muonline) = %.16e\n", lqps::norm(muonline));
  DisplayInfo(cname, fname, "muonline =\n%s\n", lqps::show(muonline).c_str());
  DisplayInfo(cname, fname, "linearFit[0] = %s\n", lqps::show(linearFit(muonline,bs[0])).c_str());
  DisplayInfo(cname, fname, "linearFit[1] = %s\n", lqps::show(linearFit(muonline,bs[1])).c_str());
  DisplayInfo(cname, fname, "linearFit[2] = %s\n", lqps::show(linearFit(muonline,bs[2])).c_str());
}

void lblMuonPart()
{
  TIMER("lblMuonPart");
  lqps::Coordinate totalSite({ 6, 6, 6, 128 });
  lqps::Geometry geo; geo.init(totalSite, 1);
  DisplayInfo(cname, fname, "geo =\n%s\n", lqps::show(geo).c_str());
  std::array<double,4> momtwist({ 0.0, 0.0, 0.0, 0.0 });
  const double mass = 0.1;
  const int tsnk = geo.totalSite(3)/4*3;
  const int tsrc = geo.totalSite(3)/4;
  std::array<lqps::SpinMatrix,3> bs;
  lblMagneticMomentSpinMatrix(bs, geo, tsnk, tsrc, mass, momtwist);
  DisplayInfo(cname, fname, "bs[0] =\n%s\n", lqps::show(bs[0]).c_str());
  DisplayInfo(cname, fname, "bs[1] =\n%s\n", lqps::show(bs[1]).c_str());
  DisplayInfo(cname, fname, "bs[2] =\n%s\n", lqps::show(bs[2]).c_str());
  lqps::Coordinate xg1({ 0, 2, 4, geo.totalSite(3)/2 - 4 });
  lqps::Coordinate xg2({ 1, 0, 2, geo.totalSite(3)/2 + 0 });
  lqps::Coordinate xg3({ 0, 3, 5, geo.totalSite(3)/2 + 4 });
  for (int mu1 = 0; mu1 < 4; ++mu1) {
    for (int mu2 = 0; mu2 < 4; ++mu2) {
      for (int mu3 = 0; mu3 < 4; ++mu3) {
        lblShowMuonPartPointSrc(geo, tsnk, tsrc, bs, xg1, mu1, xg2, mu2, xg3, mu3, mass, momtwist);
      }
    }
  }
}

int main(int argc, char* argv[])
{
  lqps::Coordinate lsizeNode({ 1, 1, 1, 2 });
  lqps::begin(&argc, &argv, lsizeNode);
  lblMuonPart();
  Timer::display();
  lqps::end();
  return 0;
}
