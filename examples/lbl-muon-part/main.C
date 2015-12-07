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

void lblMuonPart()
{
  TIMER("lblMuonPart");
  lqps::Coordinate totalSite({ 2, 2, 2, 16 });
  lqps::Geometry geo; geo.init(totalSite, 1);
  DisplayInfo(cname, fname, "geo=\n%s\n", lqps::show(geo).c_str());
  std::array<double,4> momtwist({ 0.0, 0.0, 0.0, 0.0 });
  lqps::Coordinate xg1({ 0, 0, 0, 0 });
  lqps::Coordinate xg2({ 0, 0, 0, 0 });
  lqps::Coordinate xg3({ 0, 0, 0, 0 });
  lqps::QedGaugeField egf1; egf1.init(geo);
  lqps::QedGaugeField egf2; egf2.init(geo);
  lqps::QedGaugeField egf3; egf3.init(geo);
  lqps::setPointSourcePlusM(egf1, 1.0, xg1, 0);
  lqps::setPointSourcePlusM(egf2, 1.0, xg2, 0);
  lqps::setPointSourcePlusM(egf3, 1.0, xg3, 0);
  //
  //
  //
  lqps::QedGaugeField f1; f1.init(geo);
  lqps::QedGaugeField f2; f2.init(geo);
  lqps::QedGaugeField f3; f3.init(geo);
  lqps::setZero(f1);
  lqps::setZero(f2);
  lqps::setZero(f3);
  setField(f1);
  f2 = f1;
  DisplayInfo(cname, fname, "norm(f1) = %.16E\n", lqps::norm(f1));
  lqps::fftComplexField(f1, true);
  f1 *= 1.0 / sqrt((double)geo.totalVolume());
  DisplayInfo(cname, fname, "norm(f1) = %.16E\n", lqps::norm(f1));
  f3 = f2;
  f3 -= f1;
  DisplayInfo(cname, fname, "norm(f3) = %.16E\n", lqps::norm(f3));
  lqps::fftComplexField(f1, false);
  f1 *= 1.0 / sqrt((double)geo.totalVolume());
  DisplayInfo(cname, fname, "norm(f1) = %.16E\n", lqps::norm(f1));
  f3 = f2;
  f3 -= f1;
  DisplayInfo(cname, fname, "norm(f3) = %.16E\n", lqps::norm(f3));
  f3 = f2;
  lqps::propPhotonInvert(f3, { 0.0, 0.0, 0.0, 0.0 });
  DisplayInfo(cname, fname, "norm(f3) = %.16E\n", lqps::norm(f3));
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
