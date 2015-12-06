#include <lqps/lqps.h>

#include <mpi.h>

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
  lqps::Coordinate totalSite({ 2, 2, 2, 4 });
  lqps::Geometry geo; geo.init(totalSite, 1);
  DisplayInfo(cname, fname, "geo=\n%s\n", lqps::show(geo).c_str());
  lqps::QedGaugeField f1; f1.init(geo);
  lqps::QedGaugeField f2; f2.init(geo);
  lqps::QedGaugeField f3; f3.init(geo);
  lqps::setZero(f1);
  lqps::setZero(f2);
  lqps::setZero(f3);
  setField(f1);
  f2 = f1;
  DisplayInfo(cname, fname, "norm(f1) = %.16E\n", norm(f1));
  lqps::fftComplexField(f1, true);
  f1 *= 1.0 / sqrt((double)(geo.localVolume() * geo.geon.numNode));
  DisplayInfo(cname, fname, "norm(f1) = %.16E\n", norm(f1));
  f3 = f2;
  f3 -= f1;
  DisplayInfo(cname, fname, "norm(f3) = %.16E\n", norm(f3));
  lqps::fftComplexField(f1, false);
  f1 *= 1.0 / sqrt((double)(geo.localVolume() * geo.geon.numNode));
  DisplayInfo(cname, fname, "norm(f1) = %.16E\n", norm(f1));
  f3 = f2;
  f3 -= f1;
  DisplayInfo(cname, fname, "norm(f3) = %.16E\n", norm(f3));
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
