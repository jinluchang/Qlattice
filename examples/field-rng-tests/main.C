#include <qlat/qlat.h>

#include <iostream>
#include <complex>

using namespace qlat;
using namespace std;

void coordinateHalf(Coordinate& xh, const Coordinate& x)
{
  for (int mu = 0; mu < DIM; ++mu) {
    xh[mu] = x[mu] / 2;
  }
}

void test1()
{
  TIMER("test1");
  Coordinate totalSite(16, 16, 16, 32);
  Geometry geo;
  geo.init(totalSite, 1);
  Coordinate totalSiteHalf; coordinateHalf(totalSiteHalf, totalSite);
  Geometry geoHalf;
  geoHalf.init(totalSiteHalf, 1);
  const int seed = 1231;
  const int type = 1;
  const int traj = 1;
  RngField rf;
  rf.init(geo, seed, type, traj);
  FieldM<Complex,1> af;
  FieldM<double,1> sumf;
  FieldM<double,1> sigma2f;
  af.init(geoHalf);
  sumf.init(geoHalf);
  sigma2f.init(geoHalf);
  double gsum = 0.0;
  double gsigma2 = 0.0;
  const int Ni = 2*2*2*2;
  const int Ntake = 1;
  const int Nb = geo.totalVolume() / Ni;
  const int Ntraj = 16;
  for (long traj = 0; traj < Ntraj; ++traj) {
    setZero(af);
    setZero(sumf);
    setZero(sigma2f);
    for (long index = 0; index < geo.localVolume(); ++index) {
      Coordinate x; geo.coordinateFromIndex(x, index);
      Coordinate xh; coordinateHalf(xh, x);
      RngState& rs = rf.getElem(x);
      af.getElem(xh) += polar(1.0, uRandGen(rs, PI, -PI));
    }
    for (long index = 0; index < geoHalf.localVolume(); ++index) {
      Coordinate x; geoHalf.coordinateFromIndex(x, index);
      Complex& a = af.getElem(x);
      sumf.getElem(x) += norm(a);
      sigma2f.getElem(x) += sqr(norm(a));
    }
    double sum;
    fieldGlbSumDouble(Vector<double>(sum), sumf);
    double sigma2;
    fieldGlbSumDouble(Vector<double>(sigma2), sigma2f);
    gsum += sum / Nb;
    gsigma2 += sqr(sum / Nb);
    if (0 == getIdNode()) {
      cout << "traj     : " << traj << endl;
      cout << "Expected : " << Ni * Ntake << endl;
      cout << "Mean     : " << sum / Nb << endl;
      cout << "Var      : " << sqrt(sigma2 / Nb - sqr(sum / Nb)) / sqrt(Nb-1) << endl;
    }
  }
  if (0 == getIdNode()) {
    cout << "# Final" << endl;
    cout << "Expected : " << Ni * Ntake << endl;
    cout << "Mean     : " << gsum / Ntraj << endl;
    cout << "Var      : " << sqrt(gsigma2 / Ntraj - sqr(gsum / Ntraj)) / sqrt(Ntraj-1) << endl;
  }
}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  test1();
  Timer::display();
  end();
  return 0;
}
