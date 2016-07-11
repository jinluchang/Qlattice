#include <qlat/qlat.h>

#include <iostream>
#include <complex>

using namespace std;
using namespace qlat;

void test1()
{
  TIMER("test1");
  Geometry geo;
  Coordinate totalSite(16, 16, 16, 32);
  geo.init(totalSite, 1);
  RngField rf;
  const int seed = 1231;
  const int type = 1;
  const int traj = 1;
  rf.init(geo, seed, type, traj);
  for (long traj = 0; traj < 16; ++traj) {
    Complex a = 0.0;
    for (long index = 0; index < geo.localVolume(); ++index) {
      Coordinate x; geo.coordinateFromIndex(x, index);
      RngState& rs = rf.getElem(x);
      a += polar(1.0, uRandGen(rs, PI, -PI));
    }
    double sum = norm(a);
    double sigma2 = sqr(norm(a));
    glbSum(sum);
    glbSum(sigma2);
    const int Ni = geo.localVolume();
    const int Ntake = 1;
    const int Nb = geo.geon.numNode;
    if (0 == getIdNode()) {
      cout << "traj     : " << traj << endl;
      cout << "Expected : " << Ni * Ntake << endl;
      cout << "Mean     : " << sum / Nb << endl;
      cout << "Var      : " << sqrt(sigma2 / Nb - sqr(sum / Nb)) / sqrt (Nb) << endl;
    }
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
