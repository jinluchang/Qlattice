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
}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  test1();
  Timer::display();
  end();
  return 0;
}
