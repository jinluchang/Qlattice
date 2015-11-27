#include <iostream>

#include <mpi.h>

#include <lqps/lqps.h>

const char* cname = "Main";

void lblMuonPart()
{
  TIMER("lblMuonPart");
  lqps::Coordinate totalSite({ 2, 2, 2, 4 });
  lqps::Geometry geo; geo.init(totalSite, 1);
  DisplayInfo(cname, fname, "geo=\n%s\n", lqps::show(geo).c_str());
}

int main(int argc, char* argv[])
{
  lqps::Coordinate lsizeNode({ 1, 1, 1, 2 });
  lqps::start(&argc, &argv, lsizeNode);
  lblMuonPart();
  Timer::display();
  lqps::end();
  return 0;
}
