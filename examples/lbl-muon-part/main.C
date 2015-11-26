#include <iostream>

#include <mpi.h>

#include <lqps/lqps.h>

void lblMuonPart()
{
  TIMER("lblMuonPart");
  std::cout << "hello world!" << std::endl;
}

int main(int argc, char* argv[])
{
  int lsizeNode[4] = { 1, 1, 1, 2 };
  lqps::start(&argc, &argv, lsizeNode);
  lblMuonPart();
  Timer::display();
  lqps::end();
  return 0;
}
