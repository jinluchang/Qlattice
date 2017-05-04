#include <qlat/qlat.h>

#include <iostream>
#include <complex>
#include <vector>

using namespace qlat;
using namespace std;

void test_io()
{
  TIMER("test_io");
  // Coordinate total_site(48, 48, 48, 96);
  Coordinate total_site(16, 16, 16, 32);
  crc32_check();
}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  test_io();
  Timer::display();
  end();
  return 0;
}
