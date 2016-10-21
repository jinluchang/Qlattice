#include <qlat/qlat.h>

#include <iostream>
#include <complex>
#include <vector>

using namespace qlat;
using namespace std;

void test1()
{
  TIMER("test1");
  const size_t size = 2 * 1024 * 1024;
  std::vector<double> dataSend(size), dataRecv(size);
  for (int i = 0; i < 16; ++i) {
    TIMER_VERBOSE_FLOPS("get_data"); // transferred bi-direction added in unit of Bytes
    timer.flops += size * sizeof(double) * 2 * get_num_node() * 4; // 2: two direction, 4: four transfers
    get_data_dir(Vector<double>(dataRecv), Vector<double>(dataSend), 0);
    get_data_dir(Vector<double>(dataSend), Vector<double>(dataRecv), 0);
    get_data_dir(Vector<double>(dataRecv), Vector<double>(dataSend), 0);
    get_data_dir(Vector<double>(dataSend), Vector<double>(dataRecv), 0);
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
