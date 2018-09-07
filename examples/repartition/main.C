#include <qlat/qlat.h>

#include <iostream>
#include <complex>
#include <vector>

using namespace qlat;
using namespace std;

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  if (argc < 2) {
    displayln_info("usage: ./repartition 1x1x1x8 location1 location2 ...");
    exit(-1);
  }
  const Coordinate new_size_node = read_coordinate(argv[1]);
  displayln_info("new_size_node: " + show(new_size_node));
  std::vector<std::string> fns(argc-2);
  for (long i = 0; i < (long)fns.size(); ++i) {
    fns[i] = remove_trailing_slashes(argv[2+i]);
    displayln_info(ssprintf("fns[%5d/%d] = '%s'", i, fns.size(), fns[i].c_str()));
  }
  displayln_info("Start to repartition...");
  for (long i = 0; i < (long)fns.size(); ++i) {
    displayln_info(ssprintf("fns[%5d/%d] = '%s'", i, fns.size(), fns[i].c_str()));
    dist_repartition(new_size_node, fns[i]);
  }
  Timer::display();
  end();
  return 0;
}
