#include <qlat/qlat.h>

#include <complex>
#include <iostream>
#include <vector>

using namespace qlat;
using namespace std;

int main(int argc, char* argv[])
{
  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1, 1));
  get_output_level() = 0;
  begin(&argc, &argv, size_node_list);
  get_output_level() = UINT64_MAX;
  if (argc != 2) {
    displayln_info("usage: ./latio-show filename");
    exit(-1);
  }
  const std::string fn(argv[1]);
  LatData ld;
  ld.load(fn);
  display_info(show(ld));
  get_output_level() = 0;
  end();
  get_output_level() = UINT64_MAX;
  return 0;
}
