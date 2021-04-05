#include <qlat/qlat.h>

#include <complex>
#include <iostream>
#include <vector>

using namespace qlat;
using namespace std;

int main(int argc, char* argv[])
{
  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1, 2));
  size_node_list.push_back(Coordinate(1, 1, 1, 4));
  size_node_list.push_back(Coordinate(1, 1, 1, 8));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
  begin(&argc, &argv, size_node_list);
  if (argc < 2) {
    displayln_info("usage: ./repartition 1x1x1x8 location1 location2 ...");
    exit(-1);
  }
  const Coordinate new_size_node = read_coordinate(argv[1]);
  displayln_info("new_size_node: " + show(new_size_node));
  std::vector<std::string> fns(argc - 2);
  for (long i = 0; i < (long)fns.size(); ++i) {
    fns[i] = remove_trailing_slashes(argv[2 + i]);
  }
  {
    double time = get_time();
    glb_sum(time);
    const RngState rs_permute(show(time));
    random_permute(fns, rs_permute);
  }
  for (long i = 0; i < (long)fns.size(); ++i) {
    displayln_info(
        ssprintf("fns[%5d/%d] = '%s'", i, fns.size(), fns[i].c_str()));
  }
  install_qhandle_sig();
  displayln_info("Start to repartition...");
  for (long i = 0; i < (long)fns.size(); ++i) {
    if (is_sigterm_received() > 0) {
      Timer::display();
      displayln_info("quit: because sigint received.");
      end();
      return 0;
    }
    check_time_limit();
    Timer::autodisplay();
    TIMER_VERBOSE("repartition-iter");
    displayln_info(
        ssprintf("fns[%5d/%d] = '%s'", i, fns.size(), fns[i].c_str()));
    if (does_file_exist_sync_node(fns[i] + "/metadata.txt")) {
      eigen_system_repartition(new_size_node, fns[i]);
    } else if (is_field(fns[i]) or is_dist_field(fns[i])) {
      dist_repartition(new_size_node, fns[i]);
    } else {
      displayln_info("Cannot repartition this data: '" + fns[i] + "'.");
    }
  }
  Timer::display();
  end();
  return 0;
}
