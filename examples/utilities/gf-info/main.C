#include <qlat/qlat.h>

#include <complex>
#include <iostream>
#include <vector>

int main(int argc, char* argv[])
{
  using namespace qlat;
  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1, 2));
  size_node_list.push_back(Coordinate(1, 1, 1, 4));
  size_node_list.push_back(Coordinate(1, 1, 1, 8));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
  begin(&argc, &argv, size_node_list);
  if (argc < 2) {
    displayln_info("usage: ./gf-info fn");
    displayln_info("usage: ./gf-info fn wilson/iwasaki max_iter delta_t n_steps");
    displayln_info("usage: ./gf-info fn ape max_iter");
    exit(-1);
  }
  const std::string fn = remove_trailing_slashes(argv[1]);
  std::string type;
  if (argc >= 3) {
    type = argv[2];
  }
  int max_iter = 10;
  if (argc >= 4) {
    max_iter = read_long(argv[3]);
  }
  double delta_t = 1.0;
  if (argc >= 5) {
    delta_t = read_double(argv[4]);
  }
  long n_steps = 100;
  if (argc >= 6) {
    n_steps = read_long(argv[5]);
  }
  displayln_info(ssprintf("fn = '%s'", fn.c_str()));
  GaugeField gf, gf0;
  load_gauge_field(gf0, fn);
  if (type == "" or type == "wilson") {
    gf = gf0;
    displayln_info("Wilson flow info");
    display_gauge_field_info_table_with_wilson_flow("", gf, delta_t, n_steps, max_iter);
  }
  if (type == "" or type == "iwasaki") {
    gf = gf0;
    displayln_info("Wilson flow info with c1=-0.331 (Iwasaki action)");
    display_gauge_field_info_table_with_wilson_flow("", gf, delta_t, n_steps, max_iter,
                                                    -0.331);
  }
  if (type == "" or type == "ape") {
    gf = gf0;
    displayln_info("APE smearing (alpha=0.5) info");
    display_gauge_field_info_table_with_ape_smear("", gf, 0.5, max_iter);
  }
  Timer::display();
  end();
  return 0;
}
