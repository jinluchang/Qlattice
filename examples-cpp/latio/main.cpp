#include "qlat-setup.h"

namespace qlat
{  //

inline void demo_c()
{
  const Long n_tsep = 4;
  const Long n_op = 2;
  LatData ld;
  ld.info.push_back(lat_dim_number("tsep", 0, n_tsep - 1));
  ld.info.push_back(lat_dim_number("op", 0, n_op - 1));
  ld.info.push_back(lat_dim_string("type", make_array<std::string>("a", "b")));
  ld.info.push_back(lat_dim_re_im());
  lat_data_alloc(ld);
  set_zero(ld);
  for (int i = 0; i < n_tsep; ++i) {
    for (int j = 0; j < n_op; ++j) {
      Vector<ComplexD> v = lat_data_cget(ld, make_array<Long>(i, j));
      v[0] = (ComplexD)i + ii * (ComplexD)j;
      v[1] = (ComplexD)j + ii * (ComplexD)i;
    }
  }
  ld.save("results-data-c.lat");
  LatData ld1;
  ld1.load("results-data-c.lat");
  displayln_info(show(ld1));
}

inline void demo_r()
{
  const Long n_tsep = 4;
  const Long n_op = 2;
  LatData ld;
  ld.info.push_back(lat_dim_number("tsep", 0, n_tsep - 1));
  ld.info.push_back(lat_dim_number("op", 0, n_op - 1));
  ld.info.push_back(lat_dim_string("type", make_array<std::string>("a", "b")));
  lat_data_alloc(ld);
  set_zero(ld);
  for (int i = 0; i < n_tsep; ++i) {
    for (int j = 0; j < n_op; ++j) {
      Vector<double> v = lat_data_get(ld, make_array<Long>(i, j));
      v[0] = i;
      v[1] = j;
    }
  }
  ld.save("results-data-r.lat");
  LatData ld1;
  ld1.load("results-data-r.lat");
  displayln_info(show(ld1));
}

}  // namespace qlat

int main(int argc, char* argv[])
{
  using namespace qlat;
  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1, 1));
  begin(&argc, &argv, size_node_list);
  display_geometry_node();
  setup_log_idx();
  setup();
  demo_c();
  demo_r();
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
