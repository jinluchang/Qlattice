#include <qutils/lat-io.h>

namespace qlat
{  //

inline void demo_c()
{
  const long n_tsep = 4;
  const long n_op = 2;
  LatData ld;
  ld.info.push_back(lat_dim_number("tsep", 0, n_tsep - 1));
  ld.info.push_back(lat_dim_number("op", 0, n_op - 1));
  ld.info.push_back(lat_dim_string("type", make_array<std::string>("a", "b")));
  ld.info.push_back(lat_dim_re_im());
  lat_data_alloc(ld);
  set_zero(ld);
  for (int i = 0; i < n_tsep; ++i) {
    for (int j = 0; j < n_op; ++j) {
      Vector<Complex> v = lat_data_cget(ld, make_array<long>(i, j));
      v[0] = (Complex)i + ii * (Complex)j;
      v[1] = (Complex)j + ii * (Complex)i;
    }
  }
  ld.save("results-data-c.lat");
  LatData ld1;
  ld1.load("results-data-c.lat");
  displayln_info(show(ld1));
}

inline void demo_r()
{
  const long n_tsep = 4;
  const long n_op = 2;
  LatData ld;
  ld.info.push_back(lat_dim_number("tsep", 0, n_tsep - 1));
  ld.info.push_back(lat_dim_number("op", 0, n_op - 1));
  ld.info.push_back(lat_dim_string("type", make_array<std::string>("a", "b")));
  lat_data_alloc(ld);
  set_zero(ld);
  for (int i = 0; i < n_tsep; ++i) {
    for (int j = 0; j < n_op; ++j) {
      Vector<double> v = lat_data_get(ld, make_array<long>(i, j));
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
  demo_c();
  demo_r();
  displayln_info("CHECK: finished successfully.");
  return 0;
}
