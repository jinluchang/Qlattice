#include <qutils/lat-io.h>

namespace qlat
{  //

inline void demo()
{
  const long n_tsep = 4;
  const long n_op = 2;
  LatData ld;
  ld.info.push_back(lat_dim_number("tsep", 0, n_tsep - 1));
  ld.info.push_back(lat_dim_number("op", 0, n_op - 1));
  ld.info.push_back(lat_dim_re_im());
  lat_data_alloc(ld);
  set_zero(ld);
  LatData ld1;
  ld1 = ld;
  for (int i = 0; i < n_tsep; ++i) {
    Vector<Complex> v = lat_data_complex_get(ld, make_array<long>(i));
    for (int j = 0; j < n_op; ++j) {
      v[j] = (Complex)i + ii * (Complex)j;
    }
  }
  ld.save("results-data.lat");
  ld1.load("results-data.lat");
  ld1.save("results-data1.lat");
  print(ld1);
}

}  // namespace qlat

int main(int argc, char* argv[])
{
  demo();
  return 0;
}
