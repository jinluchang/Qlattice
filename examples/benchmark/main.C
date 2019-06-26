#include <qlat/qlat.h>

#include <complex>
#include <iostream>
#include <vector>

using namespace qlat;
using namespace std;

void test_get_data()
{
  TIMER("test_get_data");
  // const size_t size = 2 * 1024 * 1024;
  const size_t size = 2 * 1024;
  std::vector<double> data_send(size), data_recv(size);
  for (int i = 0; i < 16; ++i) {
    TIMER_VERBOSE_FLOPS(
        "get_data");  // transferred bi-direction added in unit of Bytes
    timer.flops += size * sizeof(double) * 2 * get_num_node() *
                   4;  // 2: two direction, 4: four transfers
    get_data_dir(Vector<double>(data_recv), Vector<double>(data_send), 0);
    get_data_dir(Vector<double>(data_send), Vector<double>(data_recv), 0);
    get_data_dir(Vector<double>(data_recv), Vector<double>(data_send), 0);
    get_data_dir(Vector<double>(data_send), Vector<double>(data_recv), 0);
  }
}

void test_gf_fft()
{
  TIMER("test_gf_fft");
  // Coordinate total_site(48, 48, 48, 96);
  // Coordinate total_site(16, 16, 16, 32);
  Coordinate total_site(4, 4, 4, 8);
  RngState rs(getGlobalRngState(), "test_fft");
  Geometry geo;
  geo.init(total_site, 1);
  GaugeField gf;
  gf.init(geo);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    RngState rsi(rs, get_id_node() * geo.local_volume() + index);
    Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = gf.get_elems(xl);
    for (int m = 0; m < v.size(); ++m) {
      ColorMatrix& cm = v[m];
      for (int i = 0; i < NUM_COLOR; ++i) {
        for (int j = 0; j < NUM_COLOR; ++j) {
          cm(i, j) = Complex(gRandGen(rsi), gRandGen(rsi));
        }
      }
    }
  }
  const long flops = get_data_size(gf) * get_num_node();
  for (int i = 0; i < 16; ++i) {
    TIMER_VERBOSE_FLOPS("fft-gf");
    timer.flops += flops * 2;
    fft_complex_field(gf, true);
    fft_complex_field(gf, false);
  }
}

void test_fft()
{
  TIMER("test_fft");
  // Coordinate total_site(48, 48, 48, 96);
  // Coordinate total_site(16, 16, 16, 32);
  Coordinate total_site(4, 4, 4, 8);
  RngState rs(getGlobalRngState(), "test_fft");
  Geometry geo;
  geo.init(total_site, 1);
  FieldM<Complex, 12> f, ff;
  f.init(geo);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    RngState rsi(rs, get_id_node() * geo.local_volume() + index);
    Coordinate xl = geo.coordinate_from_index(index);
    Vector<Complex> v = f.get_elems(xl);
    for (int m = 0; m < v.size(); ++m) {
      v[m] = Complex(gRandGen(rsi), gRandGen(rsi));
    }
  }
  ff = f;
  const long flops = get_data_size(f) * get_num_node();
  const double total_volume = product(geo.total_site());
  for (int i = 0; i < 16; ++i) {
    TIMER_VERBOSE_FLOPS("fft");
    timer.flops += flops * 2;
    fft_complex_field(f, true);
    fft_complex_field(f, false);
    f *= 1.0 / total_volume;
  }
  ff -= f;
  displayln_info(fname + ssprintf(": qnorm(f) = %24.17E ; qnorm(diff) = %24.17E", qnorm(f), qnorm(ff)));
}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  test_get_data();
  test_gf_fft();
  test_fft();
  Timer::display();
  end();
  return 0;
}
