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

template <int N>
void test_fft()
{
  TIMER("test_fft");
  // Coordinate total_site(48, 48, 48, 96);
  // Coordinate total_site(32, 32, 32, 64);
  // Coordinate total_site(16, 16, 16, 32);
  Coordinate total_site(4, 4, 4, 8);
  RngState rs(getGlobalRngState(), "test_fft");
  Geometry geo;
  geo.init(total_site, 1);
  FieldM<Complex, N> f, ff;
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
  const double total_volume = geo.total_volume();
  displayln_info(fname + ssprintf(": %d", N));
  Timer::max_call_times_for_always_show_info() = 0;
  Timer::display();
  Timer::reset();
  for (int i = 0; i < 16; ++i) {
    TIMER_VERBOSE_FLOPS("fft");
    timer.flops += flops * 2;
    fft_complex_field(f, true);
    fft_complex_field(f, false);
    f *= 1.0 / total_volume;
  }
  Timer::display();
  displayln_info(fname + ssprintf(": total_site=%s N=%d fft", show(total_site).c_str(), N));
  Timer::reset();
  for (int i = 0; i < 16; ++i) {
    TIMER_VERBOSE_FLOPS("fft-spatial");
    timer.flops += flops * 2;
    fft_complex_field_spatial(f, true);
    fft_complex_field_spatial(f, false);
    f *= 1.0 / (total_site[0] * total_site[1] * total_site[2]);
  }
  Timer::display();
  displayln_info(fname + ssprintf(": total_site=%s N=%d fft-spatial", show(total_site).c_str(), N));
  Timer::reset();
  for (int i = 0; i < 16; ++i) {
    TIMER_VERBOSE_FLOPS("fft-dir-x");
    timer.flops += flops * 2;
    fft_complex_field_dir(f, 0, true);
    fft_complex_field_dir(f, 0, false);
    f *= 1.0 / total_site[0];
  }
  Timer::display();
  displayln_info(fname + ssprintf(": total_site=%s N=%d fft-x", show(total_site).c_str(), N));
  Timer::reset();
  for (int i = 0; i < 16; ++i) {
    TIMER_VERBOSE_FLOPS("fft-dir-y");
    timer.flops += flops * 2;
    fft_complex_field_dir(f, 1, true);
    fft_complex_field_dir(f, 1, false);
    f *= 1.0 / total_site[1];
  }
  Timer::display();
  displayln_info(fname + ssprintf(": total_site=%s N=%d fft-y", show(total_site).c_str(), N));
  Timer::reset();
  for (int i = 0; i < 16; ++i) {
    TIMER_VERBOSE_FLOPS("fft-dir-z");
    timer.flops += flops * 2;
    fft_complex_field_dir(f, 2, true);
    fft_complex_field_dir(f, 2, false);
    f *= 1.0 / total_site[2];
  }
  Timer::display();
  displayln_info(fname + ssprintf(": %d", N));
  displayln_info(fname + ssprintf(": total_site=%s N=%d fft-z", show(total_site).c_str(), N));
  Timer::reset();
  for (int i = 0; i < 16; ++i) {
    TIMER_VERBOSE_FLOPS("fft-dir-t");
    timer.flops += flops * 2;
    fft_complex_field_dir(f, 3, true);
    fft_complex_field_dir(f, 3, false);
    f *= 1.0 / total_site[3];
  }
  Timer::display();
  displayln_info(fname + ssprintf(": total_site=%s N=%d fft-t", show(total_site).c_str(), N));
  Timer::reset();
  ff -= f;
  displayln_info(fname + ssprintf(": qnorm(f) = %24.17E ; qnorm(diff) = %24.17E", qnorm(f), qnorm(ff)));
  fft_complex_field(f, true);
  displayln_info(fname + ssprintf(": crc32 of fft(f) = %06X.", field_dist_crc32(f)));
  fft_complex_field(f, false);
}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  test_get_data();
  test_gf_fft();
  test_fft<1>();
  test_fft<2>();
  test_fft<4>();
  test_fft<8>();
  test_fft<12>();
  test_fft<16>();
  test_fft<20>();
  test_fft<32>();
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
