#include <qlat/qlat.h>

#include <complex>
#include <iostream>
#include <qlat/vector_utils/utils_FFT_GPU.h>

using namespace qlat;

void simple_tests()
{
  TIMER_VERBOSE("simple_smear");
  RngState rs(get_global_rng_state(), fname);
  // const Coordinate total_site(16, 16, 16, 32);
  // const Coordinate total_site(8, 8, 8, 8);
  const Coordinate total_site(4, 4, 4, 8);
  Geometry geo;
  geo.init(total_site);

  qlat::FieldM<qlat::ComplexD, 12> src;src.init(geo);
  set_g_rand_double(src, RngState(rs, "prop-0.1"));

  std::vector< qlat::FieldM<qlat::ComplexD, 12> > cpuF;cpuF.resize(1);
  std::vector< qlat::FieldM<qlat::ComplexD, 12> > gpuF;gpuF.resize(1);
  cpuF[0] = src;
  gpuF[0] = src;

  {
    TIMER_VERBOSE("test-fft");
    qlat::fft_complex_field_spatial(cpuF[0], false);
    fft_fieldM(gpuF, false, false);
    displayln_info(ssprintf("CHECK: Consistency: orig qnorm: %.10E ; fft qnorm %.10E ; new fft qnorm: %.10E",
                            qnorm(src), qnorm(cpuF[0]), qnorm(gpuF[0])));
  }

}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  get_global_rng_state() = RngState(get_global_rng_state(), "qlat-fft-tests");
  simple_tests();
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
