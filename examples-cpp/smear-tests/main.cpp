#include <qlat/qlat.h>

#include <complex>
#include <iostream>
#include <qlat/vector_utils/utils_smear_vecs.h>

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
  GaugeField gf;
  gf.init(geo);
  set_g_rand_color_matrix_field(gf, RngState(rs, "gf-0.1"), 0.1);
  gf_show_info(gf, 1);
  GaugeField gf1;
  //

  Propagator4d prop_src;
  Propagator4d prop_qlat;
  Propagator4d prop_vec;
  prop_src.init(geo);
  prop_qlat.init(geo);
  prop_vec.init(geo);
  set_g_rand(prop_src, RngState(rs, "prop-0.1"));
  prop_qlat = prop_src;
  prop_vec  = prop_src;

  {
    TIMER_VERBOSE("test-field_shift");

    int nsmear   =  10;
    double width = 2.0;

    CoordinateD mom;
    bool smear_in_time_dir = false;
    for(int i=0;i<4;i++){mom[i] = 0.0;}

    const double coef =  3.0*width*width/(2*nsmear);
    prop_smear_qlat_convension(prop_vec, gf, coef, nsmear, mom, smear_in_time_dir);

    set_left_expanded_gauge_field(gf1, gf);
    prop_smear(prop_qlat, gf1, coef, nsmear, mom, smear_in_time_dir);

    ////prop_vec -= prop_qlat;
    displayln_info(ssprintf("CHECK: Consistency: orig qnorm: %.10E ; smear qnorm %.10E ; new smear qnorm: %.10E",
                            qnorm(prop_src), qnorm(prop_qlat), qnorm(prop_vec)));
  }

}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  get_global_rng_state() = RngState(get_global_rng_state(), "qlat-smear-tests");
  simple_tests();
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
