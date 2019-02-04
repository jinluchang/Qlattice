#include <qlat/qlat.h>

#include <complex>
#include <iostream>

using namespace qlat;
using namespace std;

void simple_tests()
{
  TIMER_VERBOSE("simple_tests");
  RngState rs(get_global_rng_state(), fname);
  // const Coordinate total_site(16, 16, 16, 32);
  // const Coordinate total_site(8, 8, 8, 8);
  const Coordinate total_site(4, 4, 4, 8);
  Geometry geo;
  geo.init(total_site, 1);
  GaugeField gf;
  gf.init(geo);
  set_unit(gf);
  gf_show_info(gf);
  set_g_rand_color_matrix_field(gf, RngState(rs, "gf-0.1"), 0.1);
  gf_show_info(gf, 1);
  set_g_rand_color_matrix_field(gf, RngState(rs, "gf-0.3"), 0.3);
  gf_show_info(gf, 1);
  GaugeTransform gt;
  gt.init(geo);
  set_g_rand_color_matrix_field(gt, RngState(rs, "gt-1.0"), 1.0);
  gf_apply_gauge_transformation(gf, gf, gt);
  gf_show_info(gf, 1);
  gf_apply_gauge_transformation(gf, gf, gt);
  gf_show_info(gf);
  unitarize(gf);
  gf_show_info(gf, 1);
  make_tree_gauge_transformation(gt, gf);
  gf_apply_gauge_transformation(gf, gf, gt);
  gf_show_info(gf);
  GaugeField gfs;
  gf_ape_smear(gfs, gf, 0.1);
  gf_show_info(gfs, 1);
  gf_hyp_smear(gfs, gf, 0.1, 0.0, 0.3);
  gf_show_info(gfs);
  gf_hyp_smear(gfs, gf, 0.75, 0.6, 0.3);
  gf_show_info(gfs, 1);
  //
  GaugeField gf1, gf2;
  const Coordinate c1(1, 2, 3, 4);
  const Coordinate c2(7, 3, 1, 8);
  field_shift(gf1, gf, c1);
  field_shift(gf2, gf1, c2);
  field_shift(gf1, gf, c1 + c2);
  gf2 -= gf1;
  displayln_info(ssprintf("orig norm: %E ; shift norm %E ; diff norm: %E",
                          norm(gf), norm(gf1), norm(gf2)));
}

void show_matrix()
{
  const std::array<SpinMatrix, 16>& cps_gms =
      SpinMatrixConstants::get_cps_gms();
  for (int i = 0; i < 16; ++i) {
    displayln_info(ssprintf("cps_gms[%d] =\n", i) + show(cps_gms[i]));
  }
}

int main(int argc, char* argv[])
{
  Timer::max_function_name_length_shown() = 50;
  begin(&argc, &argv);
  get_global_rng_state() = RngState(get_global_rng_state(), "qcd-utils-tests");
  show_matrix();
  simple_tests();
  end();
  Timer::display();
  return 0;
}
