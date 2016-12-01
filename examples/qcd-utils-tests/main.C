#include <qlat/qlat.h>

#include <iostream>
#include <complex>

using namespace qlat;
using namespace std;

void simple_tests()
{
  TIMER_VERBOSE("simple_tests");
  RngState rs(get_global_rng_state(), fname);
  // const Coordinate total_site(16, 16, 16, 32);
  const Coordinate total_site(8, 8, 8, 8);
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
  make_apply_gauge_transformation(gf, gf, gt);
  gf_show_info(gf);
  make_apply_gauge_transformation(gf, gf, gt);
  gf_show_info(gf);
  unitarize(gf);
  gf_show_info(gf);
  make_tree_gauge_transformation(gt, gf);
  make_apply_gauge_transformation(gf, gf, gt);
  gf_show_info(gf);
  GaugeField gfs;
  gf_ape_smear(gfs, gf, 0.1);
  gf_show_info(gfs, 1);
  gf_hyp_smear(gfs, gf, 0.1, 0.0, 0.3);
  gf_show_info(gfs);
  gf_hyp_smear(gfs, gf, 0.75, 0.6, 0.3);
  gf_show_info(gfs, 1);
}

int main(int argc, char* argv[])
{
  Timer::max_call_times_for_always_show_info() = 0;
  begin(&argc, &argv);
  get_global_rng_state() = RngState(get_global_rng_state(), "qcd-utils-tests");
  simple_tests();
  end();
  Timer::display();
  return 0;
}
