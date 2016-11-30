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
  displayln_info(shows("plaq : ") + show(gf_avg_plaq(gf)));
  displayln_info(shows("trace: ") + show(gf_avg_link_trace(gf)));
  set_g_rand_color_matrix_field(gf, RngState(rs, "gf-0.1"), 0.1);
  displayln_info(shows("plaq : ") + show(gf_avg_plaq(gf)));
  displayln_info(shows("trace: ") + show(gf_avg_link_trace(gf)));
  set_g_rand_color_matrix_field(gf, RngState(rs, "gf-0.3"), 0.3);
  displayln_info(shows("plaq : ") + show(gf_avg_plaq(gf)));
  displayln_info(shows("trace: ") + show(gf_avg_link_trace(gf)));
  GaugeTransform gt;
  gt.init(geo);
  set_g_rand_color_matrix_field(gt, RngState(rs, "gt-1.0"), 1.0);
  make_apply_gauge_transformation(gf, gf, gt);
  displayln_info(shows("plaq : ") + show(gf_avg_plaq(gf)));
  displayln_info(shows("trace: ") + show(gf_avg_link_trace(gf)));
  make_apply_gauge_transformation(gf, gf, gt);
  displayln_info(shows("plaq : ") + show(gf_avg_plaq(gf)));
  displayln_info(shows("trace: ") + show(gf_avg_link_trace(gf)));
  unitarize(gf);
  displayln_info(shows("plaq : ") + show(gf_avg_plaq(gf)));
  displayln_info(shows("trace: ") + show(gf_avg_link_trace(gf)));
  make_tree_gauge_transformation(gt, gf);
  make_apply_gauge_transformation(gf, gf, gt);
  displayln_info(shows("plaq : ") + show(gf_avg_plaq(gf)));
  displayln_info(shows("trace: ") + show(gf_avg_link_trace(gf)));
  gf_ape_smear(gf, gf, 0.1);
  displayln_info(shows("plaq : ") + show(gf_avg_plaq(gf)));
  displayln_info(shows("trace: ") + show(gf_avg_link_trace(gf)));
  gf_hyp_smear(gf, gf, 0.75, 0.6, 0.3);
  displayln_info(shows("plaq : ") + show(gf_avg_plaq(gf)));
  displayln_info(shows("trace: ") + show(gf_avg_link_trace(gf)));
}

int main(int argc, char* argv[])
{
  begin(&argc, &argv);
  get_global_rng_state() = RngState(get_global_rng_state(), "qcd-utils-tests");
  simple_tests();
  end();
  Timer::display();
  return 0;
}
