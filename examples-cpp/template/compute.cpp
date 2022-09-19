#include "qlat-setup.h"

namespace qlat
{  //

typedef InverterDomainWall Inverter;

// -----------------------------------------------------------------------------------

bool compute_traj_do(const std::string& job_tag, const int traj)
{
  TIMER_VERBOSE("compute_traj_do");
  displayln_info(fname +
                 ssprintf(": job_tag='%s' ; traj=%d", job_tag.c_str(), traj));
  const std::string job_path = get_job_path(job_tag, traj);
  const RngState rs = RngState("seed").split(job_tag).split(traj);
  const Coordinate total_site = get_total_site(job_tag);
  Geometry geo;
  geo.init(total_site, 1);
  //
  GaugeField gf;
  load_configuration(gf, job_tag, traj);
  gf_show_info(gf);
  qassert(is_matching_geo(gf.geo(), geo));
  //
  GaugeTransform gt;
  gt.init(geo);
  set_g_rand_color_matrix_field(gt, rs.split("gt-rs"), 1.0);
  //
  LatData ld;
  ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
  ld.info.push_back(lat_dim_re_im());
  lat_data_alloc(ld);
  set_zero(ld);
  LatData ld1;
  ld1 = ld;
  qassert(is_matching(ld, ld1));
  for (int t = 0; t < total_site[3]; ++t) {
    lat_data_complex_get(ld, make_array<int>())[t] = t;
    lat_data_complex_get(ld1, make_array(t))[0] = t;
  }
  if (0 == get_id_node()) {
    print(ld);
    print(ld1);
  }
  qassert(is_matching(ld, ld1));
  qassert(ld.res == ld1.res);
  lat_data_save_info(job_path + "/data.lat", ld);
  //
  benchmark_deflate(geo, 8, Coordinate(2, 2, 2, 2), 100, 50,
                    rs.split("benchmark_deflate"));
  //
  const FermionAction fa = get_fermion_actions(job_tag)[0];
  LowModes lm;
  // load_or_compute_low_modes(lm, get_low_modes_path(job_tag, traj), gf, fa,
  // get_lanc_arg(job_tag));
  Inverter inv;
  setup_inverter(inv, gf, fa, lm);
  //
  Propagator4d prop;
  const int tslice = 0 % total_site[3];
  set_wall_src_propagator(prop, GaugeTransformInverter<Inverter>(inv, gt),
                          tslice, CoordinateD());
  displayln_info(ssprintf("prop qnorm = %24.17E", qnorm(prop)));
  //
  displayln_info(
      ssprintf("prop qnorm = %24.17E", qnorm(prop.get_elem(Coordinate()))));
  //
  find_max_eigen_value_hermop_sym2(
      inv, RngState(rs, "find_max_eigen_value_hermop_sym2"), 5);
  //
  qtouch_info(job_path + "/checkpoint.txt");
  return false;
}

// -----------------------------------------------------------------------------------

}  // namespace qlat
