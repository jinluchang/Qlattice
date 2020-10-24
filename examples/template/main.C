#include "qlat-setup.h"

namespace qlat
{  //

typedef InverterDomainWall Inverter;

// -----------------------------------------------------------------------------------

inline bool compute_traj_do(const std::string& job_tag, const int traj)
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
  qassert(is_matching_geo(gf.geo, geo));
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
  ld.save(job_path + "/data.lat");
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
  displayln_info(ssprintf("prop qnorm = %24.17E", qnorm(prop.get_elem(Coordinate()))));
  //
  find_max_eigen_value_hermop_sym2(
      inv, RngState(rs, "find_max_eigen_value_hermop_sym2"), 5);
  //
  qtouch_info(job_path + "/checkpoint.txt");
  return false;
}

inline bool compute_traj(const std::string& job_tag, const int traj)
{
  setup(job_tag);
  TIMER_VERBOSE("compute_traj");
  displayln_info(fname + ssprintf(": Checking '%s'.",
                                  get_job_path(job_tag, traj).c_str()));
  if (does_file_exist_sync_node(get_job_path(job_tag, traj) +
                                "/checkpoint.txt")) {
    displayln_info(fname + ssprintf(": Finished '%s'.",
                                    get_job_path(job_tag, traj).c_str()));
    return false;
  }
  if (get_config_fn(job_tag, traj) == "") {
    displayln_info(fname + ssprintf(": No config '%s'.",
                                    get_job_path(job_tag, traj).c_str()));
    return false;
  }
  // if (does_file_exist_sync_node(get_low_modes_path(job_tag, traj) +
  //                               "/checkpoint")) {
  //   displayln_info(fname + ssprintf(": No low modes '%s'.",
  //                                   get_job_path(job_tag, traj).c_str()));
  //   return false;
  // }
  if (not obtain_lock(get_job_path(job_tag, traj) + "-lock")) {
    displayln_info(fname + ssprintf(": Cannot obtain lock '%s'.",
                                    get_job_path(job_tag, traj).c_str()));
    return true;
  }
  displayln_info(fname + ssprintf(": Start computing '%s'.",
                                  get_job_path(job_tag, traj).c_str()));
  setup(job_tag, traj);
  const bool is_failed = compute_traj_do(job_tag, traj);
  release_lock();
  return is_failed;
}

inline bool compute(const std::string& job_tag)
{
  bool is_failed = false;
  const std::vector<int> trajs = get_trajs(job_tag);
  for (int i = 0; i < (int)trajs.size(); ++i) {
    const int traj = trajs[i];
    is_failed = compute_traj(job_tag, traj) or is_failed;
    if (get_total_time() > 1.0) {
      Timer::display();
    }
    Timer::reset();
  }
  return is_failed;
}

// -----------------------------------------------------------------------------------

}  // namespace qlat

int main(int argc, char* argv[])
{
  using namespace qlat;
  const std::array<std::string, 1> job_tags =
      make_array<std::string>("free-4nt8");
  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(2,2,2,4));
  size_node_list.push_back(Coordinate(1,1,1,4));
  begin(&argc, &argv, size_node_list);
  display_geometry_node();
  setup_log_idx();
  setup();
  for (int k = 0; k < (int)job_tags.size(); ++k) {
    const std::string& job_tag = job_tags[k];
    if (not compute(job_tag)) {
      displayln_info(ssprintf("program finished successfully for '%s'.", job_tag.c_str()));
    }
  }
  end();
  return 0;
}
