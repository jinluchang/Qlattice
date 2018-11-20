#include "setup.h"

QLAT_START_NAMESPACE


// -----------------------------------------------------------------------------------

inline bool compute_traj(const std::string& job_tag, const int traj)
{
  TIMER_VERBOSE("compute_traj");
  const std::string job_path = get_job_path(job_tag, traj);
  qmkdir_info(job_path);
  qmkdir_sync_node(job_path + "/logs");
  switch_monitor_file_info(job_path + ssprintf("/logs/%010ld.txt", get_log_idx()));
  displayln_info(fname + ssprintf(": job_tag='%s' ; traj=%d", job_tag.c_str(), traj));
  //
  const RngState rs = RngState("seed").split(job_tag).split(traj);
  //
  const Geometry geo = get_geo(job_tag);
  //
  GaugeField gf;
  load_configuration(gf, job_tag, traj);
  gf_show_info(gf);
  qassert(is_matching_geo(gf.geo, geo));
  //
  const FermionAction fa = get_fermion_actions(job_tag)[0];
  LowModes lm;
  // load_or_compute_low_modes(lm, get_low_modes_path(job_tag, traj), gf, fa, get_lanc_arg(job_tag));
  InverterDomainWall inv;
  if (lm.initialized) {
    setup_inverter(inv, gf, fa, lm);
  } else {
    setup_inverter(inv, gf, fa);
  }
  //
  GaugeTransform gt;
  gt.init(geo);
  set_g_rand_color_matrix_field(gt, rs.split("gt-rs"), 1.0);
  //
  Propagator4d prop;
  const int tslice = 0;
  set_wall_src_propagator(prop, GaugeTransformInverter<InverterDomainWall>(inv, gt), tslice, CoordinateD());
  displayln_info(ssprintf("prop norm = %24.17E", norm(prop)));
  //
  qtouch_info(job_path + "/checkpoint.txt");
  return false;
}

inline bool compute(const std::string& job_tag)
{
  TIMER_VERBOSE("compute");
  setup(job_tag);
  const std::vector<int> trajs = get_trajs(job_tag);
  for (int i = 0; i < (int)trajs.size(); ++i) {
    const int traj = trajs[i];
    if (does_file_exist_sync_node(get_job_path(job_tag, traj) + "/checkpoint.txt")) {
      continue;
    }
    // if (does_file_exist_sync_node(get_low_modes_path(job_tag, traj) + "/checkpoint")) {
    //   continue;
    // }
    if (obtain_lock(get_job_path(job_tag, traj) + "-lock")) {
      compute_traj(job_tag, traj);
      release_lock();
    }
    Timer::display();
  }
  return false;
}

// -----------------------------------------------------------------------------------

QLAT_END_NAMESPACE

int main(int argc, char* argv[])
{
  using namespace qlat;
  const std::string job_tag = "free-4nt8";
  begin(&argc, &argv);
  setup_log_idx();
  if (not compute(job_tag)) {
    displayln_info("program finished successfully.");
  }
  Timer::display();
  end();
  return 0;
}
