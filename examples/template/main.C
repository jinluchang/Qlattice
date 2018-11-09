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
  const Geometry geo = get_geo(job_tag);
  GaugeField gf;
  load_configuration(gf, job_tag, traj);
  gf_show_info(gf);
  qassert(is_matching_geo(geo, gf.geo));
  // TODO
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
