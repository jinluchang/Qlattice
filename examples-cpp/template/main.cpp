#include "qlat-setup.h"

namespace qlat
{  //

bool compute_traj_do(const std::string& job_tag, const int traj);

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
  const array<std::string, 1> job_tags =
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
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
