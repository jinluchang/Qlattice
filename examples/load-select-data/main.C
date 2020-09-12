#include "compute-wall-src-info.h"
#include "compute-wall-src-prop-norm-ratio.h"
#include "compute-check-prop.h"
#include "data-load.h"

namespace qlat
{  //

inline void compute_traj(const std::string& job_tag, const int traj)
{
  setup(job_tag);
  TIMER_VERBOSE("compute_traj");
  // SADJUST ME
  // check_all_prop_psrc_exact(job_tag, traj);
  // check_prop_data(job_tag, traj);
  compute_wall_src_info(job_tag, traj, 0);
  compute_wall_src_info(job_tag, traj, 1);
  compute_wall_src_prop_norm_ratio(job_tag, traj);
  //
  clear_all_data_cache();
}

inline void compute(const std::string& job_tag)
{
  TIMER_VERBOSE("compute");
  const std::vector<int> trajs = get_data_trajs(job_tag);
  for (int i = 0; i < (int)trajs.size(); ++i) {
    const int traj = trajs[i];
    compute_traj(job_tag, traj);
    if (get_total_time() > 1.0) {
      Timer::display();
    }
    Timer::reset();
  }
}

}  // namespace qlat

int main(int argc, char* argv[])
{
  using namespace qlat;
  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1, 1));
  size_node_list.push_back(Coordinate(1, 1, 1, 2));
  size_node_list.push_back(Coordinate(1, 1, 1, 4));
  size_node_list.push_back(Coordinate(1, 1, 1, 8));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
  size_node_list.push_back(Coordinate(1, 1, 2, 16));
  size_node_list.push_back(Coordinate(1, 2, 2, 16));
  size_node_list.push_back(Coordinate(2, 2, 2, 16));
  begin(&argc, &argv, size_node_list);
  setup_log_idx();
  setup();
  qmkdir_info(ssprintf("data"));
  std::vector<std::string> job_tags;
  // SADJUST ME
  job_tags.push_back("24D");
  job_tags.push_back("32D");
  job_tags.push_back("24DH");
  job_tags.push_back("32Dfine");
  job_tags.push_back("48I");
  job_tags.push_back("64I");
  //
  for (int k = 0; k < (int)job_tags.size(); ++k) {
    const std::string& job_tag = job_tags[k];
    compute(job_tag);
  }
  end();
  return 0;
}
