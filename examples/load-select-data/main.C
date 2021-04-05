#include "compute-check-prop.h"
#include "compute-chvp.h"
#include "compute-meson-chvp.h"
#include "compute-meson-snk-src.h"
#include "compute-meson-vv-meson.h"
#include "compute-meson-vv.h"
#include "compute-psel-fsel-distribution.h"
#include "compute-three-point-func.h"
#include "compute-two-point-func.h"
#include "compute-wall-src-prop-norm-ratio.h"

namespace qlat
{  //

inline void compute_traj(const std::string& job_tag, const int traj)
{
  setup(job_tag);
  TIMER_VERBOSE("compute_traj");
  displayln_info(ssprintf("Computing %s %d", job_tag.c_str(), traj));
  // SADJUST ME
  // check_all_prop_psrc_exact(job_tag, traj);
  // check_prop_data(job_tag, traj);
  // compute_wall_src_info(job_tag, traj, 0);
  // compute_wall_src_info(job_tag, traj, 1);
  // compute_wall_src_prop_norm_ratio(job_tag, traj);
  //
  // SADJUST ME
  compute_meson_snk_src(job_tag, traj);
  compute_meson_chvp(job_tag, traj);
  compute_chvp(job_tag, traj);
  compute_meson_vv_meson(job_tag, traj);
  compute_meson_vv(job_tag, traj);
  compute_two_point_func(job_tag, traj);
  compute_three_point_func(job_tag, traj);
  //
  // SADJUST ME
  compute_meson_snk_src_light(job_tag, traj);
  compute_meson_chvp_light(job_tag, traj);
  compute_chvp_light(job_tag, traj);
  compute_meson_vv_meson_light(job_tag, traj);
  compute_meson_vv_light(job_tag, traj);
  compute_two_point_func_light(job_tag, traj);
  compute_three_point_func_light(job_tag, traj);
  //
  compute_psel_fsel_distribution(job_tag, traj);
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

inline void test()
{
  TIMER_VERBOSE("test");
  compute_traj("test-4nt16", 1010);
  compute_traj("24D", 1010);
  compute_traj("24D", 1900);
  compute_traj("24D", 1010);
  compute_traj("24D", 1900);
  Timer::display();
  sync_node();
  ssleep(3.0);
  exit(0);
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
  size_node_list.push_back(Coordinate(2, 2, 4, 16));
  size_node_list.push_back(Coordinate(2, 4, 4, 16));
  size_node_list.push_back(Coordinate(4, 4, 4, 16));
  begin(&argc, &argv, size_node_list);
  display_geometry_node();
  setup_log_idx();
  setup();
  qmkdir_info(ssprintf("analysis"));
  //
  // test();
  //
  std::vector<std::string> job_tags;
  // SADJUST ME
  job_tags.push_back("test-4nt16");
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
