#pragma once

#include "data-load.h"

namespace qlat
{  //

inline std::string get_three_point_func_path(const std::string& job_tag,
                                           const int traj)
{
  return ssprintf("analysis/lat-three-point/%s/results=%d", job_tag.c_str(),
                  traj);
}

inline void compute_three_point_func_type(const std::string& job_tag,
                                          const int traj, const int type1,
                                          const int type2, const int type3)
{
  check_sigterm();
  check_time_limit();
  Timer::autodisplay();
  TIMER_VERBOSE("compute_three_point_func_type");
  const std::string path = get_three_point_func_path(job_tag, traj);
  const std::string path_three_point =
      path + ssprintf("/three-point-%d-%d-%d.lat", type1, type2, type3);
  if (not does_file_exist_sync_node(path_three_point)) {
    const WallSrcProps& wsp1 = get_wall_src_props(job_tag, traj, type1);
    const WallSrcProps& wsp2 = get_wall_src_props(job_tag, traj, type2);
    const WallSrcProps& wsp3 = get_wall_src_props(job_tag, traj, type3);
    const FieldSelection& fsel = get_field_selection(job_tag, traj);
    const Coordinate total_site = get_total_site(job_tag);
    const int yt_measurement_sparsity = 1;
    qassert(total_site[3] % yt_measurement_sparsity == 0);
    RngState rs = RngState("three-point-func")
                      .split("yt-start-time-slice")
                      .split(job_tag)
                      .split(traj);
    const int yt_measurement_start = rand_gen(rs) % yt_measurement_sparsity;
    const LatData ld_three_point_func = contract_three_point_function(
        wsp1, wsp2, wsp3, fsel, yt_measurement_sparsity, yt_measurement_start);
    lat_data_save_info(path_three_point, ld_three_point_func);
  }
}

inline void compute_three_point_func(const std::string& job_tag, const int traj)
{
  check_sigterm();
  Timer::autodisplay();
  TIMER_VERBOSE("compute_three_point_func");
  const std::string path = get_three_point_func_path(job_tag, traj);
  const std::string path_checkpoint = path + "/checkpoint.txt";
  if (does_file_exist_sync_node(path_checkpoint)) {
    return;
  }
  if (not(check_wall_src_props(job_tag, traj, 0) and
          check_wall_src_props(job_tag, traj, 1))) {
    return;
  }
  check_sigterm();
  check_time_limit();
  if (not obtain_lock(
          ssprintf("lock-three-point-func-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  setup(job_tag, traj);
  qmkdir_info("analysis/lat-three-point");
  qmkdir_info(ssprintf("analysis/lat-three-point/%s", job_tag.c_str()));
  qmkdir_info(path);
  compute_three_point_func_type(job_tag, traj, 0, 0, 0);
  compute_three_point_func_type(job_tag, traj, 0, 0, 1);
  compute_three_point_func_type(job_tag, traj, 1, 1, 0);
  compute_three_point_func_type(job_tag, traj, 1, 1, 1);
  compute_three_point_func_type(job_tag, traj, 0, 1, 0);
  compute_three_point_func_type(job_tag, traj, 0, 1, 1);
  compute_three_point_func_type(job_tag, traj, 1, 0, 0);
  compute_three_point_func_type(job_tag, traj, 1, 0, 1);
  qtouch_info(path_checkpoint);
  release_lock();
}

inline void compute_three_point_func_light(const std::string& job_tag, const int traj)
{
  check_sigterm();
  Timer::autodisplay();
  TIMER_VERBOSE("compute_three_point_func_light");
  const std::string path = get_three_point_func_path(job_tag, traj);
  const std::string path_checkpoint = path + "/checkpoint.txt";
  if (does_file_exist_sync_node(path_checkpoint)) {
    return;
  }
  if (does_file_exist_sync_node(path + "/three-point-0-0-0.lat")) {
    return;
  }
  if (not check_wall_src_props(job_tag, traj, 0)) {
    return;
  }
  check_sigterm();
  check_time_limit();
  if (not obtain_lock(
          ssprintf("lock-three-point-func-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  setup(job_tag, traj);
  qmkdir_info("analysis/lat-three-point");
  qmkdir_info(ssprintf("analysis/lat-three-point/%s", job_tag.c_str()));
  qmkdir_info(path);
  compute_three_point_func_type(job_tag, traj, 0, 0, 0);
  release_lock();
}

}  // namespace qlat
