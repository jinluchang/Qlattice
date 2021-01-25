#pragma once

#include "data-load.h"

namespace qlat
{  //

inline std::string get_two_point_func_path(const std::string& job_tag,
                                           const int traj)
{
  return ssprintf("analysis/lat-two-point/%s/results=%d", job_tag.c_str(),
                  traj);
}

inline void compute_two_point_func_type(const std::string& job_tag,
                                        const int traj, const int type1,
                                        const int type2)
{
  check_sigterm();
  Timer::autodisplay();
  TIMER_VERBOSE("compute_two_point_func_type");
  const WallSrcProps& wsp1 = get_wall_src_props(job_tag, traj, type1);
  const WallSrcProps& wsp2 = get_wall_src_props(job_tag, traj, type2);
  const FieldSelection& fsel = get_field_selection(job_tag, traj);
  LatData ld_two_point_func, ld_two_point_wall_snk_func,
      ld_two_point_wall_snk_func_sparse_corrected;
  const std::string path = get_two_point_func_path(job_tag, traj);
  const std::string path_two_point =
      path + ssprintf("/two-point-%d-%d.lat", type1, type2);
  const std::string path_two_point_wall_snk =
      path + ssprintf("/two-point-wall-snk-%d-%d.lat", type1, type2);
  const std::string path_two_point_wall_snk_sparse_corrected =
      path +
      ssprintf("/two-point-wall-snk-sparse-corrected-%d-%d.lat", type1, type2);
  if (does_file_exist_sync_node(path_two_point_wall_snk)) {
    ld_two_point_wall_snk_func = lat_data_load_info(path_two_point_wall_snk);
  } else {
    ld_two_point_wall_snk_func =
        contract_two_point_wall_snk_function(wsp1, wsp2, fsel);
    lat_data_save_info(path_two_point_wall_snk, ld_two_point_wall_snk_func);
  }
  if (does_file_exist_sync_node(path_two_point)) {
    ld_two_point_func = lat_data_load_info(path_two_point);
  } else {
    ld_two_point_func = contract_two_point_function(wsp1, wsp2, fsel);
    lat_data_save_info(path_two_point, ld_two_point_func);
  }
  if (does_file_exist_sync_node(path_two_point_wall_snk_sparse_corrected)) {
    ld_two_point_wall_snk_func_sparse_corrected =
        lat_data_load_info(path_two_point_wall_snk_sparse_corrected);
  } else {
    ld_two_point_wall_snk_func_sparse_corrected =
        contract_two_point_wall_snk_function(ld_two_point_wall_snk_func,
                                             ld_two_point_func, fsel);
    lat_data_save_info(path_two_point_wall_snk_sparse_corrected,
                       ld_two_point_wall_snk_func_sparse_corrected);
  }
}

inline void compute_two_point_func(const std::string& job_tag, const int traj)
{
  check_sigterm();
  Timer::autodisplay();
  TIMER_VERBOSE("compute_two_point_func");
  const std::string path = get_two_point_func_path(job_tag, traj);
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
          ssprintf("lock-two-point-func-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  setup(job_tag, traj);
  qmkdir_info("analysis/lat-two-point");
  qmkdir_info(ssprintf("analysis/lat-two-point/%s", job_tag.c_str()));
  qmkdir_info(path);
  compute_two_point_func_type(job_tag, traj, 0, 0);
  compute_two_point_func_type(job_tag, traj, 0, 1);
  compute_two_point_func_type(job_tag, traj, 1, 0);
  compute_two_point_func_type(job_tag, traj, 1, 1);
  qtouch_info(path_checkpoint);
  release_lock();
}

inline void compute_two_point_func_light(const std::string& job_tag, const int traj)
{
  check_sigterm();
  Timer::autodisplay();
  TIMER_VERBOSE("compute_two_point_func_light");
  const std::string path = get_two_point_func_path(job_tag, traj);
  const std::string path_checkpoint = path + "/checkpoint.txt";
  if (does_file_exist_sync_node(path_checkpoint)) {
    return;
  }
  if (does_file_exist_sync_node(
          path + "/two-point-wall-snk-sparse-corrected-0-0.lat")) {
    return;
  }
  if (not check_wall_src_props(job_tag, traj, 0)) {
    return;
  }
  check_sigterm();
  check_time_limit();
  if (not obtain_lock(
          ssprintf("lock-two-point-func-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  setup(job_tag, traj);
  qmkdir_info("analysis/lat-two-point");
  qmkdir_info(ssprintf("analysis/lat-two-point/%s", job_tag.c_str()));
  qmkdir_info(path);
  compute_two_point_func_type(job_tag, traj, 0, 0);
  release_lock();
}

}  // namespace qlat
