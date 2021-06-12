#pragma once

#include "data-load.h"

namespace qlat
{  //

inline std::string get_meson_snk_src(const std::string& job_tag, const int traj)
{
  return ssprintf("analysis/lat-meson-snk-src/%s/results=%d", job_tag.c_str(),
                  traj);
}

inline void compute_meson_snk_src_type(const std::string& job_tag,
                                        const int traj, const int type1,
                                        const int type2)
{
  check_sigterm();
  Timer::autodisplay();
  const std::string path = get_meson_snk_src(job_tag, traj) +
                           ssprintf("/meson-snk-src-%d-%d.lat", type1, type2);
  if (does_file_exist_sync_node(path)) {
    return;
  }
  TIMER_VERBOSE("compute_meson_snk_src_type");
  const WallSrcProps& wsp1 = get_wall_src_props(job_tag, traj, type1);
  const WallSrcProps& wsp2 = get_wall_src_props(job_tag, traj, type2);
  const Coordinate total_site = get_total_site(job_tag);
  const LatData ld_1_2 = contract_meson_snk_src(wsp1, wsp2, total_site);
  lat_data_save_info(path, ld_1_2);
}

inline void compute_meson_snk_src(const std::string& job_tag, const int traj)
{
  check_sigterm();
  Timer::autodisplay();
  TIMER_VERBOSE("compute_meson_snk_src");
  const std::string path = get_meson_snk_src(job_tag, traj);
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
          ssprintf("lock-meson-snk-src-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  setup(job_tag, traj);
  qmkdir_info("analysis/lat-meson-snk-src");
  qmkdir_info(ssprintf("analysis/lat-meson-snk-src/%s", job_tag.c_str()));
  qmkdir_info(path);
  compute_meson_snk_src_type(job_tag, traj, 0, 0);
  compute_meson_snk_src_type(job_tag, traj, 0, 1);
  compute_meson_snk_src_type(job_tag, traj, 1, 0);
  compute_meson_snk_src_type(job_tag, traj, 1, 1);
  qtouch_info(path_checkpoint);
  release_lock();
}

inline void compute_meson_snk_src_light(const std::string& job_tag, const int traj)
{
  check_sigterm();
  Timer::autodisplay();
  TIMER_VERBOSE("compute_meson_snk_src_light");
  const std::string path = get_meson_snk_src(job_tag, traj);
  const std::string path_checkpoint = path + "/checkpoint.txt";
  if (does_file_exist_sync_node(path_checkpoint)) {
    return;
  }
  if (does_file_exist_sync_node(path + "/meson-snk-src-0-0.lat")) {
    return;
  }
  if (not check_wall_src_props(job_tag, traj, 0)) {
    return;
  }
  check_sigterm();
  check_time_limit();
  if (not obtain_lock(
          ssprintf("lock-meson-snk-src-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  setup(job_tag, traj);
  qmkdir_info("analysis/lat-meson-snk-src");
  qmkdir_info(ssprintf("analysis/lat-meson-snk-src/%s", job_tag.c_str()));
  qmkdir_info(path);
  compute_meson_snk_src_type(job_tag, traj, 0, 0);
  release_lock();
}

}  // namespace qlat
