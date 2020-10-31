#pragma once

#include "compute-utils.h"
#include "data-load.h"
#include "qlat/contract-field.h"

namespace qlat
{  //

inline std::string get_meson_vv_path(const std::string& job_tag, const int traj)
{
  return ssprintf("analysis/field-meson-vv/%s/results=%d", job_tag.c_str(),
                  traj);
}

inline void compute_meson_vv_type(const std::string& job_tag, const int traj,
                                  const int type1, const int type2,
                                  const int type3)
{
  check_sigint();
  check_time_limit();
  Timer::autodisplay();
  TIMER_VERBOSE("compute_meson_vv_type");
  // TODO
}

inline void compute_meson_vv(const std::string& job_tag, const int traj)
{
  check_sigint();
  Timer::autodisplay();
  TIMER_VERBOSE("compute_meson_vv");
  const std::string path = get_meson_vv_path(job_tag, traj);
  const std::string path_checkpoint = path + "/checkpoint.txt";
  if (does_file_exist_sync_node(path_checkpoint)) {
    return;
  }
  if (not(check_wall_src_props(job_tag, traj, 0) and
          check_wall_src_props(job_tag, traj, 1))) {
    return;
  }
  check_sigint();
  check_time_limit();
  if (not obtain_lock(ssprintf("lock-meson-vv-%s-%d", job_tag.c_str(), traj))) {
    return;
  }
  setup(job_tag, traj);
  qmkdir_info("analysis/field-meson-vv");
  qmkdir_info(ssprintf("analysis/field-meson-vv/%s", job_tag.c_str()));
  qmkdir_info(path);
  // TODO
  qtouch_info(path_checkpoint);
  release_lock();
}

}  // namespace qlat
