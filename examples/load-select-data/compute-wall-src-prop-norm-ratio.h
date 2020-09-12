#pragma once

#include "compute-wall-src-info.h"
#include "data-load.h"

namespace qlat
{  //

inline std::string get_wall_src_prop_norm_ratio_path(const std::string& job_tag,
                                                     const int traj)
{
  return ssprintf("analysis/wall-src-prop-norm-ratio/%s/results=%d",
                  job_tag.c_str(), traj);
}

inline void compute_wall_src_prop_norm_ratio(const std::string& job_tag,
                                             const int traj)
{
  const std::string rpath = get_wall_src_prop_norm_ratio_path(job_tag, traj);
  if (does_file_exist_sync_node(rpath + "/checkpoint.txt")) {
    return;
  }
  bool is_complete = true;
  for (int type = 0; type < 2; ++type) {
    const std::string type_tag = type == 0 ? "light" : "strange";
    const std::string rpath_e =
        rpath + ssprintf("/norm-exact-%s.lat", type_tag.c_str());
    const std::string rpath_s =
        rpath + ssprintf("/norm-sloppy-%s.lat", type_tag.c_str());
    if (does_file_exist_sync_node(rpath_e) and
        does_file_exist_sync_node(rpath_s)) {
      continue;
    }
    if (check_prop_wsrc(job_tag, traj, type) and
        check_wall_src_info(job_tag, traj, type)) {
      check_sigint();
      check_time_limit();
      if (not obtain_lock(ssprintf("lock-wall-src-prop-norm-raio-%s-%d-%d",
                                   job_tag.c_str(), traj, type))) {
        continue;
      }
      setup(job_tag, traj);
      TIMER_VERBOSE("compute_wall_src_prop_norm_ratio");
      qmkdir_info("analysis");
      qmkdir_info("analysis/wall-src-prop-norm-ratio");
      qmkdir_info(
          ssprintf("analysis/wall-src-prop-norm-ratio/%s", job_tag.c_str()));
      qmkdir_info(rpath);
      const std::vector<WallInfo>& wis = get_wall_src_info(job_tag, traj, type);
      LatData ld_e, ld_s;
      int count = 0;
      for (int i = 0; i < (int)wis.size(); ++i) {
        const WallInfo& wi = wis[i];
        const int tslice = wi.tslice;
        qassert(wi.type == type);
        if (wi.accuracy == 2) {
          count += 1;
          const SelProp& s_prop_e =
              get_prop_wsrc(job_tag, traj, tslice, type, 2);
          const SelProp& s_prop_s =
              get_prop_wsrc(job_tag, traj, tslice, type, 1);
          const FieldSelection& fsel = get_field_selection(job_tag, traj);
          ld_e += contract_pion(s_prop_e, fsel, tslice);
          ld_s += contract_pion(s_prop_s, fsel, tslice);
        }
      }
      ld_e *= 1.0 / (double)count;
      ld_s *= 1.0 / (double)count;
      if (get_id_node() == 0) {
        ld_e.save(rpath_e);
        ld_s.save(rpath_s);
      }
      release_lock();
    } else {
      is_complete = false;
    }
  }
  if (is_complete) {
    qtouch_info(rpath + "/checkpoint.txt");
  }
}

}  // namespace qlat
