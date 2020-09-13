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

inline LatData get_wall_src_prop_norm_ratio(const std::string& job_tag,
                                            const int type)
// return the time (relative to the source) dependence sqrt of the ratio of the
// sqrt(norm_of_sloppy / norm_of_exact)
{
  TIMER_VERBOSE("get_wall_src_prop_norm_ratio");
  const std::string type_tag = type == 0 ? "light" : "strange";
  LatData ld_acc_e, ld_acc_s;
  const std::vector<int> trajs = get_data_trajs(job_tag);
  for (int i = 0; i < (int)trajs.size(); ++i) {
    const int traj = trajs[i];
    const std::string rpath = get_wall_src_prop_norm_ratio_path(job_tag, traj);
    const std::string rpath_e =
        rpath + ssprintf("/norm-exact-%s.lat", type_tag.c_str());
    const std::string rpath_s =
        rpath + ssprintf("/norm-sloppy-%s.lat", type_tag.c_str());
    if (get_does_file_exist(rpath_e) and get_does_file_exist(rpath_s)) {
      LatData ld_e, ld_s;
      ld_e.load(rpath_e);
      ld_s.load(rpath_s);
      ld_acc_e += ld_e;
      ld_acc_s += ld_s;
    }
  }
  const Coordinate total_site = get_total_site(job_tag);
  LatData ld_ratio = ld_acc_s;
  Vector<Complex> ld_ratio_v =
      lat_data_complex_get(ld_ratio, make_array<int>());
  const Vector<Complex> ld_acc_e_v =
      lat_data_complex_get(ld_acc_e, make_array<int>());
  const Vector<Complex> ld_acc_s_v =
      lat_data_complex_get(ld_acc_s, make_array<int>());
  qassert(ld_ratio_v.size() == total_site[3]);
  qassert(ld_acc_e_v.size() == total_site[3]);
  qassert(ld_acc_s_v.size() == total_site[3]);
  for (long i = 0; i < (long)ld_ratio_v.size(); ++i) {
    const Complex s = ld_acc_s_v[i] + ld_acc_s_v[mod(-i, total_site[3])];
    const Complex e = ld_acc_e_v[i] + ld_acc_e_v[mod(-i, total_site[3])];
    ld_ratio_v[i] = std::sqrt(s / e);
  }
  return ld_ratio;
}

}  // namespace qlat
