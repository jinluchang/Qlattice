#pragma once

#include "data-load-base.h"
#include "data-load-wsrc-info.h"

namespace qlat
{  //

typedef Cache<std::string, WallSrcProps> WallSrcPropsCache;

inline WallSrcPropsCache& get_wall_src_props_cache()
{
  static WallSrcPropsCache cache("WallSrcPropsCache", 6, 2);
  return cache;
}

inline const WallSrcProps& get_wall_src_props(const std::string& job_tag,
                                              const int traj, const int type)
{
  const std::string key = ssprintf("%s,%d,%d,wsp", job_tag.c_str(), traj, type);
  WallSrcPropsCache& cache = get_wall_src_props_cache();
  if (not cache.has(key)) {
    TIMER_VERBOSE("get_wall_src_props");
    display_fields_wsrc(job_tag, traj, type);
    WallSrcProps& wsp = cache[key];
    const Coordinate total_site = get_total_site(job_tag);
    wsp.initialized = true;
    wsp.sloppy.resize(total_site[3]);
    wsp.exact.resize(total_site[3]);
    wsp.sloppy_point_snk.resize(total_site[3]);
    wsp.exact_point_snk.resize(total_site[3]);
    const std::vector<WallInfo>& wis = get_wall_src_info(job_tag, traj, type);
    for (int i = 0; i < (int)wis.size(); ++i) {
      const WallInfo& wi = wis[i];
      if (wi.accuracy == 1) {
        wsp.sloppy[wi.tslice] =
            get_prop_wsrc(job_tag, traj, wi.tslice, wi.type, wi.accuracy);
        wsp.sloppy_point_snk[wi.tslice] =
            get_psel_prop_wsrc(job_tag, traj, wi.tslice, wi.type, wi.accuracy);
      } else if (wi.accuracy == 2) {
        wsp.exact[wi.tslice] =
            get_prop_wsrc(job_tag, traj, wi.tslice, wi.type, wi.accuracy);
        wsp.exact_point_snk[wi.tslice] =
            get_psel_prop_wsrc(job_tag, traj, wi.tslice, wi.type, wi.accuracy);
      } else {
        qassert(false);
      }
    }
    const int num_exact = 2;
    const GaugeTransform& gt = get_gauge_transform(job_tag, traj);
    const PointSelection& psel = get_point_selection(job_tag, traj);
    const FieldSelection& fsel = get_field_selection(job_tag, traj);
    refresh_wsp(wsp, num_exact, gt, psel, fsel);
  }
  return cache[key];
}

inline bool check_wall_src_props(const std::string& job_tag, const int traj,
                                 const int type)
{
  TIMER_VERBOSE("check_wall_src_props");
  return check_prop_wsrc(job_tag, traj, type) and
         check_gauge_transform(job_tag, traj) and
         check_sparse_parameters(job_tag, traj);
}

}  // namespace qlat
