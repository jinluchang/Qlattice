#pragma once

#include "data-load.h"

namespace qlat
{  //

inline void show_pion_corr(const PselProp& ps_prop, const SelProp& s_prop,
                           const std::string& job_tag, const int traj,
                           const int tslice_src)
{
  TIMER_VERBOSE("show_pion_corr");
  qassert(is_initialized(ps_prop));
  qassert(is_initialized(s_prop));
  const PointSelection& psel = get_point_selection(job_tag, traj);
  const FieldSelection& fsel = get_field_selection(job_tag, traj);
  const Geometry& geo = fsel.f_rank.geo();
  const LatData ld_ps = contract_pion(ps_prop, tslice_src, geo, psel);
  const LatData ld_s = contract_pion(s_prop, tslice_src, fsel);
  displayln_info(fname + ssprintf(": PselProp pion corr"));
  display_info(show(ld_ps));
  displayln_info(fname + ssprintf(": SelProp pion corr"));
  display_info(show(ld_s));
}

inline void check_all_prop_psrc_exact(const std::string& job_tag,
                                      const int traj)
{
  if (check_prop_psrc_exact(job_tag, traj)) {
    setup(job_tag, traj);
    TIMER_VERBOSE("check_all_prop_psrc_exact");
    //
    qassert(check_sparse_parameters(job_tag, traj));
    qassert(check_point_src_info(job_tag, traj));
    //
    const std::vector<PointInfo>& pis = get_point_src_info(job_tag, traj);
    //
    qassert(not does_file_exist_sync_node(
        get_prop_psrc_exact_path(job_tag, traj), "CHECK-FILE"));
    //
    for (int i = 0; i < (int)pis.size(); ++i) {
      check_sigterm();
      const PointInfo& pi = pis[i];
      const Coordinate& xg = pi.xg;
      const int type = pi.type;
      const int accuracy = pi.accuracy;
      if (accuracy == 2) {
        qassert(get_does_prop_psrc_exact_exist(job_tag, traj, xg, type));
        const PselProp& ps_prop =
            get_psel_prop_psrc_exact(job_tag, traj, xg, type);
        const SelProp& s_prop = get_prop_psrc_exact(job_tag, traj, xg, type);
        show_pion_corr(ps_prop, s_prop, job_tag, traj, xg[3]);
      }
    }
  }
}

inline bool check_all_prop_psrc(const std::string& job_tag, const int traj)
{
  TIMER_VERBOSE("check_all_prop_psrc");
  const std::vector<PointInfo>& pis = get_point_src_info(job_tag, traj);
  // qassert(not does_file_exist_sync_node(get_prop_psrc_path(job_tag, traj, 0),
  //                                       "CHECK-FILE"));
  // qassert(not does_file_exist_sync_node(get_prop_psrc_path(job_tag, traj, 1),
  //                                       "CHECK-FILE"));
  for (int i = 0; i < (int)pis.size(); ++i) {
    check_sigterm();
    const PointInfo& pi = pis[i];
    const Coordinate& xg = pi.xg;
    const int type = pi.type;
    const int accuracy = pi.accuracy;
    if (not get_does_prop_psrc_exist(job_tag, traj, xg, type, accuracy)) {
      return false;
    }
    if (accuracy == 2) {
      const PselProp& ps_prop =
          get_psel_prop_psrc(job_tag, traj, xg, type, accuracy);
      const SelProp& s_prop = get_prop_psrc(job_tag, traj, xg, type, accuracy);
      show_pion_corr(ps_prop, s_prop, job_tag, traj, xg[3]);
    }
  }
  return true;
}

inline bool check_all_prop_wsrc(const std::string& job_tag, const int traj,
                                const int type)
{
  TIMER_VERBOSE("check_all_prop_wsrc");
  // qassert(not does_file_exist_sync_node(get_prop_wsrc_path(job_tag, traj,
  // type),
  //                                       "CHECK-FILE"));
  const Coordinate total_site = get_total_site(job_tag);
  const std::string path = get_prop_wsrc_path(job_tag, traj, type);
  const std::string path_psel = get_psel_prop_wsrc_path(job_tag, traj, type);
  long count_exact = 0;
  for (int tslice = 0; tslice < total_site[3]; ++tslice) {
    check_sigterm();
    if (not get_does_prop_wsrc_exist(job_tag, traj, tslice, type, 1)) {
      return false;
    }
    if (get_does_prop_wsrc_exist(job_tag, traj, tslice, type, 2)) {
      count_exact += 1;
      const PselProp& ps_prop_sloppy =
          get_psel_prop_wsrc(job_tag, traj, tslice, type, 1);
      const SelProp& s_prop_sloppy =
          get_prop_wsrc(job_tag, traj, tslice, type, 1);
      show_pion_corr(ps_prop_sloppy, s_prop_sloppy, job_tag, traj, tslice);
      const PselProp& ps_prop_exact =
          get_psel_prop_wsrc(job_tag, traj, tslice, type, 2);
      const SelProp& s_prop_exact =
          get_prop_wsrc(job_tag, traj, tslice, type, 2);
      show_pion_corr(ps_prop_exact, s_prop_exact, job_tag, traj, tslice);
    }
  }
  return 1 <= count_exact and count_exact <= 2;
}

inline void check_prop_data(const std::string& job_tag, const int traj)
{
  TIMER_VERBOSE("check_prop_data");
  if (check_prop_psrc(job_tag, traj, 0) and check_prop_psrc(job_tag, traj, 1)) {
    setup(job_tag, traj);
    //
    qassert(check_sparse_parameters(job_tag, traj));
    qassert(check_point_src_info(job_tag, traj));
    qassert(check_point_distribution(job_tag));
    //
    const PointDistribution& pd = get_point_distribution(job_tag);
    const Coordinate total_site = get_total_site(job_tag);
    qassert(weight_from_pd(pd, total_site, Coordinate(), Coordinate()) == 1.0);
    //
    qassert(check_all_prop_psrc(job_tag, traj));
  }
  if (check_prop_wsrc(job_tag, traj, 0)) {
    setup(job_tag, traj);
    //
    qassert(check_sparse_parameters(job_tag, traj));
    qassert(check_gauge_transform(job_tag, traj));
    //
    qassert(check_all_prop_wsrc(job_tag, traj, 0));
  }
  if (check_prop_wsrc(job_tag, traj, 1)) {
    setup(job_tag, traj);
    //
    qassert(check_sparse_parameters(job_tag, traj));
    qassert(check_gauge_transform(job_tag, traj));
    //
    qassert(check_all_prop_wsrc(job_tag, traj, 1));
  }
}

}  // namespace qlat
