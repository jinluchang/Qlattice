#pragma once

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
      check_sigterm();
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
          ld_e += contract_pion(s_prop_e, tslice, fsel);
          ld_s += contract_pion(s_prop_s, tslice, fsel);
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

inline LatData lat_data_cratio(const LatData& ld1, const LatData& ld2)
{
  LatData ld_ratio = ld1;
  Vector<Complex> ld_ratio_v = lat_data_cget(ld_ratio);
  const Vector<Complex> ld1_v = lat_data_cget_const(ld1);
  const Vector<Complex> ld2_v = lat_data_cget_const(ld2);
  qassert(ld1_v.size() == ld2_v.size());
  for (long i = 0; i < (long)ld_ratio_v.size(); ++i) {
    const Complex c1 = ld1_v[i];
    const Complex c2 = ld2_v[i];
    ld_ratio_v[i] = c1 / c2;
  }
  return ld_ratio;
}

inline LatData lat_data_cratio_sqrt(const LatData& ld1, const LatData& ld2)
{
  LatData ld_ratio = ld1;
  Vector<Complex> ld_ratio_v = lat_data_cget(ld_ratio);
  const Vector<Complex> ld1_v = lat_data_cget_const(ld1);
  const Vector<Complex> ld2_v = lat_data_cget_const(ld2);
  qassert(ld1_v.size() == ld2_v.size());
  for (long i = 0; i < (long)ld_ratio_v.size(); ++i) {
    const Complex c1 = ld1_v[i];
    const Complex c2 = ld2_v[i];
    ld_ratio_v[i] = std::sqrt(c1 / c2);
  }
  return ld_ratio;
}

inline LatData get_wall_src_prop_norm_ratio_compute(const std::string& job_tag,
                                                    const int type)
// return the time (relative to the source) dependence sqrt of the ratio of the
// sqrt(norm_of_sloppy / norm_of_exact)
{
  TIMER_VERBOSE("get_wall_src_prop_norm_ratio_compute");
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
      if (get_id_node() == 0) {
        ld_e.load(rpath_e);
        ld_s.load(rpath_s);
      }
      bcast(ld_e);
      bcast(ld_s);
      ld_acc_e += ld_e;
      ld_acc_s += ld_s;
    }
  }
  const Coordinate total_site = get_total_site(job_tag);
  LatData ld_ratio = ld_acc_s;
  Vector<Complex> ld_ratio_v = lat_data_cget(ld_ratio);
  const Vector<Complex> ld_acc_e_v = lat_data_cget_const(ld_acc_e);
  const Vector<Complex> ld_acc_s_v = lat_data_cget_const(ld_acc_s);
  qassert(ld_ratio_v.size() == total_site[3]);
  qassert(ld_acc_e_v.size() == total_site[3]);
  qassert(ld_acc_s_v.size() == total_site[3]);
  for (long i = 0; i < (long)ld_ratio_v.size(); ++i) {
    const Complex s = ld_acc_s_v[i] + ld_acc_s_v[mod(-i, total_site[3])];
    const Complex e = ld_acc_e_v[i] + ld_acc_e_v[mod(-i, total_site[3])];
    ld_ratio_v[i] = std::sqrt(s / e);
  }
  const std::string ratio_path =
      ssprintf("analysis/wall-src-prop-norm-ratio/%s/ratio-%s.lat",
               job_tag.c_str(), type_tag.c_str());
  ld_ratio.save(ratio_path);
  return ld_ratio;
}

typedef Cache<std::string, LatData> WallSrcPropNormRatioCache;

inline WallSrcPropNormRatioCache& get_wall_src_prop_norm_ratio_cache()
{
  static WallSrcPropNormRatioCache cache("WallSrcPropNormRatioCache", 8, 2);
  return cache;
}

inline LatData& get_wall_src_prop_norm_ratio(const std::string& job_tag,
                                             const int type)
{
  const std::string key = ssprintf("%s,%d,wspnr", job_tag.c_str(), type);
  WallSrcPropNormRatioCache& cache = get_wall_src_prop_norm_ratio_cache();
  if (not cache.has(key)) {
    TIMER_VERBOSE("get_wall_src_prop_norm_ratio");
    cache[key] = get_wall_src_prop_norm_ratio_compute(job_tag, type);
  }
  return cache[key];
}

inline void scale_prop_wsrc_with_ratio(SelProp& s_prop, const int tslice_src,
                                       const LatData& ld,
                                       const FieldSelection& fsel)
{
  TIMER_VERBOSE("s_prop,ld");
  const Geometry& geo = s_prop.geo();
  const Coordinate total_site = geo.total_site();
  const Vector<Complex> ldv = lat_data_cget_const(ld);
#pragma omp parallel for
  for (long idx = 0; idx < fsel.n_elems; ++idx) {
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const int tsep = mod(xg[3] - tslice_src, total_site[3]);
    s_prop.get_elem(idx) *= 1.0 / ldv[tsep];
  }
}

inline void test_wall_src_norm_ratio()
{
  TIMER_VERBOSE("test_wall_src_norm_ratio");
  // const std::string job_tag = "48I";
  // const int traj = 970;
  // const int tslice = 10;
  // const std::string job_tag = "32Dfine";
  // const int traj = 1000;
  // const int tslice = 21;
  const std::string job_tag = "24D";
  const int traj = 2500;
  const int tslice = 33;
  const LatData& ld_light = get_wall_src_prop_norm_ratio(job_tag, 0);
  displayln_info("48I light ratio");
  displayln_info(show(ld_light));
  const LatData& ld_strange = get_wall_src_prop_norm_ratio(job_tag, 1);
  displayln_info("48I strange ratio");
  displayln_info(show(ld_strange));
  const LatData& ld_scale = ld_strange;
  const FieldSelection& fsel = get_field_selection(job_tag, traj);
  const SelProp& s_prop_l_e = get_prop_wsrc(job_tag, traj, tslice, 0, 2);
  const SelProp& s_prop_l_s = get_prop_wsrc(job_tag, traj, tslice, 0, 1);
  const SelProp& s_prop_e = get_prop_wsrc(job_tag, traj, tslice, 1, 2);
  const SelProp& s_prop_s = get_prop_wsrc(job_tag, traj, tslice, 1, 1);
  SelProp s_prop_s_scaled;
  s_prop_s_scaled = s_prop_s;
  scale_prop_wsrc_with_ratio(s_prop_s_scaled, tslice, ld_scale, fsel);
  // s_prop_s_scaled = s_prop_l_s;
  // scale_prop_wsrc_with_ratio(
  //     s_prop_s_scaled, tslice,
  //     lat_data_cratio_sqrt(contract_pion(s_prop_l_s, tslice, fsel),
  //                          contract_pion(s_prop_e, tslice, fsel)), fsel);
  SelProp s_prop_diff;
  s_prop_diff = s_prop_s_scaled;
  s_prop_diff -= s_prop_e;
  {
    displayln_info("48I pion exact");
    displayln_info(show(contract_pion(s_prop_e, tslice, fsel)));
    displayln_info("48I pion sloppy ratio");
    displayln_info(
        show(lat_data_cratio(contract_pion(s_prop_s, tslice, fsel),
                             contract_pion(s_prop_e, tslice, fsel))));
    displayln_info("48I pion sloppy scaled ratio");
    displayln_info(
        show(lat_data_cratio(contract_pion(s_prop_s_scaled, tslice, fsel),
                             contract_pion(s_prop_e, tslice, fsel))));
    displayln_info("48I pion diff ratio");
    displayln_info(
        show(lat_data_cratio(contract_pion(s_prop_diff, tslice, fsel),
                             contract_pion(s_prop_e, tslice, fsel))));
    displayln_info("48I pion wall exact");
    displayln_info(show(contract_pion_wall_snk(s_prop_e, tslice, fsel)));
    displayln_info("48I pion wall sloppy ratio");
    displayln_info(
        show(lat_data_cratio(contract_pion_wall_snk(s_prop_s, tslice, fsel),
                             contract_pion_wall_snk(s_prop_e, tslice, fsel))));
    displayln_info("48I pion wall sloppy scaled ratio");
    displayln_info(show(
        lat_data_cratio(contract_pion_wall_snk(s_prop_s_scaled, tslice, fsel),
                        contract_pion_wall_snk(s_prop_e, tslice, fsel))));
    displayln_info("48I pion wall diff ratio");
    displayln_info(
        show(lat_data_cratio(contract_pion_wall_snk(s_prop_diff, tslice, fsel),
                             contract_pion_wall_snk(s_prop_e, tslice, fsel))));
  }
  {
    displayln_info("48I kaon exact");
    displayln_info(show(contract_kaon(s_prop_l_e, s_prop_e, tslice, fsel)));
    displayln_info("48I kaon sloppy ratio");
    displayln_info(show(
        lat_data_cratio(contract_kaon(s_prop_l_s, s_prop_s, tslice, fsel),
                        contract_kaon(s_prop_l_e, s_prop_e, tslice, fsel))));
    displayln_info("48I kaon sloppy scaled ratio");
    displayln_info(show(lat_data_cratio(
        contract_kaon(s_prop_l_s, s_prop_s_scaled, tslice, fsel),
        contract_kaon(s_prop_l_e, s_prop_e, tslice, fsel))));
    displayln_info("48I kaon diff ratio");
    displayln_info(show(
        lat_data_cratio(contract_kaon(s_prop_l_s, s_prop_diff, tslice, fsel),
                        contract_kaon(s_prop_l_e, s_prop_e, tslice, fsel))));
    displayln_info("48I kaon wall exact");
    displayln_info(
        show(contract_kaon_wall_snk(s_prop_l_e, s_prop_e, tslice, fsel)));
    displayln_info("48I kaon wall sloppy ratio");
    displayln_info(show(lat_data_cratio(
        contract_kaon_wall_snk(s_prop_l_s, s_prop_s, tslice, fsel),
        contract_kaon_wall_snk(s_prop_l_e, s_prop_e, tslice, fsel))));
    displayln_info("48I kaon wall sloppy scaled ratio");
    displayln_info(show(lat_data_cratio(
        contract_kaon_wall_snk(s_prop_l_s, s_prop_s_scaled, tslice, fsel),
        contract_kaon_wall_snk(s_prop_l_e, s_prop_e, tslice, fsel))));
    displayln_info("48I kaon wall diff ratio");
    displayln_info(show(lat_data_cratio(
        contract_kaon_wall_snk(s_prop_l_s, s_prop_diff, tslice, fsel),
        contract_kaon_wall_snk(s_prop_l_e, s_prop_e, tslice, fsel))));
  }
  exit(0);
}

}  // namespace qlat
