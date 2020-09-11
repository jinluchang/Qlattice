#pragma once

#include <qlat/qlat.h>

#include "data-paths.h"
#include "psrc-distribution.h"
#include "psrc-sample.h"

namespace qlat
{  //

typedef Cache<std::string, bool> DoesFileExistCache;

inline DoesFileExistCache& get_does_file_exist_cache()
{
  static DoesFileExistCache cache("DoesFileExistCache", 4096, 1024);
  return cache;
}

inline bool get_does_file_exist(const std::string& fn)
{
  const std::string key = "fn:" + fn;
  DoesFileExistCache& cache = get_does_file_exist_cache();
  if (not cache.has(key)) {
    TIMER_VERBOSE("get_does_file_exist");
    cache[key] = does_file_exist_sync_node(fn);
  }
  return cache[key];
}

typedef Cache<std::string, PointSelection> PointSelectionCache;

typedef Cache<std::string, FieldSelection> FieldSelectionCache;

typedef Cache<std::string, std::vector<PointInfo> > PointSrcInfoCache;

typedef Cache<std::string, PointDistribution> PointDistributionCache;

typedef Cache<std::string, GaugeTransform> GaugeTransformCache;

typedef Cache<std::string, PselProp> PselPropCache;

typedef Cache<std::string, SelProp> SelPropCache;

inline PointSelectionCache& get_point_selection_cache()
{
  static PointSelectionCache cache("PointSelectionCache", 8, 2);
  return cache;
}

inline FieldSelectionCache& get_field_selection_cache()
{
  static FieldSelectionCache cache("FieldSelectionCache", 8, 2);
  return cache;
}

inline PointSrcInfoCache& get_point_src_info_cache()
{
  static PointSrcInfoCache cache("PointSrcInfoCache", 8, 2);
  return cache;
}

inline PointDistributionCache& get_point_distribution_cache()
{
  static PointDistributionCache cache("PointDistributionCache", 4, 1);
  return cache;
}

inline GaugeTransformCache& get_gauge_transform_cache()
{
  static GaugeTransformCache cache("GaugeTransformCache", 8, 2);
  return cache;
}

inline PselPropCache& get_psel_prop_cache()
{
  static PselPropCache cache("PselPropCache", 2048, 128);
  return cache;
}

inline SelPropCache& get_prop_psrc_cache()
{
  static SelPropCache cache("PropPsrcCache", 128, 8);
  return cache;
}

inline SelPropCache& get_prop_wsrc_cache()
{
  static SelPropCache cache("PropWsrcCache", 256, 8);
  return cache;
}

inline void clear_all_prop_cache()
{
  TIMER_VERBOSE("clear_all_prop_cache");
  get_prop_psrc_cache().clear();
  get_prop_wsrc_cache().clear();
}

inline void clear_all_data_cache()
{
  TIMER_VERBOSE("clear_all_cache");
  get_point_selection_cache().clear();
  get_field_selection_cache().clear();
  get_point_src_info_cache().clear();
  get_gauge_transform_cache().clear();
  get_psel_prop_cache().clear();
  clear_all_prop_cache();
  // ADJUST ME
  // get_point_distribution_cache().clear();
  // get_does_file_exist_cache().clear();
  //
}

inline PointSelection& get_point_selection(const std::string& job_tag,
                                           const int traj)
{
  const std::string key = ssprintf("%s,%d,psel", job_tag.c_str(), traj);
  PointSelectionCache& cache = get_point_selection_cache();
  if (not cache.has(key)) {
    TIMER_VERBOSE("get_point_selection");
    const std::string fn_point_selection =
        get_point_selection_path(job_tag, traj);
    qassert(get_does_file_exist(fn_point_selection));
    cache[key] = load_point_selection_info(fn_point_selection);
  }
  return cache[key];
}

inline FieldSelection& get_field_selection(const std::string& job_tag,
                                           const int traj)
{
  const std::string key = ssprintf("%s,%d,fsel", job_tag.c_str(), traj);
  FieldSelectionCache& cache = get_field_selection_cache();
  if (not cache.has(key)) {
    TIMER_VERBOSE("get_field_selection");
    const std::string fn_field_selection =
        get_field_selection_path(job_tag, traj);
    qassert(get_does_file_exist(fn_field_selection));
    const Coordinate total_site = get_total_site(job_tag);
    const long spatial_vol = total_site[0] * total_site[1] * total_site[2];
    const long n_per_tslice = spatial_vol / 16;
    read_field_selection(cache[key], fn_field_selection, n_per_tslice);
  }
  return cache[key];
}

inline std::vector<PointInfo>& get_point_src_info(const std::string& job_tag,
                                                  const int traj)
{
  const std::string key = ssprintf("%s,%d,pis", job_tag.c_str(), traj);
  PointSrcInfoCache& cache = get_point_src_info_cache();
  if (not cache.has(key)) {
    TIMER_VERBOSE("get_point_src_info");
    const std::string fn = get_point_src_info_path(job_tag, traj);
    qassert(get_does_file_exist(fn));
    cache[key] = load_lbl_pis_info(fn);
  }
  return cache[key];
}

inline PointDistribution& get_point_distribution(const std::string& job_tag)
{
  const std::string key = ssprintf("%s,pd", job_tag.c_str());
  PointDistributionCache& cache = get_point_distribution_cache();
  if (not cache.has(key)) {
    TIMER_VERBOSE("get_point_src_info");
    const std::string fn = get_point_distribution_path(job_tag);
    const Coordinate total_site = get_total_site(job_tag);
    qassert(get_does_file_exist(fn));
    load_pd(cache[key], total_site, fn);
  }
  return cache[key];
}

inline GaugeTransform& get_gauge_transform(const std::string& job_tag,
                                           const int traj)
{
  const std::string key = ssprintf("%s,%d,gt", job_tag.c_str(), traj);
  GaugeTransformCache& cache = get_gauge_transform_cache();
  if (not cache.has(key)) {
    TIMER_VERBOSE("get_gauge_transform");
    const std::string fn = get_gauge_transform_path(job_tag, traj);
    qassert(get_does_file_exist(fn));
    read_field_double(cache[key], fn);
  }
  return cache[key];
}

inline long load_prop(PselProp& ps_prop, SelProp& s_prop,
                      const std::string& path, const std::string& fn,
                      const PointSelection& psel, const FieldSelection& fsel)
{
  TIMER_VERBOSE("load_prop(ps_prop,s_prop,path,fn,psel,fsel)");
  Propagator4d prop;
  const long total_bytes = read_field_double_from_float(prop, path, fn);
  if (total_bytes > 0) {
    set_selected_points(ps_prop, prop, psel);
    set_selected_field(s_prop, prop, fsel);
  }
  return total_bytes;
}

inline void load_prop(PselProp& ps_prop, const std::string& path_psel,
                      const std::string& fn)
{
  TIMER_VERBOSE("load_prop(ps_prop,path,fn)");
  const std::string path_full = path_psel + "/" + fn + ".lat";
  if (get_does_file_exist(path_full)) {
    load_selected_points_complex(ps_prop, path_full);
  } else {
    ps_prop.init();
  }
}

inline std::string get_prop_wsrc_key(const std::string& job_tag, const int traj,
                                     const int tslice, const int type,
                                     const int accuracy)
{
  return ssprintf("%d,%d,%d,%s,%d,wsrc", tslice, type, accuracy,
                  job_tag.c_str(), traj);
}

inline std::string get_prop_psrc_key(const std::string& job_tag, const int traj,
                                     const Coordinate& xg, const int type,
                                     const int accuracy)
{
  return ssprintf("%s,%d,%d,%s,%d,psrc", show(xg).c_str(), type, accuracy,
                  job_tag.c_str(), traj);
}

inline bool& is_check_prop_consistency()
{
  static bool b = true;
  return b;
}

inline PselProp& get_psel_prop_psrc(const std::string& job_tag, const int traj,
                                    const Coordinate& xg, const int type,
                                    const int accuracy)
{
  const std::string key = get_prop_psrc_key(job_tag, traj, xg, type, accuracy);
  PselPropCache& ps_cache = get_psel_prop_cache();
  if (not ps_cache.has(key)) {
    TIMER_VERBOSE("get_psel_prop_psrc");
    PselProp& ps_prop = ps_cache[key];
    const std::string path_psel = get_psel_prop_psrc_path(job_tag, traj, type);
    const std::string fn = get_psrc_tag(xg, type, accuracy);
    load_prop(ps_prop, path_psel, fn);
  }
  return ps_cache[key];
}

inline SelProp& get_prop_psrc(const std::string& job_tag, const int traj,
                              const Coordinate& xg, const int type,
                              const int accuracy)
{
  const std::string key = get_prop_psrc_key(job_tag, traj, xg, type, accuracy);
  SelPropCache& s_cache = get_prop_psrc_cache();
  if (not s_cache.has(key)) {
    TIMER_VERBOSE("get_prop_psrc");
    PselPropCache& ps_cache = get_psel_prop_cache();
    PselProp& ps_prop = ps_cache[key];
    SelProp& s_prop = s_cache[key];
    const std::string path = get_prop_psrc_path(job_tag, traj, type);
    const std::string fn = get_psrc_tag(xg, type, accuracy);
    const PointSelection& psel = get_point_selection(job_tag, traj);
    const FieldSelection& fsel = get_field_selection(job_tag, traj);
    load_prop(ps_prop, s_prop, path, fn, psel, fsel);
    if (is_check_prop_consistency()) {
      TIMER_VERBOSE("check_prop_consistency");
      const std::string path_psel =
          get_psel_prop_psrc_path(job_tag, traj, type);
      PselProp ps_prop_load;
      load_prop(ps_prop_load, path_psel, fn);
      const double qnorm_ps_prop = qnorm(ps_prop_load);
      ps_prop_load -= ps_prop;
      const double qnorm_ps_prop_diff = qnorm(ps_prop_load);
      displayln_info(fname +
                     ssprintf(": ps_prop diff qnorm = %lf. ps_prop qnorm = %lf",
                              qnorm_ps_prop_diff, qnorm_ps_prop));
      qassert(qnorm_ps_prop_diff == 0.0);
    }
  }
  return s_cache[key];
}

inline PselProp& get_psel_prop_wsrc(const std::string& job_tag, const int traj,
                                    const int tslice, const int type,
                                    const int accuracy)
{
  const std::string key =
      get_prop_wsrc_key(job_tag, traj, tslice, type, accuracy);
  PselPropCache& ps_cache = get_psel_prop_cache();
  if (not ps_cache.has(key)) {
    TIMER_VERBOSE("get_psel_prop_wsrc");
    PselProp& ps_prop = ps_cache[key];
    const std::string path_psel = get_psel_prop_wsrc_path(job_tag, traj, type);
    const std::string fn = get_wsrc_tag(tslice, type, accuracy);
    load_prop(ps_prop, path_psel, fn);
  }
  return ps_cache[key];
}

inline SelProp& get_prop_wsrc(const std::string& job_tag, const int traj,
                              const int tslice, const int type,
                              const int accuracy)
{
  const std::string key =
      get_prop_wsrc_key(job_tag, traj, tslice, type, accuracy);
  SelPropCache& s_cache = get_prop_wsrc_cache();
  if (not s_cache.has(key)) {
    TIMER_VERBOSE("get_prop_wsrc");
    PselPropCache& ps_cache = get_psel_prop_cache();
    PselProp& ps_prop = ps_cache[key];
    SelProp& s_prop = s_cache[key];
    const std::string path = get_prop_wsrc_path(job_tag, traj, type);
    const std::string fn = get_wsrc_tag(tslice, type, accuracy);
    const PointSelection& psel = get_point_selection(job_tag, traj);
    const FieldSelection& fsel = get_field_selection(job_tag, traj);
    load_prop(ps_prop, s_prop, path, fn, psel, fsel);
    if (is_check_prop_consistency()) {
      TIMER_VERBOSE("check_prop_consistency");
      const std::string path_psel =
          get_psel_prop_wsrc_path(job_tag, traj, type);
      PselProp ps_prop_load;
      load_prop(ps_prop_load, path_psel, fn);
      const double qnorm_ps_prop = qnorm(ps_prop_load);
      ps_prop_load -= ps_prop;
      const double qnorm_ps_prop_diff = qnorm(ps_prop_load);
      displayln_info(fname +
                     ssprintf(": ps_prop diff qnorm = %lf. ps_prop qnorm = %lf",
                              qnorm_ps_prop_diff, qnorm_ps_prop));
      qassert(qnorm_ps_prop_diff == 0.0);
    }
  }
  return s_cache[key];
}

inline PselProp& get_psel_prop_psrc_exact(const std::string& job_tag,
                                          const int traj, const Coordinate& xg,
                                          const int type,
                                          const int accuracy = 2)
{
  const std::string key =
      "exact:" + get_prop_psrc_key(job_tag, traj, xg, type, accuracy);
  PselPropCache& ps_cache = get_psel_prop_cache();
  if (not ps_cache.has(key)) {
    TIMER_VERBOSE("get_psel_prop_psrc");
    PselProp& ps_prop = ps_cache[key];
    const std::string path_psel = get_psel_prop_psrc_exact_path(job_tag, traj);
    const std::string fn = get_psrc_tag(xg, type, accuracy);
    load_prop(ps_prop, path_psel, fn);
  }
  return ps_cache[key];
}

inline SelProp& get_prop_psrc_exact(const std::string& job_tag, const int traj,
                                    const Coordinate& xg, const int type,
                                    const int accuracy = 2)
{
  const std::string key =
      "exact:" + get_prop_psrc_key(job_tag, traj, xg, type, accuracy);
  SelPropCache& s_cache = get_prop_psrc_cache();
  if (not s_cache.has(key)) {
    TIMER_VERBOSE("get_prop_psrc");
    PselPropCache& ps_cache = get_psel_prop_cache();
    PselProp& ps_prop = ps_cache[key];
    SelProp& s_prop = s_cache[key];
    const std::string path = get_prop_psrc_exact_path(job_tag, traj);
    const std::string fn = get_psrc_tag(xg, type, accuracy);
    const PointSelection& psel = get_point_selection(job_tag, traj);
    const FieldSelection& fsel = get_field_selection(job_tag, traj);
    load_prop(ps_prop, s_prop, path, fn, psel, fsel);
    if (is_check_prop_consistency()) {
      TIMER_VERBOSE("check_prop_consistency");
      const std::string path_psel =
          get_psel_prop_psrc_exact_path(job_tag, traj);
      PselProp ps_prop_load;
      load_prop(ps_prop_load, path_psel, fn);
      const double qnorm_ps_prop = qnorm(ps_prop_load);
      ps_prop_load -= ps_prop;
      const double qnorm_ps_prop_diff = qnorm(ps_prop_load);
      displayln_info(fname +
                     ssprintf(": ps_prop diff qnorm = %lf. ps_prop qnorm = %lf",
                              qnorm_ps_prop_diff, qnorm_ps_prop));
      qassert(qnorm_ps_prop_diff == 0.0);
    }
  }
  return s_cache[key];
}

inline bool get_does_psel_prop_psrc_exist(const std::string& job_tag,
                                          const int traj, const Coordinate& xg,
                                          const int type, const int accuracy)
{
  const std::string key =
      "psel-prop-psrc:" + get_prop_psrc_key(job_tag, traj, xg, type, accuracy);
  DoesFileExistCache& cache = get_does_file_exist_cache();
  if (not cache.has(key)) {
    TIMER_VERBOSE("get_does_psel_prop_psrc_exist");
    const std::string path_psel = get_psel_prop_psrc_path(job_tag, traj, type);
    const std::string fn = get_psrc_tag(xg, type, accuracy);
    const std::string path_full = path_psel + "/" + fn + ".lat";
    cache[key] = get_does_file_exist(path_full);
  }
  return cache[key];
}

inline bool get_does_prop_psrc_exist(const std::string& job_tag, const int traj,
                                     const Coordinate& xg, const int type,
                                     const int accuracy)
{
  const std::string key =
      "prop-psrc:" + get_prop_psrc_key(job_tag, traj, xg, type, accuracy);
  DoesFileExistCache& cache = get_does_file_exist_cache();
  if (not cache.has(key)) {
    TIMER_VERBOSE("get_does_prop_psrc_exist");
    const std::string path = get_prop_psrc_path(job_tag, traj, type);
    const std::string fn = get_psrc_tag(xg, type, accuracy);
    cache[key] = does_file_exist_sync_node(path, fn);
  }
  return cache[key];
}

inline bool get_does_psel_prop_wsrc_exist(const std::string& job_tag,
                                          const int traj, const int tslice,
                                          const int type, const int accuracy)
{
  const std::string key =
      "psel-prop-wsrc:" +
      get_prop_wsrc_key(job_tag, traj, tslice, type, accuracy);
  DoesFileExistCache& cache = get_does_file_exist_cache();
  if (not cache.has(key)) {
    TIMER_VERBOSE("get_does_psel_prop_wsrc_exist");
    const std::string path_psel = get_psel_prop_wsrc_path(job_tag, traj, type);
    const std::string fn = get_wsrc_tag(tslice, type, accuracy);
    const std::string path_full = path_psel + "/" + fn + ".lat";
    cache[key] = get_does_file_exist(path_full);
  }
  return cache[key];
}

inline bool get_does_prop_wsrc_exist(const std::string& job_tag, const int traj,
                                     const int tslice, const int type,
                                     const int accuracy)
{
  const std::string key =
      "prop-wsrc:" + get_prop_wsrc_key(job_tag, traj, tslice, type, accuracy);
  DoesFileExistCache& cache = get_does_file_exist_cache();
  if (not cache.has(key)) {
    TIMER_VERBOSE("get_does_prop_wsrc_exist");
    const std::string path = get_prop_wsrc_path(job_tag, traj, type);
    const std::string fn = get_wsrc_tag(tslice, type, accuracy);
    cache[key] = does_file_exist_sync_node(path, fn);
  }
  return cache[key];
}

inline bool get_does_psel_prop_psrc_exact_exist(const std::string& job_tag,
                                                const int traj,
                                                const Coordinate& xg,
                                                const int type,
                                                const int accuracy = 2)
{
  const std::string key = "psel-prop-psrc-exact:" +
                          get_prop_psrc_key(job_tag, traj, xg, type, accuracy);
  DoesFileExistCache& cache = get_does_file_exist_cache();
  if (not cache.has(key)) {
    TIMER_VERBOSE("get_does_psel_prop_psrc_exist");
    const std::string path_psel = get_psel_prop_psrc_exact_path(job_tag, traj);
    const std::string fn = get_psrc_tag(xg, type, accuracy);
    const std::string path_full = path_psel + "/" + fn + ".lat";
    cache[key] = get_does_file_exist(path_full);
  }
  return cache[key];
}

inline bool get_does_prop_psrc_exact_exist(const std::string& job_tag,
                                           const int traj, const Coordinate& xg,
                                           const int type,
                                           const int accuracy = 2)
{
  const std::string key =
      "prop-psrc-exact:" + get_prop_psrc_key(job_tag, traj, xg, type, accuracy);
  DoesFileExistCache& cache = get_does_file_exist_cache();
  if (not cache.has(key)) {
    TIMER_VERBOSE("get_does_prop_psrc_exist");
    const std::string path = get_prop_psrc_exact_path(job_tag, traj);
    const std::string fn = get_psrc_tag(xg, type, accuracy);
    cache[key] = does_file_exist_sync_node(path, fn);
  }
  return cache[key];
}

inline bool check_sparse_parameters(const std::string& job_tag, const int traj)
{
  TIMER_VERBOSE("check_sparse_parameters");
  const std::string fn_point_selection =
      get_point_selection_path(job_tag, traj);
  const std::string fn_field_selection =
      get_field_selection_path(job_tag, traj);
  return get_does_file_exist(fn_point_selection) and
         get_does_file_exist(fn_field_selection);
}

inline bool check_gauge_transform(const std::string& job_tag, const int traj)
{
  TIMER_VERBOSE("check_gauge_transform");
  return get_does_file_exist(get_gauge_transform_path(job_tag, traj));
}

inline bool check_point_distribution(const std::string& job_tag)
{
  TIMER_VERBOSE("check_point_distribution");
  return get_does_file_exist(get_point_distribution_path(job_tag));
}

inline bool check_point_src_info(const std::string& job_tag, const int traj)
{
  TIMER_VERBOSE("check_point_src_info");
  return get_does_file_exist(get_point_src_info_path(job_tag, traj));
}

inline bool check_prop_psrc(const std::string& job_tag, const int traj,
                            const int type)
{
  TIMER_VERBOSE("check_prop_psrc");
  return get_does_file_exist(get_prop_psrc_path(job_tag, traj, type)) and
         get_does_file_exist(get_psel_prop_psrc_path(job_tag, traj, type));
}

inline bool check_prop_wsrc(const std::string& job_tag, const int traj,
                            const int type)
{
  TIMER_VERBOSE("check_prop_wsrc");
  return get_does_file_exist(get_prop_wsrc_path(job_tag, traj, type)) and
         get_does_file_exist(get_psel_prop_wsrc_path(job_tag, traj, type));
}

inline bool check_prop_psrc_exact(const std::string& job_tag, const int traj)
{
  TIMER_VERBOSE("check_prop_psrc_exact");
  return get_does_file_exist(get_prop_psrc_exact_path(job_tag, traj)) and
         get_does_file_exist(get_psel_prop_psrc_exact_path(job_tag, traj));
}

}  // namespace qlat
