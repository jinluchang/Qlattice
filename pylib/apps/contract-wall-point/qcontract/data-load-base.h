#pragma once

#include <qlat/qlat.h>

#include "configs.h"
#include "data-paths.h"
#include "psrc-distribution.h"
#include "psrc-sample.h"
#include "qlat-setup.h"

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

typedef Cache<std::string, ShuffledBitSet> ShuffledBitSetCache;

typedef Cache<std::string, std::vector<PointInfo> > PointSrcInfoCache;

typedef Cache<std::string, PointDistribution> PointDistributionCache;

typedef Cache<std::string, GaugeTransform> GaugeTransformCache;

typedef Cache<std::string, PselProp> PselPropCache;

typedef Cache<std::string, SelProp> SelPropCache;

typedef std::vector<std::vector<int> > TypeAccuracyTable;

typedef Cache<std::string, TypeAccuracyTable> TypeAccuracyTableCache;

typedef std::map<std::string, std::vector<PointInfo> > PointInfoMap;

typedef Cache<std::string, PointInfoMap> PointInfoMapCache;

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

inline ShuffledBitSetCache& get_shuffled_bit_set_cache()
{
  static ShuffledBitSetCache cache("ShuffledBitSetCache", 8, 2);
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
  static PselPropCache cache("PselPropCache", 128, 4);
  return cache;
}

inline SelPropCache& get_prop_psrc_cache()
{
  static SelPropCache cache("PropPsrcCache", 128, 2);
  return cache;
}

inline SelPropCache& get_prop_wsrc_cache()
{
  static SelPropCache cache("PropWsrcCache", 4, 1);
  return cache;
}

inline TypeAccuracyTableCache& get_type_accuracy_table_cache()
{
  static TypeAccuracyTableCache cache("TypeAccuracyTableCache", 4, 2);
  return cache;
}

inline PointInfoMapCache& get_point_info_map_cache()
{
  static PointInfoMapCache cache("PointInfoMapCache", 4, 2);
  return cache;
}

inline const PointSelection& get_point_selection(const std::string& job_tag,
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

inline const FieldSelection& get_field_selection(const std::string& job_tag,
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

inline const ShuffledBitSet& get_shuffled_bit_set(const std::string& job_tag,
                                                  const int traj,
                                                  const std::string& path)
{
  const std::string key =
      ssprintf("%s,%d,%s,sbs", job_tag.c_str(), traj, path.c_str());
  ShuffledBitSetCache& cache = get_shuffled_bit_set_cache();
  if (not cache.has(key)) {
    TIMER_VERBOSE("get_shuffled_bit_set");
    const FieldSelection& fsel = get_field_selection(job_tag, traj);
    const PointSelection& psel = get_point_selection(job_tag, traj);
    FieldSelection fselc;
    fselc.f_rank = fsel.f_rank;
    add_field_selection(fselc.f_rank, psel);
    update_field_selection(fselc);
    update_field_selection(fselc, fsel.n_per_tslice);
    ShuffledFieldsReader& sfr = get_shuffled_fields_reader(path);
    cache[key] = mk_shuffled_bitset(fselc, sfr.new_size_node);
  }
  return cache[key];
}

inline const std::vector<PointInfo>& get_point_src_info(
    const std::string& job_tag, const int traj)
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

inline const PointDistribution& get_point_distribution(
    const std::string& job_tag)
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

inline const GaugeTransform& get_gauge_transform(const std::string& job_tag,
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

inline const TypeAccuracyTable& get_type_accuracy_table(
    const std::string& job_tag, const int traj)
{
  const std::string key = ssprintf("%s,%d,tat", job_tag.c_str(), traj);
  TypeAccuracyTableCache& cache = get_type_accuracy_table_cache();
  if (not cache.has(key)) {
    TIMER_VERBOSE("get_type_accuracy_table");
    cache[key] = mk_type_accuracy_table(get_point_src_info(job_tag, traj));
  }
  return cache[key];
}

inline const PointInfoMap& get_point_info_map(const std::string& job_tag,
                                              const int traj)
{
  const std::string key = ssprintf("%s,%d,pim", job_tag.c_str(), traj);
  PointInfoMapCache& cache = get_point_info_map_cache();
  if (not cache.has(key)) {
    TIMER_VERBOSE("get_point_info_map");
    PointInfoMap& pim = cache[key];
    const std::vector<PointInfo>& pis = get_point_src_info(job_tag, traj);
    for (long i = 0; i < (long)pis.size(); ++i) {
      const PointInfo& pi = pis[i];
      const std::string tag = ssprintf("%s,%d", show(pi.xg).c_str(), pi.type);
      pim[tag].push_back(pi);
    }
  }
  return cache[key];
}

inline const std::vector<PointInfo>& get_point_src_info(
    const std::string& job_tag, const int traj, const Coordinate& xg,
    const int type)
{
  const PointInfoMap& pim = get_point_info_map(job_tag, traj);
  const std::string tag = ssprintf("%s,%d", show(xg).c_str(), type);
  if (pim.count(tag) > 0) {
    return pim.at(tag);
  } else {
    static std::vector<PointInfo> pis_empty;
    return pis_empty;
  }
}

inline void display_fields_psrc(const std::string& job_tag, const int traj,
                                const int type)
{
  TIMER_VERBOSE("display_fields_psrc");
  const std::string path = get_prop_psrc_path(job_tag, traj, type);
  const std::vector<std::string> fns = list_fields(path);
  displayln_info(fname + ssprintf(": path='%s'.", path.c_str()));
  display_info(show_list(fns));
}

inline void display_fields_wsrc(const std::string& job_tag, const int traj,
                                const int type)
{
  TIMER_VERBOSE("display_fields_wsrc");
  const std::string path = get_prop_wsrc_path(job_tag, traj, type);
  const std::vector<std::string> fns = list_fields(path);
  displayln_info(fname + ssprintf(": path='%s'.", path.c_str()));
  display_info(show_list(fns));
}

inline long load_prop(PselProp& ps_prop, SelProp& s_prop,
                      const std::string& path, const std::string& fn,
                      const PointSelection& psel, const FieldSelection& fsel)
{
  TIMER_VERBOSE("load_prop(ps_prop,s_prop,path,fn,psel,fsel)");
  displayln_info(
      fname +
      ssprintf(": WARNING: obsolete. Try to use the sbs version instead."));
  Propagator4d prop;
  const long total_bytes = read_field_double_from_float(prop, path, fn);
  if (total_bytes > 0) {
    set_selected_points(ps_prop, prop, psel);
    set_selected_field(s_prop, prop, fsel);
  }
  return total_bytes;
}

inline long load_prop(PselProp& ps_prop, SelProp& s_prop,
                      const std::string& path, const std::string& fn,
                      const PointSelection& psel, const FieldSelection& fsel,
                      const ShuffledBitSet& sbs)
{
  TIMER_VERBOSE("load_prop(ps_prop,s_prop,path,fn,psel,fsel,sbs)");
  SelProp sprop;
  const long total_bytes = read_field_double_from_float(sprop, path, fn, sbs);
  if (total_bytes > 0) {
    set_selected_points(ps_prop, sprop, psel, sbs.fsel);
    set_selected_field(s_prop, sprop, fsel, sbs.fsel);
  } else {
    qassert(false);
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
    qassert(false);
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

inline const PselProp& get_psel_prop_psrc(const std::string& job_tag,
                                          const int traj, const Coordinate& xg,
                                          const int type, const int accuracy)
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

inline const SelProp& get_prop_psrc(const std::string& job_tag, const int traj,
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
    const ShuffledBitSet& sbs = get_shuffled_bit_set(job_tag, traj, path);
    load_prop(ps_prop, s_prop, path, fn, psel, fsel, sbs);
    if (is_check_prop_consistency()) {
      TIMER_VERBOSE("check_prop_consistency");
      const std::string path_psel =
          get_psel_prop_psrc_path(job_tag, traj, type);
      PselProp ps_prop_load;
      load_prop(ps_prop_load, path_psel, fn);
      const double qnorm_ps_prop = qnorm(ps_prop_load);
      ps_prop_load -= ps_prop;
      const double qnorm_ps_prop_diff = qnorm(ps_prop_load);
      displayln_info(fname + ssprintf(": ps_prop diff qnorm = %24.17E. ps_prop "
                                      "qnorm = %24.17E or %24.17E.",
                                      qnorm_ps_prop_diff, qnorm_ps_prop,
                                      qnorm(ps_prop)));
      qassert(qnorm_ps_prop_diff <= 1e-10 * qnorm_ps_prop);
    }
  }
  return s_cache[key];
}

inline const PselProp& get_psel_prop_wsrc(const std::string& job_tag,
                                          const int traj, const int tslice,
                                          const int type, const int accuracy)
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

inline const SelProp& get_prop_wsrc(const std::string& job_tag, const int traj,
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
    const ShuffledBitSet& sbs = get_shuffled_bit_set(job_tag, traj, path);
    load_prop(ps_prop, s_prop, path, fn, psel, fsel, sbs);
    if (is_check_prop_consistency()) {
      TIMER_VERBOSE("check_prop_consistency");
      const std::string path_psel =
          get_psel_prop_wsrc_path(job_tag, traj, type);
      PselProp ps_prop_load;
      load_prop(ps_prop_load, path_psel, fn);
      const double qnorm_ps_prop = qnorm(ps_prop_load);
      ps_prop_load -= ps_prop;
      const double qnorm_ps_prop_diff = qnorm(ps_prop_load);
      displayln_info(fname + ssprintf(": ps_prop diff qnorm = %24.17E. ps_prop "
                                      "qnorm = %24.17E or %24.17E.",
                                      qnorm_ps_prop_diff, qnorm_ps_prop,
                                      qnorm(ps_prop)));
      qassert(qnorm_ps_prop_diff <= 1e-10 * qnorm_ps_prop);
    }
  }
  return s_cache[key];
}

inline const PselProp& get_psel_prop_psrc_exact(const std::string& job_tag,
                                                const int traj,
                                                const Coordinate& xg,
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

inline const SelProp& get_prop_psrc_exact(const std::string& job_tag,
                                          const int traj, const Coordinate& xg,
                                          const int type,
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
    const ShuffledBitSet& sbs = get_shuffled_bit_set(job_tag, traj, path);
    load_prop(ps_prop, s_prop, path, fn, psel, fsel, sbs);
    if (is_check_prop_consistency()) {
      TIMER_VERBOSE("check_prop_consistency");
      const std::string path_psel =
          get_psel_prop_psrc_exact_path(job_tag, traj);
      PselProp ps_prop_load;
      load_prop(ps_prop_load, path_psel, fn);
      const double qnorm_ps_prop = qnorm(ps_prop_load);
      ps_prop_load -= ps_prop;
      const double qnorm_ps_prop_diff = qnorm(ps_prop_load);
      displayln_info(fname + ssprintf(": ps_prop diff qnorm = %24.17E. ps_prop "
                                      "qnorm = %24.17E or %24.17E.",
                                      qnorm_ps_prop_diff, qnorm_ps_prop,
                                      qnorm(ps_prop)));
      qassert(qnorm_ps_prop_diff <= 1e-10 * qnorm_ps_prop);
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
  // ADJUST ME
  if (job_tag == "48I" and type == 1) {
    return false;
  }
  //
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
