#pragma once

#include <qlat/qcd.h>
#include <qlat/selected-field.h>
#include <qlat/selected-points.h>

namespace qlat
{  //

typedef SelectedField<WilsonMatrix> SelProp;

typedef SelectedPoints<WilsonMatrix> PselProp;

struct WallSrcProps {
  bool initialized;
  std::vector<SelProp> sloppy;
  std::vector<SelProp> exact;
  std::vector<PselProp> sloppy_point_snk;
  std::vector<PselProp> exact_point_snk;
  //
  // below is refreshed by refresh_wall_snk_prop
  std::vector<std::vector<WilsonMatrix> >
      sloppy_wall_snk;  // need sparse correction when contracted with itself
  std::vector<std::vector<WilsonMatrix> >
      exact_wall_snk;  // need sparse correction when contracted with itself
  //
  // below is refreshed by refresh_prob
  std::vector<bool> exact_tslice_mask;  // true if has exact prop
  int num_sloppy;                       // number of sloppy props
  int num_exact;                // number of random choices (include duplicates)
  double sloppy_exact_ratio_1;  // inverse prob of one prop has exact version
  double sloppy_exact_ratio_11;  // inverse prob of two different props have
                                 // exact version
  //
  void init()
  {
    initialized = false;
    clear(sloppy);
    clear(exact);
    clear(sloppy_point_snk);
    clear(exact_point_snk);
    clear(sloppy_wall_snk);
    clear(exact_wall_snk);
    clear(exact_tslice_mask);
  }
  //
  WallSrcProps() { init(); }
};

inline bool is_initialized(const WallSrcProps& wsp) { return wsp.initialized; }

inline std::vector<WilsonMatrix> contract_wall_snk_prop(
    const SelProp& prop, const FieldSelection& fsel)
{
  TIMER_VERBOSE("contract_wall_snk_prop");
  const Geometry& geo = prop.geo;
  const Coordinate total_site = geo.total_site();
  std::vector<WilsonMatrix> gwm_ts(omp_get_max_threads() * total_site[3]);
  set_zero(gwm_ts);
#pragma omp parallel
  {
    std::vector<WilsonMatrix> wm_ts(total_site[3]);
    set_zero(wm_ts);
#pragma omp for
    for (long idx = 0; idx < (long)fsel.indices.size(); ++idx) {
      const long index = fsel.indices[idx];
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      wm_ts[xg[3]] += prop.get_elem(idx);
    }
    for (int t = 0; t < total_site[3]; ++t) {
      gwm_ts[omp_get_thread_num() * total_site[3] + t] = wm_ts[t];
    }
  }
  std::vector<WilsonMatrix> wm_ts(total_site[3]);
  set_zero(wm_ts);
  for (int i = 0; i < omp_get_max_threads(); ++i) {
    for (int t = 0; t < total_site[3]; ++t) {
      wm_ts[t] += gwm_ts[i * total_site[3] + t];
    }
  }
  glb_sum_double_vec(get_data(wm_ts));
  const double ratio = 1.0 / fsel.prob;
  for (int t = 0; t < total_site[3]; ++t) {
    wm_ts[t] *= (Complex)ratio;
  }
  return wm_ts;
}

inline void refresh_wall_snk_prop(WallSrcProps& wsp, const FieldSelection& fsel)
{
  TIMER_VERBOSE("refresh_wall_snk_prop");
  const Coordinate total_site = fsel.f_rank.geo.total_site();
  wsp.sloppy_wall_snk.resize(total_site[3]);
  wsp.exact_wall_snk.resize(total_site[3]);
  for (int i = 0; i < total_site[3]; ++i) {
    {
      const SelProp& prop = wsp.sloppy[i];
      if (prop.initialized) {
        wsp.sloppy_wall_snk[i] = contract_wall_snk_prop(prop, fsel);
      }
    }
    {
      const SelProp& prop = wsp.exact[i];
      if (prop.initialized) {
        wsp.exact_wall_snk[i] = contract_wall_snk_prop(prop, fsel);
      }
    }
  }
}

inline void refresh_prob(WallSrcProps& wsp, const Coordinate& total_site,
                         const int num_exact)
// num_exact: number of random choices (include duplicates)
{
  TIMER_VERBOSE("refresh_prob");
  wsp.exact_tslice_mask.resize(total_site[3], false);
  wsp.num_sloppy = 0;
  wsp.num_exact = 0;
  for (int i = 0; i < total_site[3]; ++i) {
    if (wsp.sloppy[i].initialized) {
      wsp.num_sloppy += 1;
    }
    if (wsp.exact[i].initialized) {
      wsp.exact_tslice_mask[i] = true;
      wsp.num_exact += 1;
    } else {
      wsp.exact_tslice_mask[i] = false;
    }
  }
  qassert(wsp.num_sloppy == total_site[3]);
  qassert(wsp.num_exact <= num_exact);
  wsp.num_exact = num_exact;
  const double prob = 1.0 / (double)wsp.num_sloppy;
  const double prob_1 = 1.0 - std::pow(1.0 - prob, wsp.num_exact);
  const double prob_00 =
      std::pow(1.0 - 2.0 / (double)wsp.num_sloppy, wsp.num_exact);
  const double prob_01 =
      std::pow(1.0 - 1.0 / (double)wsp.num_sloppy, wsp.num_exact) *
      (1.0 - std::pow(1.0 - 1.0 / (double)(wsp.num_sloppy - 1), wsp.num_exact));
  const double prob_11 = 1.0 - prob_00 - 2 * prob_01;
  wsp.sloppy_exact_ratio_1 = 1.0 / prob_1;
  wsp.sloppy_exact_ratio_11 = 1.0 / prob_11;
}

inline void refresh_prop_with_gt(WallSrcProps& wsp, const GaugeTransform& gt,
                                 const PointSelection& psel,
                                 const FieldSelection& fsel)
// need refresh_wall_snk_prop before this func
{
  TIMER_VERBOSE("refresh_prop_with_gt");
  const Coordinate total_site = fsel.f_rank.geo.total_site();
  GaugeTransform gt_inv;
  gt_invert(gt_inv, gt);
  for (int i = 0; i < total_site[3]; ++i) {
    {
      SelProp& prop = wsp.sloppy[i];
      if (prop.initialized) {
        prop_apply_gauge_transformation(prop, prop, gt_inv, fsel);
      }
    }
    {
      SelProp& prop = wsp.exact[i];
      if (prop.initialized) {
        prop_apply_gauge_transformation(prop, prop, gt_inv, fsel);
      }
    }
    {
      PselProp& prop = wsp.sloppy_point_snk[i];
      if (prop.initialized) {
        prop_apply_gauge_transformation(prop, prop, gt_inv, psel);
      }
    }
    {
      PselProp& prop = wsp.exact_point_snk[i];
      if (prop.initialized) {
        prop_apply_gauge_transformation(prop, prop, gt_inv, psel);
      }
    }
  }
}

inline void refresh_wsp(WallSrcProps& wsp, const int num_exact,
                        const GaugeTransform& gt, const PointSelection& psel,
                        const FieldSelection& fsel)
// interface function
{
  TIMER_VERBOSE("refresh_wsp");
  const Coordinate total_site = fsel.f_rank.geo.total_site();
  refresh_wall_snk_prop(wsp, fsel);
  refresh_prob(wsp, total_site, num_exact);
  refresh_prop_with_gt(wsp, gt, psel, fsel);
}

inline const SelProp& get_prop(const WallSrcProps& wsp, const int tslice, const bool exact)
{
  if (exact) {
    return wsp.exact[tslice];
  } else {
    return wsp.sloppy[tslice];
  }
}

inline const PselProp& get_psel_prop(const WallSrcProps& wsp, const int tslice, const bool exact)
{
  if (exact) {
    return wsp.exact_point_snk[tslice];
  } else {
    return wsp.sloppy_point_snk[tslice];
  }
}

}  // namespace qlat
