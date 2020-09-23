#pragma once

#include <qlat/contract-wall-src-prop.h>

namespace qlat
{  //

inline LatData mk_pion_corr_table(const Coordinate& total_site)
{
  LatData ld;
  ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
  ld.info.push_back(lat_dim_re_im());
  lat_data_alloc(ld);
  set_zero(ld);
  return ld;
}

inline LatData contract_pion(const Propagator4d& prop, const int tslice_src)
{
  TIMER_VERBOSE("contract_pion(prop,tsrc)");
  const Geometry& geo = prop.geo;
  const Coordinate total_site = geo.total_site();
  LatData ld = mk_pion_corr_table(total_site);
  Vector<Complex> ldv = lat_data_cget(ld);
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Complex val = qnorm(prop.get_elem(xl));
    const int tsep = mod(xg[3] - tslice_src, total_site[3]);
    ldv[tsep] += val;
  }
  glb_sum_lat_data(ld);
  return ld;
}

inline LatData contract_pion(const PselProp& prop, const int tslice_src,
                             const Geometry& geo, const PointSelection& psel)
{
  TIMER_VERBOSE("contract_pion(ps_prop,tsrc,geo,psel)");
  const long n_points = prop.n_points;
  qassert(n_points == (long)psel.size());
  const Coordinate total_site = geo.total_site();
  LatData ld = mk_pion_corr_table(total_site);
  Vector<Complex> ldv = lat_data_cget(ld);
  for (long idx = 0; idx < n_points; ++idx) {
    const Coordinate& xg = psel[idx];
    const Complex val = qnorm(prop.get_elem(idx));
    const int tsep = mod(xg[3] - tslice_src, total_site[3]);
    ldv[tsep] += val;
  }
  ld *= (double)product(total_site) / (double)n_points;
  return ld;
}

inline LatData contract_pion(const SelProp& prop, const int tslice_src,
                             const FieldSelection& fsel)
{
  TIMER_VERBOSE("contract_pion(s_prop,tsrc,fsel)");
  const Geometry& geo = prop.geo;
  const Coordinate total_site = geo.total_site();
  LatData ld = mk_pion_corr_table(total_site);
  Vector<Complex> ldv = lat_data_cget(ld);
  for (long idx = 0; idx < fsel.n_elems; ++idx) {
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Complex val = qnorm(prop.get_elem(idx));
    const int tsep = mod(xg[3] - tslice_src, total_site[3]);
    ldv[tsep] += val;
  }
  glb_sum_lat_data(ld);
  ld *= 1.0 / fsel.prob;
  return ld;
}

inline LatData contract_kaon(const SelProp& prop1, const SelProp& prop2,
                             const int tslice_src, const FieldSelection& fsel)
{
  TIMER_VERBOSE("contract_kaon(s_prop1,s_prop2,fsel)");
  const Geometry& geo = prop1.geo;
  const Coordinate total_site = geo.total_site();
  LatData ld = mk_pion_corr_table(total_site);
  Vector<Complex> ldv = lat_data_cget(ld);
  for (long idx = 0; idx < fsel.n_elems; ++idx) {
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Complex val =
        matrix_trace(prop1.get_elem(idx), matrix_adjoint(prop2.get_elem(idx)));
    const int tsep = mod(xg[3] - tslice_src, total_site[3]);
    ldv[tsep] += val;
  }
  glb_sum_lat_data(ld);
  ld *= 1.0 / fsel.prob;
  return ld;
}

inline LatData contract_pion_wall_snk(const SelProp& prop, const int tslice_src,
                                      const FieldSelection& fsel)
// is already sparse corrected
{
  TIMER_VERBOSE("contract_pion_wall_snk(s_prop,fsel)");
  const Geometry& geo = prop.geo;
  const Coordinate total_site = geo.total_site();
  const std::vector<WilsonMatrix> wm_ts = contract_wall_snk_prop(prop, fsel);
  qassert((int)wm_ts.size() == total_site[3]);
  LatData ld = mk_pion_corr_table(total_site);
  Vector<Complex> ldv = lat_data_cget(ld);
  for (int t = 0; t < total_site[3]; ++t) {
    const int tsep = mod(t - tslice_src, total_site[3]);
    ldv[tsep] = qnorm(wm_ts[t]);
  }
  LatData ld_ps = contract_pion(prop, tslice_src, fsel);
  ld_ps *= 1.0 - 1.0 / fsel.prob;
  ld += ld_ps;
  return ld;
}

inline LatData contract_kaon_wall_snk(const SelProp& prop1,
                                      const SelProp& prop2,
                                      const int tslice_src,
                                      const FieldSelection& fsel)
// is already sparse corrected
{
  TIMER_VERBOSE("contract_kaon_wall_snk(s_prop1,s_prop2,tsrc,fsel)");
  const Geometry& geo = prop1.geo;
  const Coordinate total_site = geo.total_site();
  const std::vector<WilsonMatrix> wm1_ts = contract_wall_snk_prop(prop1, fsel);
  const std::vector<WilsonMatrix> wm2_ts = contract_wall_snk_prop(prop2, fsel);
  qassert((int)wm1_ts.size() == total_site[3]);
  qassert((int)wm2_ts.size() == total_site[3]);
  LatData ld = mk_pion_corr_table(total_site);
  Vector<Complex> ldv = lat_data_cget(ld);
  for (int t = 0; t < total_site[3]; ++t) {
    const int tsep = mod(t - tslice_src, total_site[3]);
    ldv[tsep] = matrix_trace(wm1_ts[t], matrix_adjoint(wm2_ts[t]));
  }
  LatData ld_ps = contract_kaon(prop1, prop2, tslice_src, fsel);
  ld_ps *= 1.0 - 1.0 / fsel.prob;
  ld += ld_ps;
  return ld;
}

// -----------------------------------------------------------------------------------

inline LatData mk_two_point_table(const Coordinate& total_site)
{
  LatData ld;
  ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
  ld.info.push_back(lat_dim_number("op-src", 0, 15));
  ld.info.push_back(lat_dim_number("op-snk", 0, 15));
  ld.info.push_back(lat_dim_re_im());
  lat_data_alloc(ld);
  set_zero(ld);
  return ld;
}

inline LatData contract_two_point_function(const SelProp& prop1,
                                           const SelProp& prop2,
                                           const int tslice,
                                           const FieldSelection& fsel)
// m_ts[tsep][op_src][op_snk] = trace( (\sum_x prop1(x) gms[op_src] gamma5
// prop2(x)^\dagger gamma5) gms[op_snk] ) 0 <= tsep < total_site[3]
{
  TIMER_VERBOSE("contract_two_point_function");
  const std::array<SpinMatrix, 16>& gms = SpinMatrixConstants::get_cps_gms();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const Geometry& geo = prop1.geo;
  const Coordinate total_site = geo.total_site();
  std::vector<std::array<WilsonMatrix, 16> > gwm_ts(omp_get_max_threads() *
                                                    total_site[3]);
  set_zero(gwm_ts);
#pragma omp parallel
  {
    std::vector<std::array<WilsonMatrix, 16> > wm_ts(total_site[3]);
    set_zero(wm_ts);
#pragma omp for
    for (long idx = 0; idx < (long)fsel.indices.size(); ++idx) {
      const long index = fsel.indices[idx];
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      const int tsep = mod(xg[3] - tslice, total_site[3]);
      const WilsonMatrix wm = prop1.get_elem(idx);
      const WilsonMatrix wmd =
          gamma5 * (WilsonMatrix)matrix_adjoint(prop2.get_elem(idx)) * gamma5;
      for (int op_src = 0; op_src < 16; ++op_src) {
        wm_ts[tsep][op_src] += wm * gms[op_src] * wmd;
      }
    }
    for (int t = 0; t < total_site[3]; ++t) {
      gwm_ts[omp_get_thread_num() * total_site[3] + t] = wm_ts[t];
    }
  }
  for (int i = 1; i < omp_get_max_threads(); ++i) {
    for (int t = 0; t < total_site[3]; ++t) {
      for (int op_src = 0; op_src < 16; ++op_src) {
        gwm_ts[t][op_src] += gwm_ts[i * total_site[3] + t][op_src];
      }
    }
  }
  std::vector<std::array<Complex, 16 * 16> > m_ts(total_site[3]);
  set_zero(m_ts);
#pragma omp parallel for
  for (int t = 0; t < total_site[3]; ++t) {
    for (int op_src = 0; op_src < 16; ++op_src) {
      for (int op_snk = 0; op_snk < 16; ++op_snk) {
        m_ts[t][op_src * 16 + op_snk] =
            matrix_trace(gwm_ts[t][op_src], gms[op_snk]);
      }
    }
  }
  glb_sum_double_vec(get_data(m_ts));
  LatData ld = mk_two_point_table(total_site);
  set_zero(ld);
  for (int tsep = 0; tsep < total_site[3]; ++tsep) {
    Vector<Complex> m_src_snk = lat_data_complex_get(ld, make_array(tsep));
    for (int k = 0; k < 16 * 16; ++k) {
      m_src_snk[k] += m_ts[tsep][k];
    }
  }
  ld *= 1.0 / fsel.prob;
  return ld;
}

inline LatData contract_two_point_wall_snk_function(
    const std::vector<WilsonMatrix>& prop1,
    const std::vector<WilsonMatrix>& prop2, const int tslice,
    const Coordinate& total_site)
// m_ts[tsep][op_src][op_snk] = trace( prop1[t] gms[op_src] gamma5
// prop2[t]^\dagger gamma5 gms[op_snk] ) 0 <= tsep < total_site[3]
{
  TIMER_VERBOSE("contract_two_point_wall_snk_function");
  const std::array<SpinMatrix, 16>& gms = SpinMatrixConstants::get_cps_gms();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  qassert((int)prop1.size() == total_site[3]);
  qassert((int)prop2.size() == total_site[3]);
  std::vector<std::array<Complex, 16 * 16> > m_ts(total_site[3]);
  set_zero(m_ts);
#pragma omp parallel for
  for (int t = 0; t < total_site[3]; ++t) {
    const WilsonMatrix& wm = prop1[mod(tslice + t, total_site[3])];
    const WilsonMatrix wmd =
        gamma5 *
        (WilsonMatrix)matrix_adjoint(prop2[mod(tslice + t, total_site[3])]) *
        gamma5;
    for (int op_src = 0; op_src < 16; ++op_src) {
      const WilsonMatrix wm_t = wm * gms[op_src] * wmd;
      for (int op_snk = 0; op_snk < 16; ++op_snk) {
        m_ts[t][op_src * 16 + op_snk] = matrix_trace(wm_t, gms[op_snk]);
      }
    }
  }
  LatData ld = mk_two_point_table(total_site);
  set_zero(ld);
  for (int tsep = 0; tsep < total_site[3]; ++tsep) {
    Vector<Complex> m_src_snk = lat_data_complex_get(ld, make_array(tsep));
    for (int k = 0; k < 16 * 16; ++k) {
      m_src_snk[k] += m_ts[tsep][k];
    }
  }
  return ld;
}

inline LatData contract_two_point_wall_snk_function(
    const LatData& ld_two_point_wall_snk_func, const LatData& ld_two_point_func,
    const FieldSelection& fsel)
// perform sparse correction
{
  TIMER_VERBOSE("contract_two_point_wall_snk_function");
  LatData ld = ld_two_point_func;
  ld *= 1.0 - 1.0 / fsel.prob;
  ld += ld_two_point_wall_snk_func;
  return ld;
}

inline LatData contract_two_point_wall_snk_function(const SelProp& prop1,
                                                    const SelProp& prop2,
                                                    const int tslice,
                                                    const FieldSelection& fsel)
// with sparse correction
// m_ts[tsep][op_src][op_snk] = trace( prop1[t] gms[op_src] gamma5
// prop2[t]^\dagger gamma5 gms[op_snk] ) 0 <= tsep < total_site[3]
{
  TIMER_VERBOSE("contract_two_point_wall_snk_function");
  const Geometry& geo = fsel.f_rank.geo;
  const Coordinate total_site = geo.total_site();
  const LatData ld_two_point_func =
      contract_two_point_function(prop1, prop2, tslice, fsel);
  const std::vector<WilsonMatrix> wm1_ts = contract_wall_snk_prop(prop1, fsel);
  const std::vector<WilsonMatrix> wm2_ts = contract_wall_snk_prop(prop2, fsel);
  const LatData ld_two_point_wall_snk_func =
      contract_two_point_wall_snk_function(wm1_ts, wm2_ts, tslice, total_site);
  return contract_two_point_wall_snk_function(ld_two_point_wall_snk_func,
                                              ld_two_point_func, fsel);
}

inline LatData contract_two_point_function(const WallSrcProps& wsp1,
                                           const WallSrcProps& wsp2,
                                           const FieldSelection& fsel)
{
  TIMER_VERBOSE("contract_two_point_function(wsp)");
  qassert(wsp1.sloppy_exact_ratio_1 == wsp2.sloppy_exact_ratio_1);
  const Geometry& geo = fsel.f_rank.geo;
  const Coordinate total_site = geo.total_site();
  LatData ld = mk_two_point_table(total_site);
  for (int tslice = 0; tslice < total_site[3]; ++tslice) {
    Timer::autodisplay();
    qassert(is_initialized(wsp1.sloppy[tslice]));
    qassert(is_initialized(wsp2.sloppy[tslice]));
    const LatData ld_0 = contract_two_point_function(
        wsp1.sloppy[tslice], wsp2.sloppy[tslice], tslice, fsel);
    ld += ld_0;
    qassert(wsp1.exact_tslice_mask[tslice] == wsp2.exact_tslice_mask[tslice]);
    if (wsp1.exact_tslice_mask[tslice]) {
      const LatData ld_1 = contract_two_point_function(
          wsp1.exact[tslice], wsp1.exact[tslice], tslice, fsel);
      ld += wsp1.sloppy_exact_ratio_1 * (ld_1 - ld_0);
    }
  }
  ld *= 1.0 / (double)total_site[3];
  return ld;
}

inline LatData contract_two_point_wall_snk_function(const WallSrcProps& wsp1,
                                                    const WallSrcProps& wsp2,
                                                    const FieldSelection& fsel)
// need to be sparse corrected
{
  TIMER_VERBOSE("contract_two_point_wall_snk_function(wsp)");
  qassert(wsp1.sloppy_exact_ratio_1 == wsp2.sloppy_exact_ratio_1);
  const Geometry& geo = fsel.f_rank.geo;
  const Coordinate total_site = geo.total_site();
  LatData ld = mk_two_point_table(total_site);
  for (int tslice = 0; tslice < total_site[3]; ++tslice) {
    qassert(is_initialized(wsp1.sloppy[tslice]));
    qassert(is_initialized(wsp2.sloppy[tslice]));
    const LatData ld_0 = contract_two_point_wall_snk_function(
        wsp1.sloppy_wall_snk[tslice], wsp2.sloppy_wall_snk[tslice], tslice,
        total_site);
    ld += ld_0;
    qassert(wsp1.exact_tslice_mask[tslice] == wsp2.exact_tslice_mask[tslice]);
    if (wsp1.exact_tslice_mask[tslice]) {
      const LatData ld_1 = contract_two_point_wall_snk_function(
          wsp1.exact_wall_snk[tslice], wsp1.exact_wall_snk[tslice], tslice,
          total_site);
      ld += wsp1.sloppy_exact_ratio_1 * (ld_1 - ld_0);
    }
  }
  ld *= 1.0 / (double)total_site[3];
  return ld;
}

// -----------------------------------------------------------------------------------

inline LatData mk_three_point_table(const Coordinate& total_site)
{
  LatData ld;
  ld.info.push_back(lat_dim_number("tsep", 0, total_site[3] - 1));
  ld.info.push_back(lat_dim_number("top", 0, total_site[3] - 1));
  ld.info.push_back(lat_dim_number("op", 0, 15));
  ld.info.push_back(lat_dim_re_im());
  lat_data_alloc(ld);
  set_zero(ld);
  return ld;
}

inline LatData contract_three_point_function(const SelProp& prop_a,
                                             const SelProp& prop_b,
                                             const WilsonMatrix& wm_ab,
                                             const int ta, const int tb,
                                             const FieldSelection& fsel)
// ``wm_ab'' is prop from ``tb'' to ``ta''.
// |  ->- prop_a ->- op ->- inv prop_b ->- |
// a (gamma5)                              b (gamma5)
// |            -<- wm_ab -<-              |
{
  TIMER_VERBOSE("contract_three_point_function");
  const std::array<SpinMatrix, 16>& gms = SpinMatrixConstants::get_cps_gms();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const Geometry& geo = fsel.f_rank.geo;
  const Coordinate total_site = geo.total_site();
  qassert(is_matching_geo_mult(prop_a.geo, geo));
  qassert(is_matching_geo_mult(prop_b.geo, geo));
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
      wm_ts[xg[3]] +=
          (WilsonMatrix)(prop_a.get_elem(idx) * gamma5 * wm_ab) * gamma5 *
          (gamma5 * (WilsonMatrix)matrix_adjoint(prop_b.get_elem(idx)) *
           gamma5);
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
  LatData ld = mk_three_point_table(total_site);
  const int tsep = mod(tb - ta, total_site[3]);
  for (int t = 0; t < total_site[3]; ++t) {
    const int top = mod(t - ta, total_site[3]);
    Vector<Complex> v = lat_data_complex_get(ld, make_array(tsep, top));
    for (int op = 0; op < 16; ++op) {
      v[op] = matrix_trace(wm_ts[t], gms[op]);
    }
  }
  ld *= 1.0 / fsel.prob;
  return ld;
}

inline LatData contract_three_point_function(const WallSrcProps& wsp1,
                                             const WallSrcProps& wsp2,
                                             const WallSrcProps& wsp3,
                                             const FieldSelection& fsel,
                                             const std::string& job_tag,
                                             const int traj)
// job_tag and traj only used to initialize RngState
{
  TIMER_VERBOSE("compute_three_point_function");
  const Geometry& geo = fsel.f_rank.geo;
  const Coordinate total_site = geo.total_site();
  qassert(wsp1.initialized);
  qassert(wsp2.initialized);
  qassert(wsp3.initialized);
  qassert((int)wsp1.exact_tslice_mask.size() == total_site[3]);
  qassert((int)wsp2.exact_tslice_mask.size() == total_site[3]);
  qassert((int)wsp3.exact_tslice_mask.size() == total_site[3]);
  for (int i = 0; i < total_site[3]; ++i) {
    qassert(is_initialized(wsp1.sloppy[i]));
    qassert(is_initialized(wsp2.sloppy[i]));
    qassert(is_initialized(wsp3.sloppy[i]));
    qassert(wsp1.exact_tslice_mask[i] == wsp2.exact_tslice_mask[i]);
    qassert(wsp1.exact_tslice_mask[i] == wsp3.exact_tslice_mask[i]);
  }
  const int yt_measurement_sparsity = 4;
  qassert(total_site[3] % yt_measurement_sparsity == 0);
  RngState rs = RngState("contract_three_point_function")
                    .split("yt-start-time-slice")
                    .split(job_tag)
                    .split(traj);
  const int yt_measurement_start = rand_gen(rs) % yt_measurement_sparsity;
  qassert(0 <= yt_measurement_start and
          yt_measurement_start < yt_measurement_sparsity);
  LatData ld;
  for (int tslice = yt_measurement_start; tslice < total_site[3];
       tslice += yt_measurement_sparsity) {
    Timer::autodisplay();
    TIMER_VERBOSE("compute_three_point_function-tslice");
    for (int tsep = 0; tsep < total_site[3]; ++tsep) {
      const int ta = tslice;
      const int tb = mod(ta + tsep, total_site[3]);
      LatData ld_00, ld_10, ld_01, ld_11;
      ld_00 = contract_three_point_function(wsp1.sloppy[ta], wsp2.sloppy[tb],
                                            wsp3.sloppy_wall_snk[tb][ta], ta,
                                            tb, fsel);
      ld += ld_00;
      if (ta != tb) {
        if (wsp1.exact_tslice_mask[ta]) {
          ld_10 = contract_three_point_function(wsp1.exact[ta], wsp2.sloppy[tb],
                                                wsp3.sloppy_wall_snk[tb][ta],
                                                ta, tb, fsel);
          ld += wsp1.sloppy_exact_ratio_1 * (ld_10 - ld_00);
        }
        if (wsp1.exact_tslice_mask[tb]) {
          ld_01 = contract_three_point_function(wsp1.sloppy[ta], wsp2.exact[tb],
                                                wsp3.exact_wall_snk[tb][ta], ta,
                                                tb, fsel);
          ld += wsp1.sloppy_exact_ratio_1 * (ld_01 - ld_00);
        }
      }
      if (wsp1.exact_tslice_mask[ta] and wsp1.exact_tslice_mask[tb]) {
        ld_11 = contract_three_point_function(wsp1.exact[ta], wsp2.exact[tb],
                                              wsp3.exact_wall_snk[tb][ta], ta,
                                              tb, fsel);
        if (ta == tb) {
          ld += wsp1.sloppy_exact_ratio_1 * (ld_11 - ld_00);
        } else {
          ld += wsp1.sloppy_exact_ratio_11 * (ld_11 - ld_10 - ld_01 + ld_00);
        }
      }
    }
  }
  ld *= (double)yt_measurement_sparsity / (double)total_site[3];
  return ld;
}

}  // namespace qlat
