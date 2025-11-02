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

inline LatData contract_pion(const Propagator4d& prop, const Int tslice_src)
{
  TIMER_VERBOSE("contract_pion(prop,tsrc)");
  const Geometry& geo = prop.geo();
  const Coordinate total_site = geo.total_site();
  LatData ld = mk_pion_corr_table(total_site);
  Vector<ComplexD> ldv = lat_data_cget(ld);
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const ComplexD val = qnorm(prop.get_elem(xl));
    const Int tsep = mod(xg[3] - tslice_src, total_site[3]);
    ldv[tsep] += val;
  }
  glb_sum(ld);
  return ld;
}

inline LatData contract_pion(const PselProp& prop, const Int tslice_src,
                             const Geometry& geo, const PointsSelection& psel)
{
  TIMER_VERBOSE("contract_pion(ps_prop,tsrc,geo,psel)");
  const Long n_points = prop.n_points;
  qassert(n_points == (Long)psel.size());
  const Coordinate total_site = geo.total_site();
  LatData ld = mk_pion_corr_table(total_site);
  Vector<ComplexD> ldv = lat_data_cget(ld);
  for (Long idx = 0; idx < n_points; ++idx) {
    const Coordinate& xg = psel[idx];
    const ComplexD val = qnorm(prop.get_elem(idx));
    const Int tsep = mod(xg[3] - tslice_src, total_site[3]);
    ldv[tsep] += val;
  }
  ld *= (RealD)product(total_site) / (RealD)n_points;
  return ld;
}

inline LatData contract_pion(const SelProp& prop, const Int tslice_src,
                             const FieldSelection& fsel)
{
  TIMER_VERBOSE("contract_pion(s_prop,tsrc,fsel)");
  const Geometry& geo = prop.geo();
  const Coordinate total_site = geo.total_site();
  LatData ld = mk_pion_corr_table(total_site);
  Vector<ComplexD> ldv = lat_data_cget(ld);
  for (Long idx = 0; idx < fsel.n_elems; ++idx) {
    const Long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const ComplexD val = qnorm(prop.get_elem(idx));
    const Int tsep = mod(xg[3] - tslice_src, total_site[3]);
    ldv[tsep] += val;
  }
  glb_sum(ld);
  ld *= 1.0 / get_fsel_prob(fsel);
  return ld;
}

inline LatData contract_kaon(const SelProp& prop1, const SelProp& prop2,
                             const Int tslice_src, const FieldSelection& fsel)
{
  TIMER_VERBOSE("contract_kaon(s_prop1,s_prop2,fsel)");
  const Geometry& geo = prop1.geo();
  const Coordinate total_site = geo.total_site();
  LatData ld = mk_pion_corr_table(total_site);
  Vector<ComplexD> ldv = lat_data_cget(ld);
  for (Long idx = 0; idx < fsel.n_elems; ++idx) {
    const Long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const ComplexD val =
        matrix_trace(prop1.get_elem(idx), matrix_adjoint(prop2.get_elem(idx)));
    const Int tsep = mod(xg[3] - tslice_src, total_site[3]);
    ldv[tsep] += val;
  }
  glb_sum(ld);
  ld *= 1.0 / get_fsel_prob(fsel);
  return ld;
}

inline LatData contract_pion_wall_snk(const SelProp& prop, const Int tslice_src,
                                      const FieldSelection& fsel)
// is already sparse corrected
{
  TIMER_VERBOSE("contract_pion_wall_snk(s_prop,fsel)");
  const Geometry& geo = prop.geo();
  const Coordinate total_site = geo.total_site();
  const PselProp wm_ts = contract_wall_snk_prop(prop, fsel);
  qassert(wm_ts.n_points == (Long)total_site[3]);
  LatData ld = mk_pion_corr_table(total_site);
  Vector<ComplexD> ldv = lat_data_cget(ld);
  for (Int t = 0; t < total_site[3]; ++t) {
    const Int tsep = mod(t - tslice_src, total_site[3]);
    ldv[tsep] = qnorm(wm_ts.get_elem(t));
  }
  LatData ld_ps = contract_pion(prop, tslice_src, fsel);
  ld_ps *= 1.0 - 1.0 / get_fsel_prob(fsel);
  ld += ld_ps;
  return ld;
}

inline LatData contract_kaon_wall_snk(const SelProp& prop1,
                                      const SelProp& prop2,
                                      const Int tslice_src,
                                      const FieldSelection& fsel)
// is already sparse corrected
{
  TIMER_VERBOSE("contract_kaon_wall_snk(s_prop1,s_prop2,tsrc,fsel)");
  const Geometry& geo = prop1.geo();
  const Coordinate total_site = geo.total_site();
  const PselProp wm1_ts = contract_wall_snk_prop(prop1, fsel);
  const PselProp wm2_ts = contract_wall_snk_prop(prop2, fsel);
  qassert(wm1_ts.n_points == (Long)total_site[3]);
  qassert(wm2_ts.n_points == (Long)total_site[3]);
  LatData ld = mk_pion_corr_table(total_site);
  Vector<ComplexD> ldv = lat_data_cget(ld);
  for (Int t = 0; t < total_site[3]; ++t) {
    const Int tsep = mod(t - tslice_src, total_site[3]);
    ldv[tsep] =
        matrix_trace(wm1_ts.get_elem(t), matrix_adjoint(wm2_ts.get_elem(t)));
  }
  LatData ld_ps = contract_kaon(prop1, prop2, tslice_src, fsel);
  ld_ps *= 1.0 - 1.0 / get_fsel_prob(fsel);
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
                                           const Int tslice,
                                           const FieldSelection& fsel)
// m_ts[tsep][op_src][op_snk] = trace( (\sum_x prop1(x) gms[op_src] gamma5
// prop2(x)^\dagger gamma5) gms[op_snk] ) 0 <= tsep < total_site[3]
{
  TIMER_VERBOSE("contract_two_point_function");
  const array<SpinMatrix, 16>& gms = SpinMatrixConstants::get_cps_gms();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const Geometry& geo = prop1.geo();
  const Coordinate total_site = geo.total_site();
  vector<array<WilsonMatrix, 16> > gwm_ts(omp_get_max_threads() *
                                          total_site[3]);
  set_zero(gwm_ts);
#pragma omp parallel for
  for (Long idx = 0; idx < (Long)fsel.indices.size(); ++idx) {
    const Long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Int tsep = mod(xg[3] - tslice, total_site[3]);
    const WilsonMatrix wm = prop1.get_elem(idx);
    const WilsonMatrix wmd =
        gamma5 * (WilsonMatrix)matrix_adjoint(prop2.get_elem(idx)) * gamma5;
    for (Int op_src = 0; op_src < 16; ++op_src) {
      gwm_ts[omp_get_thread_num() * total_site[3] + tsep][op_src] +=
          wm * gms[op_src] * wmd;
    }
  }
  for (Int i = 1; i < omp_get_max_threads(); ++i) {
    for (Int t = 0; t < total_site[3]; ++t) {
      for (Int op_src = 0; op_src < 16; ++op_src) {
        gwm_ts[t][op_src] += gwm_ts[i * total_site[3] + t][op_src];
      }
    }
  }
  vector<array<ComplexD, 16 * 16> > m_ts(total_site[3]);
  set_zero(m_ts);
#pragma omp parallel for
  for (Int t = 0; t < total_site[3]; ++t) {
    for (Int op_src = 0; op_src < 16; ++op_src) {
      for (Int op_snk = 0; op_snk < 16; ++op_snk) {
        m_ts[t][op_src * 16 + op_snk] =
            matrix_trace(gwm_ts[t][op_src], gms[op_snk]);
      }
    }
  }
  glb_sum(m_ts);
  LatData ld = mk_two_point_table(total_site);
  set_zero(ld);
  for (Int tsep = 0; tsep < total_site[3]; ++tsep) {
    Vector<ComplexD> m_src_snk = lat_data_cget(ld, make_array(tsep));
    for (Int k = 0; k < 16 * 16; ++k) {
      m_src_snk[k] += m_ts[tsep][k];
    }
  }
  ld *= 1.0 / get_fsel_prob(fsel);
  return ld;
}

inline LatData contract_two_point_wall_snk_function(
    const PselProp& prop1, const PselProp& prop2, const Int tslice,
    const Coordinate& total_site)
// m_ts[tsep][op_src][op_snk] = trace( prop1[t] gms[op_src] gamma5
// prop2[t]^\dagger gamma5 gms[op_snk] ) 0 <= tsep < total_site[3]
{
  TIMER_VERBOSE("contract_two_point_wall_snk_function");
  const array<SpinMatrix, 16>& gms = SpinMatrixConstants::get_cps_gms();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  qassert(prop1.n_points == (Long)total_site[3]);
  qassert(prop2.n_points == (Long)total_site[3]);
  vector<array<ComplexD, 16 * 16> > m_ts(total_site[3]);
  set_zero(m_ts);
#pragma omp parallel for
  for (Int t = 0; t < total_site[3]; ++t) {
    const WilsonMatrix& wm = prop1.get_elem(mod(tslice + t, total_site[3]));
    const WilsonMatrix wmd = gamma5 *
                             (WilsonMatrix)matrix_adjoint(prop2.get_elem(
                                 mod(tslice + t, total_site[3]))) *
                             gamma5;
    for (Int op_src = 0; op_src < 16; ++op_src) {
      const WilsonMatrix wm_t = wm * gms[op_src] * wmd;
      for (Int op_snk = 0; op_snk < 16; ++op_snk) {
        m_ts[t][op_src * 16 + op_snk] = matrix_trace(wm_t, gms[op_snk]);
      }
    }
  }
  LatData ld = mk_two_point_table(total_site);
  set_zero(ld);
  for (Int tsep = 0; tsep < total_site[3]; ++tsep) {
    Vector<ComplexD> m_src_snk = lat_data_cget(ld, make_array(tsep));
    for (Int k = 0; k < 16 * 16; ++k) {
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
  ld *= 1.0 - 1.0 / get_fsel_prob(fsel);
  ld += ld_two_point_wall_snk_func;
  return ld;
}

inline LatData contract_two_point_wall_snk_function(const SelProp& prop1,
                                                    const SelProp& prop2,
                                                    const Int tslice,
                                                    const FieldSelection& fsel)
// with sparse correction
// m_ts[tsep][op_src][op_snk] = trace( prop1[t] gms[op_src] gamma5
// prop2[t]^\dagger gamma5 gms[op_snk] ) 0 <= tsep < total_site[3]
{
  TIMER_VERBOSE("contract_two_point_wall_snk_function");
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate total_site = geo.total_site();
  const LatData ld_two_point_func =
      contract_two_point_function(prop1, prop2, tslice, fsel);
  const PselProp wm1_ts = contract_wall_snk_prop(prop1, fsel);
  const PselProp wm2_ts = contract_wall_snk_prop(prop2, fsel);
  const LatData ld_two_point_wall_snk_func =
      contract_two_point_wall_snk_function(wm1_ts, wm2_ts, tslice, total_site);
  return contract_two_point_wall_snk_function(ld_two_point_wall_snk_func,
                                              ld_two_point_func, fsel);
}

inline LatData contract_two_point_function(const WallSrcProps& wsp1,
                                           const WallSrcProps& wsp2,
                                           const FieldSelection& fsel)
{
  Timer::autodisplay();
  TIMER_VERBOSE("contract_two_point_function(wsp)");
  qassert(wsp1.sloppy_exact_ratio_1 == wsp2.sloppy_exact_ratio_1);
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate total_site = geo.total_site();
  LatData ld = mk_two_point_table(total_site);
  for (Int tslice = 0; tslice < total_site[3]; ++tslice) {
    qassert(is_initialized(wsp1.sloppy[tslice]));
    qassert(is_initialized(wsp2.sloppy[tslice]));
    const LatData ld_0 = contract_two_point_function(
        wsp1.sloppy[tslice], wsp2.sloppy[tslice], tslice, fsel);
    ld += ld_0;
    qassert(wsp1.exact_tslice_mask[tslice] == wsp2.exact_tslice_mask[tslice]);
    if (wsp1.exact_tslice_mask[tslice]) {
      const LatData ld_1 = contract_two_point_function(
          wsp1.exact[tslice], wsp2.exact[tslice], tslice, fsel);
      ld += wsp1.sloppy_exact_ratio_1 * (ld_1 - ld_0);
    }
  }
  ld *= 1.0 / (RealD)total_site[3];
  return ld;
}

inline LatData contract_two_point_wall_snk_function(const WallSrcProps& wsp1,
                                                    const WallSrcProps& wsp2,
                                                    const FieldSelection& fsel)
// need to be sparse corrected
{
  TIMER_VERBOSE("contract_two_point_wall_snk_function(wsp)");
  qassert(wsp1.sloppy_exact_ratio_1 == wsp2.sloppy_exact_ratio_1);
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate total_site = geo.total_site();
  LatData ld = mk_two_point_table(total_site);
  for (Int tslice = 0; tslice < total_site[3]; ++tslice) {
    qassert(is_initialized(wsp1.sloppy[tslice]));
    qassert(is_initialized(wsp2.sloppy[tslice]));
    const LatData ld_0 = contract_two_point_wall_snk_function(
        wsp1.sloppy_wall_snk[tslice], wsp2.sloppy_wall_snk[tslice], tslice,
        total_site);
    ld += ld_0;
    qassert(wsp1.exact_tslice_mask[tslice] == wsp2.exact_tslice_mask[tslice]);
    if (wsp1.exact_tslice_mask[tslice]) {
      const LatData ld_1 = contract_two_point_wall_snk_function(
          wsp1.exact_wall_snk[tslice], wsp2.exact_wall_snk[tslice], tslice,
          total_site);
      ld += wsp1.sloppy_exact_ratio_1 * (ld_1 - ld_0);
    }
  }
  ld *= 1.0 / (RealD)total_site[3];
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
                                             const Int ta, const Int tb,
                                             const FieldSelection& fsel)
// ``wm_ab'' is prop from ``tb'' to ``ta''.
// |  ->- prop_a ->- op ->- inv prop_b ->- |
// a (gamma5)                              b (gamma5)
// |            -<- wm_ab -<-              |
//
// prop_a (type1)
// prop_b (type2)
// wm_ab (type3)
{
  TIMER("contract_three_point_function");
  const array<SpinMatrix, 16>& gms = SpinMatrixConstants::get_cps_gms();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate total_site = geo.total_site();
  qassert(is_matching_geo(prop_a.geo(), geo));
  qassert(is_matching_geo(prop_b.geo(), geo));
  vector<WilsonMatrix> gwm_ts(omp_get_max_threads() * total_site[3]);
  set_zero(gwm_ts);
#pragma omp parallel for
  for (Long idx = 0; idx < (Long)fsel.indices.size(); ++idx) {
    const Long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    gwm_ts[omp_get_thread_num() * total_site[3] + xg[3]] +=
        (WilsonMatrix)(prop_a.get_elem(idx) * gamma5 * wm_ab) * gamma5 *
        (gamma5 * (WilsonMatrix)matrix_adjoint(prop_b.get_elem(idx)) * gamma5);
  }
  vector<WilsonMatrix> wm_ts(total_site[3]);
  set_zero(wm_ts);
  for (Int i = 0; i < omp_get_max_threads(); ++i) {
    for (Int t = 0; t < total_site[3]; ++t) {
      wm_ts[t] += gwm_ts[i * total_site[3] + t];
    }
  }
  glb_sum(wm_ts);
  LatData ld = mk_three_point_table(total_site);
  const Int tsep = mod(tb - ta, total_site[3]);
  for (Int t = 0; t < total_site[3]; ++t) {
    const Int top = mod(t - ta, total_site[3]);
    Vector<ComplexD> v = lat_data_cget(ld, make_array(tsep, top));
    for (Int op = 0; op < 16; ++op) {
      v[op] = matrix_trace(wm_ts[t], gms[op]);
    }
  }
  ld *= 1.0 / get_fsel_prob(fsel);
  return ld;
}

inline LatData contract_three_point_function(
    const WallSrcProps& wsp1, const WallSrcProps& wsp2,
    const WallSrcProps& wsp3, const FieldSelection& fsel,
    const Int yt_measurement_sparsity = 1, const Int yt_measurement_start = 0)
{
  TIMER_VERBOSE("compute_three_point_function");
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate total_site = geo.total_site();
  qassert(total_site[3] % yt_measurement_sparsity == 0);
  qassert(0 <= yt_measurement_start and
          yt_measurement_start < yt_measurement_sparsity);
  qassert(wsp1.initialized);
  qassert(wsp2.initialized);
  qassert(wsp3.initialized);
  qassert((int)wsp1.exact_tslice_mask.size() == total_site[3]);
  qassert((int)wsp2.exact_tslice_mask.size() == total_site[3]);
  qassert((int)wsp3.exact_tslice_mask.size() == total_site[3]);
  for (Int i = 0; i < total_site[3]; ++i) {
    qassert(is_initialized(wsp1.sloppy[i]));
    qassert(is_initialized(wsp2.sloppy[i]));
    qassert(is_initialized(wsp3.sloppy[i]));
    qassert(wsp1.exact_tslice_mask[i] == wsp2.exact_tslice_mask[i]);
    qassert(wsp1.exact_tslice_mask[i] == wsp3.exact_tslice_mask[i]);
  }
  LatData ld;
  for (Int tslice = yt_measurement_start; tslice < total_site[3];
       tslice += yt_measurement_sparsity) {
    Timer::autodisplay();
    TIMER_VERBOSE("compute_three_point_function-tslice");
    for (Int tsep = 0; tsep < total_site[3]; ++tsep) {
      const Int ta = tslice;
      const Int tb = mod(ta + tsep, total_site[3]);
      LatData ld_00, ld_10, ld_01, ld_11;
      ld_00 = contract_three_point_function(
          wsp1.sloppy[ta], wsp2.sloppy[tb],
          wsp3.sloppy_wall_snk[tb].get_elem(ta), ta, tb, fsel);
      ld += ld_00;
      if (ta != tb) {
        if (wsp1.exact_tslice_mask[ta]) {
          ld_10 = contract_three_point_function(
              wsp1.exact[ta], wsp2.sloppy[tb],
              wsp3.sloppy_wall_snk[tb].get_elem(ta), ta, tb, fsel);
          ld += wsp1.sloppy_exact_ratio_1 * (ld_10 - ld_00);
        }
        if (wsp1.exact_tslice_mask[tb]) {
          ld_01 = contract_three_point_function(
              wsp1.sloppy[ta], wsp2.exact[tb],
              wsp3.exact_wall_snk[tb].get_elem(ta), ta, tb, fsel);
          ld += wsp1.sloppy_exact_ratio_1 * (ld_01 - ld_00);
        }
      }
      if (wsp1.exact_tslice_mask[ta] and wsp1.exact_tslice_mask[tb]) {
        ld_11 = contract_three_point_function(
            wsp1.exact[ta], wsp2.exact[tb],
            wsp3.exact_wall_snk[tb].get_elem(ta), ta, tb, fsel);
        if (ta == tb) {
          ld += wsp1.sloppy_exact_ratio_1 * (ld_11 - ld_00);
        } else {
          ld += wsp1.sloppy_exact_ratio_11 * (ld_11 - ld_10 - ld_01 + ld_00);
        }
      }
    }
  }
  ld *= (RealD)yt_measurement_sparsity / (RealD)total_site[3];
  return ld;
}

// -----------------------------------------------------------------------------------

inline LatData mk_meson_snk_src_table(const Coordinate& total_site)
{
  LatData ld;
  ld.info.push_back(lat_dim_number("tsnk", 0, total_site[3] - 1));
  ld.info.push_back(lat_dim_number("tsrc", 0, total_site[3] - 1));
  ld.info.push_back(lat_dim_re_im());
  lat_data_alloc(ld);
  set_zero(ld);
  return ld;
}

inline ComplexD contract_meson_snk_src(const WallSrcProps& wsp1,
                                      const WallSrcProps& wsp2, const Int t_snk,
                                      const bool exact_snk, const Int t_src,
                                      const bool exact_src)
{
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const WilsonMatrix& wm1_snk_src =
      get_wsnk_prop(wsp1, t_src, exact_src).get_elem(t_snk);
  const WilsonMatrix& wm2_snk_src =
      get_wsnk_prop(wsp2, t_src, exact_src).get_elem(t_snk);
  const WilsonMatrix& wm1_src_snk =
      get_wsnk_prop(wsp1, t_snk, exact_snk).get_elem(t_src);
  const WilsonMatrix& wm2_src_snk =
      get_wsnk_prop(wsp2, t_snk, exact_snk).get_elem(t_src);
  const ComplexD v1 = matrix_trace(gamma5 * wm1_snk_src, gamma5 * wm2_src_snk);
  const ComplexD v2 = matrix_trace(gamma5 * wm2_snk_src, gamma5 * wm1_src_snk);
  return 0.5 * (v1 + qconj(v2));
}

inline ComplexD contract_meson_snk_src(const WallSrcProps& wsp1,
                                      const WallSrcProps& wsp2, const Int t_snk,
                                      const Int t_src)
{
  qassert(wsp1.exact_tslice_mask.size() == wsp2.exact_tslice_mask.size());
  qassert(0 <= t_snk and t_snk < (int)wsp1.exact_tslice_mask.size());
  qassert(0 <= t_src and t_src < (int)wsp1.exact_tslice_mask.size());
  const bool has_exact_snk = wsp1.exact_tslice_mask[t_snk];
  qassert(wsp2.exact_tslice_mask[t_snk] == has_exact_snk);
  const bool has_exact_src = wsp1.exact_tslice_mask[t_src];
  qassert(wsp2.exact_tslice_mask[t_src] == has_exact_src);
  ComplexD ret = 0.0;
  if (t_src == t_snk and has_exact_src) {
    qassert(has_exact_src == has_exact_snk);
    const RealD sloppy_exact_ratio_1 = wsp1.sloppy_exact_ratio_1;
    qassert(sloppy_exact_ratio_1 == wsp2.sloppy_exact_ratio_1);
    const ComplexD coef1 = sloppy_exact_ratio_1;
    const ComplexD coef2 = 1.0 - sloppy_exact_ratio_1;
    ret += coef1 * contract_meson_snk_src(wsp1, wsp2, t_snk, true, t_src, true);
    ret +=
        coef2 * contract_meson_snk_src(wsp1, wsp2, t_snk, false, t_src, false);
  } else if (has_exact_src and has_exact_snk) {
    const RealD sloppy_exact_ratio_11 = wsp1.sloppy_exact_ratio_11;
    const RealD sloppy_exact_ratio_1 = wsp1.sloppy_exact_ratio_1;
    qassert(sloppy_exact_ratio_1 == wsp2.sloppy_exact_ratio_1);
    const ComplexD coef1 = sloppy_exact_ratio_11;
    const ComplexD coef2 = sloppy_exact_ratio_1 - sloppy_exact_ratio_11;
    const ComplexD coef3 =
        1.0 - 2.0 * sloppy_exact_ratio_1 + sloppy_exact_ratio_11;
    ret += coef1 * contract_meson_snk_src(wsp1, wsp2, t_snk, true, t_src, true);
    ret +=
        coef2 * contract_meson_snk_src(wsp1, wsp2, t_snk, false, t_src, true);
    ret +=
        coef2 * contract_meson_snk_src(wsp1, wsp2, t_snk, true, t_src, false);
    ret +=
        coef3 * contract_meson_snk_src(wsp1, wsp2, t_snk, false, t_src, false);
  } else if (has_exact_snk or has_exact_src) {
    const RealD sloppy_exact_ratio_1 = wsp1.sloppy_exact_ratio_1;
    qassert(sloppy_exact_ratio_1 == wsp2.sloppy_exact_ratio_1);
    const ComplexD coef1 = sloppy_exact_ratio_1;
    const ComplexD coef2 = 1.0 - sloppy_exact_ratio_1;
    if (has_exact_snk and (not has_exact_src)) {
      ret +=
          coef1 * contract_meson_snk_src(wsp1, wsp2, t_snk, true, t_src, false);
    } else if ((not has_exact_snk) and has_exact_src) {
      ret +=
          coef1 * contract_meson_snk_src(wsp1, wsp2, t_snk, false, t_src, true);
    } else {
      qassert(false);
    }
    ret +=
        coef2 * contract_meson_snk_src(wsp1, wsp2, t_snk, false, t_src, false);
  } else if ((not has_exact_snk) and (not has_exact_src)) {
    ret += contract_meson_snk_src(wsp1, wsp2, t_snk, false, t_src, false);
  } else {
    qassert(false);
  }
  return ret;
}

inline LatData contract_meson_snk_src(const WallSrcProps& wsp1,
                                      const WallSrcProps& wsp2,
                                      const Coordinate& total_site)
{
  TIMER_VERBOSE("contract_meson");
  LatData ld = mk_meson_snk_src_table(total_site);
  for (Int t_snk = 0; t_snk < total_site[3]; ++t_snk) {
    for (Int t_src = 0; t_src < total_site[3]; ++t_src) {
      lat_data_cget(ld, make_array<int>(t_snk, t_src))[0] =
          contract_meson_snk_src(wsp1, wsp2, t_snk, t_src);
    }
  }
  return ld;
}

}  // namespace qlat
