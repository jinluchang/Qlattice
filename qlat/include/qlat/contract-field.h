#pragma once

#include <qlat/contract-wall-src-prop.h>

namespace qlat
{  //

template <class M>
void acc_field(Field<M>& f, const Complex coef, const SelectedField<M>& sf,
               const ShiftShufflePlan& ssp)
{
  TIMER("acc_field(f,coef,sf,ssp)");
  SelectedField<M> sf_shifted;
  field_shift(sf_shifted, sf, ssp);
  acc_field(f, coef, sf_shifted, ssp.fsel);
}

inline void contract_psel_fsel_distribution_acc(FieldM<Complex, 1>& pos,
                                                const Coordinate& xg_y,
                                                const FieldSelection& fsel,
                                                const ShiftShufflePlan& ssp)
// ssp = make_shift_shuffle_plan(fsel, -xg_y);
{
  TIMER_VERBOSE("contract_psel_fsel_distribution");
  qassert(ssp.shift == -xg_y);
  qassert(ssp.is_reflect == false);
  SelectedField<Complex> s_pos;
  s_pos.init(fsel, 1);
#pragma omp parallel for
  for (long idx = 0; idx < fsel.n_elems; ++idx) {
    s_pos.get_elem(idx) = 1.0;
  }
  qassert(fsel.prob == ssp.fsel.prob);
  const Complex coef = 1.0 / fsel.prob;
  acc_field(pos, coef, s_pos, ssp);
}

template <class M>
void rescale_field_with_psel_fsel_distribution(Field<M>& f,
                                               const FieldM<Complex, 1>& pfdist)
{
  TIMER_VERBOSE("rescale_field_with_psel_fsel_distribution");
  const Geometry& geo = f.geo();
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Complex factor = pfdist.get_elem(xl);
    Vector<M> fv = f.get_elems(xl);
    if (factor != 0.0) {
      const Complex coef = 1.0 / factor;
      for (int m = 0; m < geo.multiplicity; ++m) {
        fv[m] *= coef;
      }
    } else {
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      displayln(fname +
                ssprintf(": WARNING: pfdist[%s] = 0.", show(xg).c_str()));
      qassert(qnorm(fv) == 0.0);
    }
  }
}

API inline array<SpinMatrix, 8>& get_va_matrices()
{
  static bool initialized = false;
  static array<SpinMatrix, 8> ms;
  if (not initialized) {
    initialized = true;
    const array<SpinMatrix, 16>& gms = SpinMatrixConstants::get_cps_gms();
    ms[0] = gms[1];
    ms[1] = gms[2];
    ms[2] = gms[4];
    ms[3] = gms[8];
    ms[4] = gms[14];
    ms[5] = -gms[13];
    ms[6] = gms[11];
    ms[7] = -gms[7];
  }
  return ms;
}

// -----------------------------------------------------------------------------------

inline void field_permute_mu_nu(FieldM<Complex, 8 * 8>& f)
{
  TIMER_VERBOSE("field_permute_mu_nu");
  const Geometry& geo = f.geo();
  FieldM<Complex, 8 * 8> f0;
  f0 = f;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Vector<Complex> fv0 = f0.get_elems_const(index);
    Vector<Complex> fv = f.get_elems(index);
    for (int mu = 0; mu < 8; ++mu) {
      for (int nu = 0; nu < 8; ++nu) {
        fv[mu * 8 + nu] = fv0[nu * 8 + mu];
      }
    }
  }
}

inline void field_conjugate_mu_nu(FieldM<Complex, 8 * 8>& f)
{
  TIMER_VERBOSE("field_conjugate_mu_nu");
  const Geometry& geo = f.geo();
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Vector<Complex> fv = f.get_elems(index);
    for (int mu = 0; mu < 8; ++mu) {
      const double theta_mu = mu < 4 ? 1.0 : -1.0;
      for (int nu = 0; nu < 8; ++nu) {
        const double theta_nu = nu < 4 ? 1.0 : -1.0;
        const double coef = theta_mu * theta_nu;
        fv[mu * 8 + nu] = coef * fv[mu * 8 + nu];
      }
    }
  }
}

inline void field_complex_conjugate(Field<Complex>& f)
{
  TIMER_VERBOSE("field_complex_conjugate");
  const Geometry& geo = f.geo();
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Vector<Complex> fv = f.get_elems(index);
    for (int m = 0; m < geo.multiplicity; ++m) {
      fv[m] = qconj(fv[m]);
    }
  }
}

// -----------------------------------------------------------------------------------

inline void contract_meson_vv_unshifted_acc_x(
    Vector<Complex> v, const Complex coef, const WallSrcProps& wsp1,
    const WallSrcProps& wsp2, const WilsonMatrix& wm3_x_y,
    const Coordinate& xg_x, const long xg_x_idx, const Coordinate& xg_y,
    const long xg_y_psel_idx, const int t_wall, const bool exact)
{
  (void)xg_x;
  (void)xg_y;
  const array<SpinMatrix, 8>& va_ms = get_va_matrices();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  if (exact) {
    qassert(wsp1.exact_tslice_mask[t_wall]);
    qassert(wsp2.exact_tslice_mask[t_wall]);
  }
  const WilsonMatrix& wm1_x_tsrc =
      get_prop(wsp1, t_wall, exact).get_elem(xg_x_idx);
  const WilsonMatrix& wm2_y_tsrc =
      get_psel_prop(wsp2, t_wall, exact).get_elem(xg_y_psel_idx);
  const WilsonMatrix wm1_tsrc_x =
      gamma5 * (WilsonMatrix)matrix_adjoint(wm1_x_tsrc) * gamma5;
  const WilsonMatrix wm_y_tsrc_x = wm2_y_tsrc * gamma5 * wm1_tsrc_x;
  for (int mu = 0; mu < 8; ++mu) {
    const WilsonMatrix wm = wm_y_tsrc_x * va_ms[mu] * wm3_x_y;
    for (int nu = 0; nu < 8; ++nu) {
      const int mu_nu = 8 * mu + nu;
      v[mu_nu] += coef * matrix_trace(wm, va_ms[nu]);
    }
  }
}

inline void contract_meson_vv_unshifted_acc_x(
    Vector<Complex> v, const Complex coef, const WallSrcProps& wsp1,
    const WallSrcProps& wsp2, const WilsonMatrix& wm3_x_y,
    const Coordinate& xg_x, const long xg_x_idx, const Coordinate& xg_y,
    const long xg_y_psel_idx, const int t_wall)
// perform AMA correction for wall src props
{
  qassert(wsp1.exact_tslice_mask.size() == wsp2.exact_tslice_mask.size());
  qassert(0 <= t_wall and t_wall < (int)wsp1.exact_tslice_mask.size());
  const bool has_exact = wsp1.exact_tslice_mask[t_wall];
  qassert(wsp2.exact_tslice_mask[t_wall] == has_exact);
  if (has_exact) {
    const double sloppy_exact_ratio_1 = wsp1.sloppy_exact_ratio_1;
    qassert(sloppy_exact_ratio_1 == wsp2.sloppy_exact_ratio_1);
    const Complex coef1 = sloppy_exact_ratio_1 * coef;
    const Complex coef2 = (1.0 - sloppy_exact_ratio_1) * coef;
    contract_meson_vv_unshifted_acc_x(v, coef1, wsp1, wsp2, wm3_x_y, xg_x,
                                      xg_x_idx, xg_y, xg_y_psel_idx, t_wall,
                                      true);
    contract_meson_vv_unshifted_acc_x(v, coef2, wsp1, wsp2, wm3_x_y, xg_x,
                                      xg_x_idx, xg_y, xg_y_psel_idx, t_wall,
                                      false);
  } else {
    contract_meson_vv_unshifted_acc_x(v, coef, wsp1, wsp2, wm3_x_y, xg_x,
                                      xg_x_idx, xg_y, xg_y_psel_idx, t_wall,
                                      false);
  }
}

inline void contract_meson_vv_unshifted(
    std::vector<SelectedField<Complex> >& sfs, const WallSrcProps& wsp1,
    const WallSrcProps& wsp2, const SelProp& prop3_x_y, const Coordinate& xg_y,
    const long xg_y_psel_idx, const int tsep, const PointSelection& psel,
    const FieldSelection& fsel)
// fsel.prob is NOT accounted.
{
  TIMER_VERBOSE("contract_meson_vv_acc_unshifted");
  qassert(psel[xg_y_psel_idx] == xg_y);
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate total_site = geo.total_site();
  const int multiplicity = 8 * 8;
  clear(sfs);
  sfs.resize(2);
  for (int i = 0; i < (int)sfs.size(); ++i) {
    sfs[i].init(fsel, multiplicity);
    set_zero(sfs[i]);
  }
  SelectedField<Complex>& s_decay = sfs[0];
  SelectedField<Complex>& s_fission = sfs[1];
  const int yt = xg_y[3];
#pragma omp parallel for
  for (long idx = 0; idx < fsel.n_elems; ++idx) {
    const long xg_x_idx = idx;
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Coordinate& xg_x = xg;
    const int xt = xg[3];
    int t_src, t_snk;
    if (smod(xt - yt, total_site[3]) >= 0) {
      t_src = mod(yt - tsep, total_site[3]);
      t_snk = mod(xt + tsep, total_site[3]);
    } else {
      t_src = mod(xt - tsep, total_site[3]);
      t_snk = mod(yt + tsep, total_site[3]);
    }
    const WilsonMatrix& wm3_x_y = prop3_x_y.get_elem(idx);
    Vector<Complex> vd = s_decay.get_elems(idx);
    Vector<Complex> vf = s_fission.get_elems(idx);
    contract_meson_vv_unshifted_acc_x(vd, 1.0, wsp1, wsp2, wm3_x_y, xg_x,
                                      xg_x_idx, xg_y, xg_y_psel_idx, t_src);
    contract_meson_vv_unshifted_acc_x(vf, 1.0, wsp1, wsp2, wm3_x_y, xg_x,
                                      xg_x_idx, xg_y, xg_y_psel_idx, t_snk);
  }
  sync_node();
}

inline void contract_meson_vv_acc(
    FieldM<Complex, 8 * 8>& decay, FieldM<Complex, 8 * 8>& fission,
    const WallSrcProps& wsp1, const WallSrcProps& wsp2,
    const SelProp& prop3_x_y, const Coordinate& xg_y, const long xg_y_psel_idx,
    const int tsep, const PointSelection& psel, const FieldSelection& fsel,
    const ShiftShufflePlan& ssp)
// xg_y = psel[xg_y_psel_idx] is the point src location for prop3_x_y
// ssp = make_shift_shuffle_plan(fsel, -xg_y);
{
  TIMER_VERBOSE("contract_meson_vv_acc");
  const Geometry& geo = fsel.f_rank.geo();
  qassert(is_initialized(decay) == is_initialized(fission));
  qassert(geo == prop3_x_y.geo());
  qassert(fsel.n_elems == prop3_x_y.n_elems);
  qassert(is_initialized(wsp1));
  qassert(is_initialized(wsp2));
  qassert(is_initialized(prop3_x_y));
  qassert(psel[xg_y_psel_idx] == xg_y);
  qassert(ssp.shift == -xg_y);
  qassert(ssp.is_reflect == false);
  std::vector<SelectedField<Complex> > sfs;
  contract_meson_vv_unshifted(sfs, wsp1, wsp2, prop3_x_y, xg_y, xg_y_psel_idx,
                              tsep, psel, fsel);
  qassert(sfs.size() == 2);
  SelectedField<Complex>& s_decay = sfs[0];
  SelectedField<Complex>& s_fission = sfs[1];
  qassert(fsel.prob == ssp.fsel.prob);
  const Complex coef = 1.0 / fsel.prob;
  acc_field(decay, coef, s_decay, ssp);
  acc_field(fission, coef, s_fission, ssp);
}

// -----------------------------------------------------------------------------------

inline WilsonMatrix get_wsnk_prop_avg(const WallSrcProps& wsp, const int t_snk,
                                      const bool exact_snk, const int t_src,
                                      const bool exact_src)
{
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  return (Complex)0.5 *
         (get_wsnk_prop(wsp, t_src, exact_src).get_elem(t_snk) +
          gamma5 *
              (WilsonMatrix)matrix_adjoint(
                  get_wsnk_prop(wsp, t_snk, exact_snk).get_elem(t_src)) *
              gamma5);
}

inline void contract_meson_vv_meson_unshifted_acc_x(
    Vector<Complex> v, const Complex coef, const WallSrcProps& wsp1,
    const WallSrcProps& wsp2, const WallSrcProps& wsp3,
    const WilsonMatrix& wm4_x_y, const Coordinate& xg_x, const long xg_x_idx,
    const Coordinate& xg_y, const long xg_y_psel_idx, const int t_wall_snk,
    const bool exact_snk, const int t_wall_src, const bool exact_src)
{
  (void)xg_x;
  (void)xg_y;
  const array<SpinMatrix, 8>& va_ms = get_va_matrices();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  if (exact_src) {
    qassert(wsp1.exact_tslice_mask[t_wall_src]);
    qassert(wsp2.exact_tslice_mask[t_wall_src]);
    qassert(wsp3.exact_tslice_mask[t_wall_src]);
  }
  if (exact_snk) {
    qassert(wsp1.exact_tslice_mask[t_wall_snk]);
    qassert(wsp2.exact_tslice_mask[t_wall_snk]);
    qassert(wsp3.exact_tslice_mask[t_wall_snk]);
  }
  const WilsonMatrix& wm1_x_tsnk =
      get_prop(wsp1, t_wall_snk, exact_snk).get_elem(xg_x_idx);
  const WilsonMatrix& wm2_y_tsrc =
      get_psel_prop(wsp2, t_wall_src, exact_src).get_elem(xg_y_psel_idx);
  const WilsonMatrix wm1_tsnk_x =
      gamma5 * (WilsonMatrix)matrix_adjoint(wm1_x_tsnk) * gamma5;
  const WilsonMatrix wm3_tsrc_tsnk =
      gamma5 *
      get_wsnk_prop_avg(wsp3, t_wall_src, exact_src, t_wall_snk, exact_snk) *
      gamma5;
  const WilsonMatrix wm_y_tsrc_tsnk_x = wm2_y_tsrc * wm3_tsrc_tsnk * wm1_tsnk_x;
  for (int mu = 0; mu < 8; ++mu) {
    const WilsonMatrix wm = wm_y_tsrc_tsnk_x * va_ms[mu] * wm4_x_y;
    for (int nu = 0; nu < 8; ++nu) {
      const int mu_nu = 8 * mu + nu;
      v[mu_nu] += coef * matrix_trace(wm, va_ms[nu]);
    }
  }
}

inline void contract_meson_vv_meson_unshifted_acc_x(
    Vector<Complex> v, const Complex coef, const WallSrcProps& wsp1,
    const WallSrcProps& wsp2, const WallSrcProps& wsp3,
    const WilsonMatrix& wm4_x_y, const Coordinate& xg_x, const long xg_x_idx,
    const Coordinate& xg_y, const long xg_y_psel_idx, const int t_wall_snk,
    const int t_wall_src)
// perform AMA correction for wall src props
{
  qassert(wsp1.exact_tslice_mask.size() == wsp2.exact_tslice_mask.size());
  qassert(wsp1.exact_tslice_mask.size() == wsp3.exact_tslice_mask.size());
  qassert(0 <= t_wall_snk and t_wall_snk < (int)wsp1.exact_tslice_mask.size());
  qassert(0 <= t_wall_src and t_wall_src < (int)wsp1.exact_tslice_mask.size());
  const bool has_exact_snk = wsp1.exact_tslice_mask[t_wall_snk];
  qassert(wsp2.exact_tslice_mask[t_wall_snk] == has_exact_snk);
  qassert(wsp3.exact_tslice_mask[t_wall_snk] == has_exact_snk);
  const bool has_exact_src = wsp1.exact_tslice_mask[t_wall_src];
  qassert(wsp2.exact_tslice_mask[t_wall_src] == has_exact_src);
  qassert(wsp3.exact_tslice_mask[t_wall_src] == has_exact_src);
  if (t_wall_src == t_wall_snk and has_exact_src) {
    qassert(has_exact_src == has_exact_snk);
    const double sloppy_exact_ratio_1 = wsp1.sloppy_exact_ratio_1;
    qassert(sloppy_exact_ratio_1 == wsp2.sloppy_exact_ratio_1);
    const Complex coef1 = sloppy_exact_ratio_1 * coef;
    const Complex coef2 = (1.0 - sloppy_exact_ratio_1) * coef;
    contract_meson_vv_meson_unshifted_acc_x(v, coef1, wsp1, wsp2, wsp3, wm4_x_y,
                                            xg_x, xg_x_idx, xg_y, xg_y_psel_idx,
                                            t_wall_snk, true, t_wall_src, true);
    contract_meson_vv_meson_unshifted_acc_x(
        v, coef2, wsp1, wsp2, wsp3, wm4_x_y, xg_x, xg_x_idx, xg_y,
        xg_y_psel_idx, t_wall_snk, false, t_wall_src, false);
  } else if (has_exact_src and has_exact_snk) {
    const double sloppy_exact_ratio_11 = wsp1.sloppy_exact_ratio_11;
    const double sloppy_exact_ratio_1 = wsp1.sloppy_exact_ratio_1;
    qassert(sloppy_exact_ratio_1 == wsp2.sloppy_exact_ratio_1);
    const Complex coef1 = sloppy_exact_ratio_11 * coef;
    const Complex coef2 = (sloppy_exact_ratio_1 - sloppy_exact_ratio_11) * coef;
    const Complex coef3 =
        (1.0 - 2.0 * sloppy_exact_ratio_1 + sloppy_exact_ratio_11) * coef;
    contract_meson_vv_meson_unshifted_acc_x(v, coef1, wsp1, wsp2, wsp3, wm4_x_y,
                                            xg_x, xg_x_idx, xg_y, xg_y_psel_idx,
                                            t_wall_snk, true, t_wall_src, true);
    contract_meson_vv_meson_unshifted_acc_x(
        v, coef2, wsp1, wsp2, wsp3, wm4_x_y, xg_x, xg_x_idx, xg_y,
        xg_y_psel_idx, t_wall_snk, false, t_wall_src, true);
    contract_meson_vv_meson_unshifted_acc_x(
        v, coef2, wsp1, wsp2, wsp3, wm4_x_y, xg_x, xg_x_idx, xg_y,
        xg_y_psel_idx, t_wall_snk, true, t_wall_src, false);
    contract_meson_vv_meson_unshifted_acc_x(
        v, coef3, wsp1, wsp2, wsp3, wm4_x_y, xg_x, xg_x_idx, xg_y,
        xg_y_psel_idx, t_wall_snk, false, t_wall_src, false);
  } else if (has_exact_snk or has_exact_src) {
    const double sloppy_exact_ratio_1 = wsp1.sloppy_exact_ratio_1;
    qassert(sloppy_exact_ratio_1 == wsp2.sloppy_exact_ratio_1);
    const Complex coef1 = sloppy_exact_ratio_1 * coef;
    const Complex coef2 = (1.0 - sloppy_exact_ratio_1) * coef;
    if (has_exact_snk and (not has_exact_src)) {
      contract_meson_vv_meson_unshifted_acc_x(
          v, coef1, wsp1, wsp2, wsp3, wm4_x_y, xg_x, xg_x_idx, xg_y,
          xg_y_psel_idx, t_wall_snk, true, t_wall_src, false);
    } else if ((not has_exact_snk) and has_exact_src) {
      contract_meson_vv_meson_unshifted_acc_x(
          v, coef1, wsp1, wsp2, wsp3, wm4_x_y, xg_x, xg_x_idx, xg_y,
          xg_y_psel_idx, t_wall_snk, false, t_wall_src, true);
    } else {
      qassert(false);
    }
    contract_meson_vv_meson_unshifted_acc_x(
        v, coef2, wsp1, wsp2, wsp3, wm4_x_y, xg_x, xg_x_idx, xg_y,
        xg_y_psel_idx, t_wall_snk, false, t_wall_src, false);
  } else if ((not has_exact_snk) and (not has_exact_src)) {
    contract_meson_vv_meson_unshifted_acc_x(
        v, coef, wsp1, wsp2, wsp3, wm4_x_y, xg_x, xg_x_idx, xg_y, xg_y_psel_idx,
        t_wall_snk, false, t_wall_src, false);
  } else {
    qassert(false);
  }
}

inline void contract_meson_vv_meson_unshifted(
    std::vector<SelectedField<Complex> >& meson_vv_meson,
    const WallSrcProps& wsp1, const WallSrcProps& wsp2,
    const WallSrcProps& wsp3, const SelProp& prop4_x_y, const Coordinate& xg_y,
    const long xg_y_psel_idx, const int tsep, const PointSelection& psel,
    const FieldSelection& fsel)
// fsel.prob is NOT accounted.
{
  TIMER_VERBOSE("contract_meson_vv_meson_unshifted");
  qassert(psel[xg_y_psel_idx] == xg_y);
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate total_site = geo.total_site();
  const int multiplicity = 8 * 8;
  clear(meson_vv_meson);
  meson_vv_meson.resize(2);
  for (int i = 0; i < (int)meson_vv_meson.size(); ++i) {
    meson_vv_meson[i].init(fsel, multiplicity);
    set_zero(meson_vv_meson[i]);
  }
  SelectedField<Complex>& meson_vv_meson_forward = meson_vv_meson[0];
  SelectedField<Complex>& meson_vv_meson_backward = meson_vv_meson[1];
  const int yt = xg_y[3];
#pragma omp parallel for
  for (long idx = 0; idx < fsel.n_elems; ++idx) {
    const long xg_x_idx = idx;
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Coordinate& xg_x = xg;
    const int xt = xg[3];
    int t_src, t_snk;
    if (smod(xt - yt, total_site[3]) >= 0) {
      t_src = mod(yt - tsep, total_site[3]);
      t_snk = mod(xt + tsep, total_site[3]);
    } else {
      t_src = mod(xt - tsep, total_site[3]);
      t_snk = mod(yt + tsep, total_site[3]);
    }
    const WilsonMatrix& wm4_x_y = prop4_x_y.get_elem(idx);
    Vector<Complex> vf = meson_vv_meson_forward.get_elems(idx);
    Vector<Complex> vb = meson_vv_meson_backward.get_elems(idx);
    contract_meson_vv_meson_unshifted_acc_x(vf, 1.0, wsp1, wsp2, wsp3, wm4_x_y,
                                            xg_x, xg_x_idx, xg_y, xg_y_psel_idx,
                                            t_snk, t_src);
    contract_meson_vv_meson_unshifted_acc_x(vb, 1.0, wsp1, wsp2, wsp3, wm4_x_y,
                                            xg_x, xg_x_idx, xg_y, xg_y_psel_idx,
                                            t_src, t_snk);
  }
  sync_node();
}

inline void contract_meson_vv_meson_acc(
    FieldM<Complex, 8 * 8>& forward, FieldM<Complex, 8 * 8>& backward,
    const WallSrcProps& wsp1, const WallSrcProps& wsp2,
    const WallSrcProps& wsp3, const SelProp& prop4_x_y, const Coordinate& xg_y,
    const long xg_y_psel_idx, const int tsep, const PointSelection& psel,
    const FieldSelection& fsel, const ShiftShufflePlan& ssp)
// xg_y = psel[xg_y_psel_idx] is the point src location for prop3_x_y
// ssp = make_shift_shuffle_plan(fsel, -xg_y);
{
  TIMER_VERBOSE("contract_meson_vv_meson_acc");
  const Geometry& geo = fsel.f_rank.geo();
  qassert(is_initialized(forward) == is_initialized(backward));
  qassert(geo == prop4_x_y.geo());
  qassert(fsel.n_elems == prop4_x_y.n_elems);
  qassert(is_initialized(wsp1));
  qassert(is_initialized(wsp2));
  qassert(is_initialized(wsp3));
  qassert(is_initialized(prop4_x_y));
  qassert(psel[xg_y_psel_idx] == xg_y);
  qassert(ssp.shift == -xg_y);
  qassert(ssp.is_reflect == false);
  std::vector<SelectedField<Complex> > sfs;
  contract_meson_vv_meson_unshifted(sfs, wsp1, wsp2, wsp3, prop4_x_y, xg_y,
                                    xg_y_psel_idx, tsep, psel, fsel);
  qassert(sfs.size() == 2);
  SelectedField<Complex>& s_forward = sfs[0];
  SelectedField<Complex>& s_backward = sfs[1];
  qassert(fsel.prob == ssp.fsel.prob);
  const Complex coef = 1.0 / fsel.prob;
  acc_field(forward, coef, s_forward, ssp);
  acc_field(backward, coef, s_backward, ssp);
}

// -----------------------------------------------------------------------------------

inline void contract_chvp(SelectedField<Complex>& chvp,
                          const SelProp& prop1_x_y, const SelProp& prop2_x_y,
                          const FieldSelection& fsel)
// fsel.prob is NOT accounted.
{
  TIMER_VERBOSE("contract_chvp");
  const array<SpinMatrix, 8>& va_ms = get_va_matrices();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  qassert(fsel.n_elems == prop1_x_y.n_elems);
  qassert(fsel.n_elems == prop2_x_y.n_elems);
  const int multiplicity = 8 * 8;
  chvp.init();
  chvp.init(fsel, multiplicity);
  qacc_for(idx, fsel.n_elems, {
    const WilsonMatrix& wm1_x_y = prop1_x_y.get_elem(idx);
    const WilsonMatrix& wm2_x_y = prop2_x_y.get_elem(idx);
    const WilsonMatrix wm2_y_x =
        gamma5 * (WilsonMatrix)matrix_adjoint(wm2_x_y) * gamma5;
    Vector<Complex> chvp_v = chvp.get_elems(idx);
    for (int mu = 0; mu < 8; ++mu) {
      const WilsonMatrix wm = wm2_y_x * va_ms[mu] * wm1_x_y;
      for (int nu = 0; nu < 8; ++nu) {
        const int mu_nu = 8 * mu + nu;
        chvp_v[mu_nu] = matrix_trace(wm, va_ms[nu]);
      }
    }
  });
}

inline void contract_chvp_16(FieldM<Complex, 16>& chvp,
                             const Propagator4d& prop1_x_y,
                             const Propagator4d& prop2_x_y)
// chvp.get_elem(x, mu * 4 + nu) ==
// tr(g5_herm(prop2_x_y.get_elem(x)) * gammas[mu]
// * prop1_x_y.get_elem(x) * gammas[nu])
//
// mu: polarization at sink location x
// nu: polarization at source location y
{
  TIMER_VERBOSE("contract_chvp_16");
  const array<SpinMatrix, 4>& gammas = SpinMatrixConstants::get_cps_gammas();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const Geometry& geo = prop1_x_y.geo();
  qassert(geo.multiplicity == 1);
  qassert(is_matching_geo_mult(geo, prop2_x_y.geo()));
  qassert(prop1_x_y.geo().is_only_local);
  qassert(prop2_x_y.geo().is_only_local);
  chvp.init();
  chvp.init(geo);
  set_zero(chvp);
  qassert(chvp.geo().multiplicity == 16);
  qacc_for(index, geo.local_volume(), {
    const WilsonMatrix& wm1_x_y = prop1_x_y.get_elem(index);
    const WilsonMatrix& wm2_x_y = prop2_x_y.get_elem(index);
    const WilsonMatrix wm2_y_x =
        gamma5 * (WilsonMatrix)matrix_adjoint(wm2_x_y) * gamma5;
    Vector<Complex> chvp_v = chvp.get_elems(index);
    for (int mu = 0; mu < 4; ++mu) {
      const WilsonMatrix wm = wm2_y_x * gammas[mu] * wm1_x_y;
      for (int nu = 0; nu < 4; ++nu) {
        const int mu_nu = 4 * mu + nu;
        chvp_v[mu_nu] = matrix_trace(wm, gammas[nu]);
      }
    }
  });
}

// -----------------------------------------------------------------------------------

inline LatData meson_snk_src_shift(const LatData& ld, const int shift)
{
  TIMER_VERBOSE("meson_snk_src_shift");
  qassert(ld.info.size() == 3);
  const int t_size = ld.info[0].size;
  qassert(t_size == ld.info[1].size);
  qassert(is_lat_info_complex(ld.info));
  LatData ld_s;
  ld_s.info = ld.info;
  lat_data_alloc(ld_s);
  set_zero(ld_s);
#pragma omp parallel for collapse(2)
  for (int tsnk = 0; tsnk < t_size; ++tsnk) {
    for (int tsrc = 0; tsrc < t_size; ++tsrc) {
      const array<int, 2> idx = make_array<int>(tsnk, tsrc);
      const int tsnk_s = mod(tsnk + shift, t_size);
      const int tsrc_s = mod(tsrc + shift, t_size);
      const array<int, 2> idx_s = make_array<int>(tsnk_s, tsrc_s);
      lat_data_cget(ld_s, idx_s)[0] = lat_data_cget_const(ld, idx)[0];
    }
  }
  return ld_s;
}

inline void contract_meson_chvp_acc(FieldM<Complex, 8 * 8>& mchvp,
                                    const LatData& ld_meson_snk_src_1_2,
                                    const FieldM<Complex, 8 * 8>& chvp_3_4,
                                    const int tsep)
// ld_meson_snk_src_1_2 should already shifted to origin (t_y -> 0)
// chvp_3_4 already should shifted to origin (xg_y -> 0)
{
  TIMER_VERBOSE("contract_meson_chvp_acc");
  qassert(is_initialized(ld_meson_snk_src_1_2));
  qassert(is_initialized(chvp_3_4));
  const Geometry& geo = chvp_3_4.geo();
  const Coordinate total_site = geo.total_site();
  if (not is_initialized(mchvp)) {
    mchvp.init(geo);
    set_zero(mchvp);
  }
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const int yt = 0;
    const int xt = xg[3];
    int t_src, t_snk;
    if (smod(xt - yt, total_site[3]) >= 0) {
      t_src = mod(yt - tsep, total_site[3]);
      t_snk = mod(xt + tsep, total_site[3]);
    } else {
      t_src = mod(xt - tsep, total_site[3]);
      t_snk = mod(yt + tsep, total_site[3]);
    }
    const array<int, 2> idx = make_array<int>(t_snk, t_src);
    const Complex mss = lat_data_cget_const(ld_meson_snk_src_1_2, idx)[0];
    Vector<Complex> fv = mchvp.get_elems(xl);
    const Vector<Complex> fv0 = chvp_3_4.get_elems_const(xl);
    for (int m = 0; m < geo.multiplicity; ++m) {
      fv[m] += mss * fv0[m];
    }
  }
}

}  // namespace qlat
