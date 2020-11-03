#pragma once

#include <qlat/contract-wall-src-prop.h>

namespace qlat
{  //

inline void contract_psel_fsel_distribution(FieldM<Complex, 1>& pos,
                                            const Coordinate& xg_y,
                                            const long xg_y_psel_idx,
                                            const PointSelection& psel,
                                            const FieldSelection& fsel,
                                            const ShiftShufflePlan& ssp)
// xg_y = psel[xg_y_psel_idx] is the point src location for prop3_x_y
// ssp = make_shift_shuffle_plan(fsel, -xg_y);
{
  TIMER_VERBOSE("contract_psel_fsel_distribution");
  qassert(psel[xg_y_psel_idx] == xg_y);
  qassert(ssp.shift == -xg_y);
  qassert(ssp.is_reflect == false);
  SelectedField<Complex> s_pos;
  s_pos.init(fsel, 1);
#pragma omp parallel for
  for (long idx = 0; idx < fsel.n_elems; ++idx) {
    s_pos.get_elem(idx) = 1.0;
  }
  field_shift(s_pos, s_pos, ssp);
  qassert(fsel.prob == ssp.fsel.prob);
  const Complex coef = 1.0 / fsel.prob;
  acc_field(pos, coef, s_pos, ssp.fsel);
}

template <class M>
void rescale_field_with_psel_fsel_distribution(Field<M>& f,
                                               const FieldM<Complex, 1>& pfdist)
{
  TIMER_VERBOSE("rescale_field_with_psel_fsel_distribution");
  const Geometry& geo = f.geo;
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

inline std::array<SpinMatrix, 8>& get_va_matrices()
{
  static bool initialized = false;
  static std::array<SpinMatrix, 8> ms;
  if (not initialized) {
    initialized = true;
    const std::array<SpinMatrix, 16>& gms = SpinMatrixConstants::get_cps_gms();
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

inline void contract_meson_vv_unshifted_acc_x(
    Vector<Complex> v1, Vector<Complex> v2, const Complex coef,
    const WallSrcProps& wsp1, const WallSrcProps& wsp2,
    const WilsonMatrix& wm3_x_y, const WilsonMatrix& wm3_y_x,
    const Coordinate& xg_x, const long xg_x_idx, const Coordinate& xg_y,
    const long xg_y_psel_idx, const int t_wall, const bool exact)
{
  const std::array<SpinMatrix, 8>& va_ms = get_va_matrices();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  if (exact) {
    qassert(wsp1.exact_tslice_mask[t_wall]);
    qassert(wsp2.exact_tslice_mask[t_wall]);
  }
  const WilsonMatrix& wm1_x_tsrc =
      get_prop(wsp1, t_wall, exact).get_elem(xg_x_idx);
  const WilsonMatrix& wm2_x_tsrc =
      get_prop(wsp2, t_wall, exact).get_elem(xg_x_idx);
  const WilsonMatrix& wm1_y_tsrc =
      get_psel_prop(wsp1, t_wall, exact).get_elem(xg_y_psel_idx);
  const WilsonMatrix& wm2_y_tsrc =
      get_psel_prop(wsp2, t_wall, exact).get_elem(xg_y_psel_idx);
  const WilsonMatrix wm1_tsrc_x =
      gamma5 * (WilsonMatrix)matrix_adjoint(wm1_x_tsrc) * gamma5;
  const WilsonMatrix wm1_tsrc_y =
      gamma5 * (WilsonMatrix)matrix_adjoint(wm1_y_tsrc) * gamma5;
  const WilsonMatrix wm_y_tsrc_x = wm2_y_tsrc * gamma5 * wm1_tsrc_x;
  const WilsonMatrix wm_x_tsrc_y = wm2_x_tsrc * gamma5 * wm1_tsrc_y;
  for (int mu = 0; mu < 8; ++mu) {
    const WilsonMatrix wm_1 = wm_y_tsrc_x * va_ms[mu] * wm3_x_y;
    const WilsonMatrix wm_2 = wm_x_tsrc_y * va_ms[mu] * wm3_y_x;
    for (int nu = 0; nu < 8; ++nu) {
      const int mu_nu = 8 * mu + nu;
      v1[mu_nu] += coef * matrix_trace(wm_1, va_ms[nu]);
      v2[mu_nu] += coef * matrix_trace(wm_2, va_ms[nu]);
    }
  }
}

inline void contract_meson_vv_unshifted_acc_x(
    Vector<Complex> v1, Vector<Complex> v2, const Complex coef,
    const WallSrcProps& wsp1, const WallSrcProps& wsp2,
    const WilsonMatrix& wm3_x_y, const WilsonMatrix& wm3_y_x,
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
    contract_meson_vv_unshifted_acc_x(v1, v2, coef1, wsp1, wsp2, wm3_x_y,
                                      wm3_y_x, xg_x, xg_x_idx, xg_y,
                                      xg_y_psel_idx, t_wall, true);
    contract_meson_vv_unshifted_acc_x(v1, v2, coef2, wsp1, wsp2, wm3_x_y,
                                      wm3_y_x, xg_x, xg_x_idx, xg_y,
                                      xg_y_psel_idx, t_wall, false);
  } else {
    contract_meson_vv_unshifted_acc_x(v1, v2, coef, wsp1, wsp2, wm3_x_y,
                                      wm3_y_x, xg_x, xg_x_idx, xg_y,
                                      xg_y_psel_idx, t_wall, false);
  }
}

inline void contract_meson_vv_unshifted(
    std::vector<SelectedField<Complex> >& meson_vv, const WallSrcProps& wsp1,
    const WallSrcProps& wsp2, const SelProp& prop3_x_y, const Coordinate& xg_y,
    const long xg_y_psel_idx, const int tsep, const PointSelection& psel,
    const FieldSelection& fsel)
// fsel.prob is NOT accounted.
{
  TIMER_VERBOSE("contract_meson_vv_acc_unshifted");
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  qassert(psel[xg_y_psel_idx] == xg_y);
  const Geometry& geo = fsel.f_rank.geo;
  const Coordinate total_site = geo.total_site();
  const int multiplicity = 8 * 8;
  clear(meson_vv);
  meson_vv.resize(4);
  for (int i = 0; i < (int)meson_vv.size(); ++i) {
    meson_vv[i].init(fsel, multiplicity);
    set_zero(meson_vv[i]);
  }
  SelectedField<Complex>& meson_vv_decay_1 = meson_vv[0];
  SelectedField<Complex>& meson_vv_decay_2 = meson_vv[1];
  SelectedField<Complex>& meson_vv_fission_1 = meson_vv[2];
  SelectedField<Complex>& meson_vv_fission_2 = meson_vv[3];
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
    const WilsonMatrix wm3_y_x =
        gamma5 * (WilsonMatrix)matrix_adjoint(wm3_x_y) * gamma5;
    Vector<Complex> vd1 = meson_vv_decay_1.get_elems(idx);
    Vector<Complex> vd2 = meson_vv_decay_2.get_elems(idx);
    Vector<Complex> vf1 = meson_vv_fission_1.get_elems(idx);
    Vector<Complex> vf2 = meson_vv_fission_2.get_elems(idx);
    contract_meson_vv_unshifted_acc_x(vd1, vd2, 1.0, wsp1, wsp2, wm3_x_y,
                                      wm3_y_x, xg_x, xg_x_idx, xg_y,
                                      xg_y_psel_idx, t_src);
    contract_meson_vv_unshifted_acc_x(vf1, vf2, 1.0, wsp1, wsp2, wm3_x_y,
                                      wm3_y_x, xg_x, xg_x_idx, xg_y,
                                      xg_y_psel_idx, t_snk);
  }
  sync_node();
}

inline void contract_meson_vv_acc(
    FieldM<Complex, 8 * 8>& meson_vv_decay,
    FieldM<Complex, 8 * 8>& meson_vv_fission, const WallSrcProps& wsp1,
    const WallSrcProps& wsp2, const SelProp& prop3_x_y, const Coordinate& xg_y,
    const long xg_y_psel_idx, const int tsep, const PointSelection& psel,
    const FieldSelection& fsel, const ShiftShufflePlan& ssp,
    const ShiftShufflePlan& ssp_reflect)
// xg_y = psel[xg_y_psel_idx] is the point src location for prop3_x_y
// ssp = make_shift_shuffle_plan(fsel, -xg_y);
// ssp_reflect = make_shift_shuffle_plan(fsel, -xg_y, true);
{
  TIMER_VERBOSE("contract_meson_vv_acc");
  const Geometry& geo = fsel.f_rank.geo;
  qassert(is_initialized(meson_vv_decay));
  qassert(is_initialized(meson_vv_fission));
  qassert(geo == prop3_x_y.geo);
  qassert(fsel.n_elems == prop3_x_y.n_elems);
  qassert(is_initialized(wsp1));
  qassert(is_initialized(wsp2));
  qassert(is_initialized(prop3_x_y));
  qassert(psel[xg_y_psel_idx] == xg_y);
  qassert(ssp.shift == -xg_y);
  qassert(ssp_reflect.shift == -xg_y);
  qassert(ssp.is_reflect == false);
  qassert(ssp_reflect.is_reflect == true);
  std::vector<SelectedField<Complex> > meson_vv;
  contract_meson_vv_unshifted(meson_vv, wsp1, wsp2, prop3_x_y, xg_y,
                              xg_y_psel_idx, tsep, psel, fsel);
  qassert(meson_vv.size() == 4);
  SelectedField<Complex>& meson_vv_decay_1 = meson_vv[0];
  SelectedField<Complex>& meson_vv_decay_2 = meson_vv[1];
  SelectedField<Complex>& meson_vv_fission_1 = meson_vv[2];
  SelectedField<Complex>& meson_vv_fission_2 = meson_vv[3];
  field_shift(meson_vv_decay_1, meson_vv_decay_1, ssp);
  field_shift(meson_vv_decay_2, meson_vv_decay_2, ssp_reflect);
  field_shift(meson_vv_fission_1, meson_vv_fission_1, ssp);
  field_shift(meson_vv_fission_2, meson_vv_fission_2, ssp_reflect);
  qassert(fsel.prob == ssp.fsel.prob);
  qassert(fsel.prob == ssp_reflect.fsel.prob);
  const Complex coef = 1.0 / 2.0 / fsel.prob;
  acc_field(meson_vv_decay, coef, meson_vv_decay_1, ssp.fsel);
  acc_field(meson_vv_decay, coef, meson_vv_decay_2, ssp_reflect.fsel);
  acc_field(meson_vv_fission, coef, meson_vv_fission_1, ssp.fsel);
  acc_field(meson_vv_fission, coef, meson_vv_fission_2, ssp_reflect.fsel);
}

// -----------------------------------------------------------------------------------

inline void contract_meson_vv_meson_unshifted_acc_x(
    Vector<Complex> v1, Vector<Complex> v2, const Complex coef,
    const WallSrcProps& wsp1, const WallSrcProps& wsp2,
    const WallSrcProps& wsp3, const WilsonMatrix& wm4_x_y,
    const WilsonMatrix& wm4_y_x, const Coordinate& xg_x, const long xg_x_idx,
    const Coordinate& xg_y, const long xg_y_psel_idx, const int t_wall_snk,
    const bool exact_snk, const int t_wall_src, const bool exact_src)
{
  const std::array<SpinMatrix, 8>& va_ms = get_va_matrices();
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
  const WilsonMatrix& wm2_x_tsrc =
      get_prop(wsp2, t_wall_src, exact_src).get_elem(xg_x_idx);
  const WilsonMatrix& wm1_y_tsnk=
      get_psel_prop(wsp1, t_wall_snk, exact_snk).get_elem(xg_y_psel_idx);
  const WilsonMatrix& wm2_y_tsrc =
      get_psel_prop(wsp2, t_wall_src, exact_src).get_elem(xg_y_psel_idx);
  const WilsonMatrix wm1_tsnk_x =
      gamma5 * (WilsonMatrix)matrix_adjoint(wm1_x_tsnk) * gamma5;
  const WilsonMatrix wm1_tsnk_y =
      gamma5 * (WilsonMatrix)matrix_adjoint(wm1_y_tsnk) * gamma5;
  const WilsonMatrix wm3_tsrc_tsnk =
      gamma5 * get_wsnk_prop(wsp3, t_wall_snk, exact_snk)[t_wall_src] * gamma5;
  const WilsonMatrix wm_y_tsrc_tsnk_x = wm2_y_tsrc * wm3_tsrc_tsnk * wm1_tsnk_x;
  const WilsonMatrix wm_x_tsrc_tsnk_y = wm2_x_tsrc * wm3_tsrc_tsnk * wm1_tsnk_y;
  for (int mu = 0; mu < 8; ++mu) {
    const WilsonMatrix wm_1 = wm_y_tsrc_tsnk_x * va_ms[mu] * wm4_x_y;
    const WilsonMatrix wm_2 = wm_x_tsrc_tsnk_y * va_ms[mu] * wm4_y_x;
    for (int nu = 0; nu < 8; ++nu) {
      const int mu_nu = 8 * mu + nu;
      v1[mu_nu] += coef * matrix_trace(wm_1, va_ms[nu]);
      v2[mu_nu] += coef * matrix_trace(wm_2, va_ms[nu]);
    }
  }
}

inline void contract_meson_vv_meson_unshifted_acc_x(
    Vector<Complex> v1, Vector<Complex> v2, const Complex coef,
    const WallSrcProps& wsp1, const WallSrcProps& wsp2,
    const WallSrcProps& wsp3, const WilsonMatrix& wm4_x_y,
    const WilsonMatrix& wm4_y_x, const Coordinate& xg_x, const long xg_x_idx,
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
    contract_meson_vv_meson_unshifted_acc_x(
        v1, v2, coef1, wsp1, wsp2, wsp3, wm4_x_y, wm4_y_x, xg_x, xg_x_idx, xg_y,
        xg_y_psel_idx, t_wall_snk, true, t_wall_src, true);
    contract_meson_vv_meson_unshifted_acc_x(
        v1, v2, coef2, wsp1, wsp2, wsp3, wm4_x_y, wm4_y_x, xg_x, xg_x_idx, xg_y,
        xg_y_psel_idx, t_wall_snk, false, t_wall_src, false);
  } else if (has_exact_src and has_exact_snk) {
    const double sloppy_exact_ratio_11 = wsp1.sloppy_exact_ratio_11;
    const double sloppy_exact_ratio_1 = wsp1.sloppy_exact_ratio_1;
    qassert(sloppy_exact_ratio_1 == wsp2.sloppy_exact_ratio_1);
    const Complex coef1 = sloppy_exact_ratio_11 * coef;
    const Complex coef2 = (sloppy_exact_ratio_1 - sloppy_exact_ratio_11) * coef;
    const Complex coef3 =
        (1.0 - 2.0 * sloppy_exact_ratio_1 + sloppy_exact_ratio_11) * coef;
    contract_meson_vv_meson_unshifted_acc_x(
        v1, v2, coef1, wsp1, wsp2, wsp3, wm4_x_y, wm4_y_x, xg_x, xg_x_idx, xg_y,
        xg_y_psel_idx, t_wall_snk, true, t_wall_src, true);
    contract_meson_vv_meson_unshifted_acc_x(
        v1, v2, coef2, wsp1, wsp2, wsp3, wm4_x_y, wm4_y_x, xg_x, xg_x_idx, xg_y,
        xg_y_psel_idx, t_wall_snk, false, t_wall_src, true);
    contract_meson_vv_meson_unshifted_acc_x(
        v1, v2, coef2, wsp1, wsp2, wsp3, wm4_x_y, wm4_y_x, xg_x, xg_x_idx, xg_y,
        xg_y_psel_idx, t_wall_snk, true, t_wall_src, false);
    contract_meson_vv_meson_unshifted_acc_x(
        v1, v2, coef3, wsp1, wsp2, wsp3, wm4_x_y, wm4_y_x, xg_x, xg_x_idx, xg_y,
        xg_y_psel_idx, t_wall_snk, false, t_wall_src, false);
  } else if (has_exact_snk or has_exact_src) {
    const double sloppy_exact_ratio_1 = wsp1.sloppy_exact_ratio_1;
    qassert(sloppy_exact_ratio_1 == wsp2.sloppy_exact_ratio_1);
    const Complex coef1 = sloppy_exact_ratio_1 * coef;
    const Complex coef2 = (1.0 - sloppy_exact_ratio_1) * coef;
    if (has_exact_snk and (not has_exact_src)) {
      contract_meson_vv_meson_unshifted_acc_x(
          v1, v2, coef1, wsp1, wsp2, wsp3, wm4_x_y, wm4_y_x, xg_x, xg_x_idx,
          xg_y, xg_y_psel_idx, t_wall_snk, true, t_wall_src, false);
    } else if ((not has_exact_snk) and has_exact_src) {
      contract_meson_vv_meson_unshifted_acc_x(
          v1, v2, coef1, wsp1, wsp2, wsp3, wm4_x_y, wm4_y_x, xg_x, xg_x_idx,
          xg_y, xg_y_psel_idx, t_wall_snk, false, t_wall_src, true);
    } else {
      qassert(false);
    }
    contract_meson_vv_meson_unshifted_acc_x(
        v1, v2, coef2, wsp1, wsp2, wsp3, wm4_x_y, wm4_y_x, xg_x, xg_x_idx, xg_y,
        xg_y_psel_idx, t_wall_snk, false, t_wall_src, false);
  } else if ((not has_exact_snk) and (not has_exact_src)) {
    contract_meson_vv_meson_unshifted_acc_x(
        v1, v2, coef, wsp1, wsp2, wsp3, wm4_x_y, wm4_y_x, xg_x, xg_x_idx, xg_y,
        xg_y_psel_idx, t_wall_snk, false, t_wall_src, false);
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
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  qassert(psel[xg_y_psel_idx] == xg_y);
  const Geometry& geo = fsel.f_rank.geo;
  const Coordinate total_site = geo.total_site();
  const int multiplicity = 8 * 8;
  clear(meson_vv_meson);
  meson_vv_meson.resize(4);
  for (int i = 0; i < (int)meson_vv_meson.size(); ++i) {
    meson_vv_meson[i].init(fsel, multiplicity);
    set_zero(meson_vv_meson[i]);
  }
  SelectedField<Complex>& meson_vv_meson_forward_1 = meson_vv_meson[0];
  SelectedField<Complex>& meson_vv_meson_forward_2 = meson_vv_meson[1];
  SelectedField<Complex>& meson_vv_meson_backward_1 = meson_vv_meson[2];
  SelectedField<Complex>& meson_vv_meson_backward_2 = meson_vv_meson[3];
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
    const WilsonMatrix wm4_y_x =
        gamma5 * (WilsonMatrix)matrix_adjoint(wm4_x_y) * gamma5;
    Vector<Complex> vf1 = meson_vv_meson_forward_1.get_elems(idx);
    Vector<Complex> vf2 = meson_vv_meson_forward_2.get_elems(idx);
    Vector<Complex> vb1 = meson_vv_meson_backward_1.get_elems(idx);
    Vector<Complex> vb2 = meson_vv_meson_backward_2.get_elems(idx);
    contract_meson_vv_meson_unshifted_acc_x(vf1, vf2, 1.0, wsp1, wsp2, wsp3,
                                            wm4_x_y, wm4_y_x, xg_x, xg_x_idx,
                                            xg_y, xg_y_psel_idx, t_snk, t_src);
    contract_meson_vv_meson_unshifted_acc_x(vb1, vb2, 1.0, wsp1, wsp2, wsp3,
                                            wm4_x_y, wm4_y_x, xg_x, xg_x_idx,
                                            xg_y, xg_y_psel_idx, t_src, t_snk);
  }
  sync_node();
}

inline void contract_meson_vv_meson_acc(
    FieldM<Complex, 8 * 8>& meson_vv_meson_forward,
    FieldM<Complex, 8 * 8>& meson_vv_meson_backward, const WallSrcProps& wsp1,
    const WallSrcProps& wsp2, const WallSrcProps& wsp3,
    const SelProp& prop4_x_y, const Coordinate& xg_y, const long xg_y_psel_idx,
    const int tsep, const PointSelection& psel, const FieldSelection& fsel,
    const ShiftShufflePlan& ssp, const ShiftShufflePlan& ssp_reflect)
// xg_y = psel[xg_y_psel_idx] is the point src location for prop3_x_y
// ssp = make_shift_shuffle_plan(fsel, -xg_y);
// ssp_reflect = make_shift_shuffle_plan(fsel, -xg_y, true);
{
  TIMER_VERBOSE("contract_meson_vv_meson_acc");
  const Geometry& geo = fsel.f_rank.geo;
  qassert(is_initialized(meson_vv_meson_forward));
  qassert(is_initialized(meson_vv_meson_backward));
  qassert(geo == prop4_x_y.geo);
  qassert(fsel.n_elems == prop4_x_y.n_elems);
  qassert(is_initialized(wsp1));
  qassert(is_initialized(wsp2));
  qassert(is_initialized(wsp3));
  qassert(is_initialized(prop4_x_y));
  qassert(psel[xg_y_psel_idx] == xg_y);
  qassert(ssp.shift == -xg_y);
  qassert(ssp_reflect.shift == -xg_y);
  qassert(ssp.is_reflect == false);
  qassert(ssp_reflect.is_reflect == true);
  std::vector<SelectedField<Complex> > meson_vv_meson;
  contract_meson_vv_meson_unshifted(meson_vv_meson, wsp1, wsp2, wsp3, prop4_x_y,
                                    xg_y, xg_y_psel_idx, tsep, psel, fsel);
  qassert(meson_vv_meson.size() == 4);
  SelectedField<Complex>& meson_vv_meson_forward_1 = meson_vv_meson[0];
  SelectedField<Complex>& meson_vv_meson_forward_2 = meson_vv_meson[1];
  SelectedField<Complex>& meson_vv_meson_backward_1 = meson_vv_meson[2];
  SelectedField<Complex>& meson_vv_meson_backward_2 = meson_vv_meson[3];
  field_shift(meson_vv_meson_forward_1, meson_vv_meson_forward_1, ssp);
  field_shift(meson_vv_meson_forward_2, meson_vv_meson_forward_2, ssp_reflect);
  field_shift(meson_vv_meson_backward_1, meson_vv_meson_backward_1, ssp);
  field_shift(meson_vv_meson_backward_2, meson_vv_meson_backward_2,
              ssp_reflect);
  qassert(fsel.prob == ssp.fsel.prob);
  qassert(fsel.prob == ssp_reflect.fsel.prob);
  const Complex coef = 1.0 / 2.0 / fsel.prob;
  acc_field(meson_vv_meson_forward, coef, meson_vv_meson_forward_1, ssp.fsel);
  acc_field(meson_vv_meson_forward, coef, meson_vv_meson_forward_2,
            ssp_reflect.fsel);
  acc_field(meson_vv_meson_backward, coef, meson_vv_meson_backward_1, ssp.fsel);
  acc_field(meson_vv_meson_backward, coef, meson_vv_meson_backward_2,
            ssp_reflect.fsel);
}

}  // namespace qlat
