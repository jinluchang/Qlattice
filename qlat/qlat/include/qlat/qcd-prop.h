#pragma once

#include <qlat/qcd.h>
#include <qlat/selected-field.h>
#include <qlat/selected-points.h>
#include <qlat/field-fft.h>

namespace qlat
{  //

template <class T>
WilsonMatrixT<T> make_wilson_matrix_from_vectors(
    const array<ConstHandle<WilsonVectorT<T>>, 4 * NUM_COLOR>& cols)
{
  WilsonMatrixT<T> ret;
  for (Int i = 0; i < 4 * NUM_COLOR; ++i) {
    for (Int j = 0; j < 4 * NUM_COLOR; ++j) {
      ret(j, i) = cols[i]()(j);
    }
  }
  return ret;
}

template <class T>
void set_propagator_from_fermion_fields(
    Propagator4dT<T>& prop, const Array<FermionField4dT<T>, 4 * NUM_COLOR> ffs)
{
  TIMER_VERBOSE("set_propagator_from_fermion_fields");
  const Geometry& geo = ffs[0].geo();
  for (Int i = 0; i < 4 * NUM_COLOR; ++i) {
    qassert(geo == ffs[i].geo());
  }
  prop.init(geo_resize(geo));
  qassert(prop.geo() == geo);
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    array<ConstHandle<WilsonVectorT<T>>, 4 * NUM_COLOR> cols;
    for (Int k = 0; k < 4 * NUM_COLOR; ++k) {
      cols[k].init(ffs[k].get_elem(xl));
    }
    prop.get_elem(xl) = make_wilson_matrix_from_vectors(cols);
  }
}

template <class T>
qacc void set_wilson_matrix_col_from_vector(WilsonMatrixT<T>& wm, const Int idx,
                                            const WilsonVectorT<T>& col)
{
  for (Int j = 0; j < 4 * NUM_COLOR; ++j) {
    wm(j, idx) = col(j);
  }
}

template <class T>
void set_propagator_col_from_fermion_field(Propagator4dT<T>& prop,
                                           const Int idx,
                                           const FermionField4dT<T>& ff)
{
  TIMER("set_propagator_col_from_fermion_field");
  const Geometry& geo = ff.geo();
  qassert(geo == prop.geo());
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    set_wilson_matrix_col_from_vector(prop.get_elem(xl), idx, ff.get_elem(xl));
  }
}

template <class T>
qacc void set_wilson_vector_from_matrix_col(WilsonVectorT<T>& col,
                                            const WilsonMatrixT<T>& wm,
                                            const Int idx)
{
  for (Int j = 0; j < 4 * NUM_COLOR; ++j) {
    col(j) = wm(j, idx);
  }
}

template <class T>
void set_fermion_field_from_propagator_col(FermionField4dT<T>& ff,
                                           const Propagator4dT<T>& prop,
                                           const Int idx)
{
  TIMER("set_fermion_field_from_propagator_col");
  const Geometry& geo = prop.geo();
  ff.init(geo_resize(geo));
  qassert(geo == ff.geo());
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    set_wilson_vector_from_matrix_col(ff.get_elem(xl), prop.get_elem(xl), idx);
  }
}

template <class T>
void fermion_field_5d_from_4d(FermionField5dT<T>& ff5d,
                              const FermionField4dT<T>& ff4d, const Int upper,
                              const Int lower)
// ff5d need to be initialized
// upper componets are right handed
// lower componets are left handed
{
  TIMER("fermion_field_5d_from_4d");
  const Geometry& geo = ff5d.geo();
  set_zero(ff5d);
  const Int sizewvh = sizeof(WilsonVectorT<T>) / 2;
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    Coordinate x = geo.coordinate_from_index(index);
    memcpy((char*)&(ff5d.get_elem(x, upper)), (const char*)&(ff4d.get_elem(x)),
           sizewvh);
    memcpy((char*)&(ff5d.get_elem(x, lower)) + sizewvh,
           (const char*)&(ff4d.get_elem(x)) + sizewvh, sizewvh);
  }
}

template <class T>
void fermion_field_4d_from_5d(FermionField4dT<T>& ff4d,
                              const FermionField5dT<T>& ff5d, const Int upper,
                              const Int lower)
// upper componets are right handed
// lower componets are left handed
{
  TIMER("fermion_field_4d_from_5d");
  const Geometry& geo = ff5d.geo();
  ff4d.init(geo_resize(geo));
  set_zero(ff4d);
  const Int sizewvh = sizeof(WilsonVectorT<T>) / 2;
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    Coordinate x = geo.coordinate_from_index(index);
    memcpy((char*)&(ff4d.get_elem(x)), (const char*)&(ff5d.get_elem(x, upper)),
           sizewvh);
    memcpy((char*)&(ff4d.get_elem(x)) + sizewvh,
           (const char*)&(ff5d.get_elem(x, lower)) + sizewvh, sizewvh);
  }
}

void set_ff_vec_from_prop(std::vector<FermionField4d>& ff_vec,
                          const Propagator4d& prop);

void set_prop_from_ff_vec(Propagator4d& prop,
                          const std::vector<FermionField4d>& ff_vec);

template <class Inverter, class T>
Long invert_dwf(FermionField4dT<T>& sol, const FermionField4dT<T>& src,
                const Inverter& inv, const Int ls_ = 0)
// sol do not need to be initialized
// inv.geo() must be the geometry of the fermion field
// invert(sol5d, src5d, inv) perform the inversion
{
  TIMER_VERBOSE("invert_dwf(4d,4d,inv)");
  const Geometry& geo = src.geo();
  qassert(check_matching_geo(geo, inv.geo()));
  const Int ls = ls_ != 0 ? ls_ : inv.fa.ls;
  const Geometry geo_ls = geo_resize(inv.geo());
  FermionField5dT<T> sol5d, src5d;
  sol5d.init(geo_ls, ls);
  src5d.init(geo_ls, ls);
  fermion_field_5d_from_4d(src5d, src, 0, ls - 1);
  set_zero(sol5d);
  const Long iter = invert(sol5d, src5d, inv);
  sol.init(geo, src.multiplicity);
  fermion_field_4d_from_5d(sol, sol5d, ls - 1, 0);
  return iter;
}

template <class Inverter, class T>
void invert(Propagator4dT<T>& sol, const Propagator4dT<T>& src,
            const Inverter& inv)
// sol do not need to be initialized
// inv.geo() must be the geometry of the fermion field
// invert(4d, 4d, inv) perform the inversion
{
  TIMER_VERBOSE("invert(p4d,p4d,inv)");
  qassert(&sol != &src);
  const Geometry geo = geo_resize(src.geo());
  sol.init(geo);
  FermionField4dT<T> ff_sol, ff_src;
  for (Int j = 0; j < 4 * NUM_COLOR; ++j) {
    set_fermion_field_from_propagator_col(ff_src, src, j);
    invert(ff_sol, ff_src, inv);
    set_propagator_col_from_fermion_field(sol, j, ff_sol);
  }
}

void set_point_src_fermion_field(FermionField4dT<RealD>& ff,
                                 const Coordinate& xg, const Int cs,
                                 const ComplexD& value = 1.0);

void set_point_src_fermion_field(FermionField4dT<RealF>& ff,
                                 const Coordinate& xg, const Int cs,
                                 const ComplexD& value = 1.0);

void set_point_src(Propagator4dT<RealD>& prop, const Geometry& geo_input,
                   const Coordinate& xg, const ComplexD& value = 1.0);

void set_point_src(Propagator4dT<RealF>& prop, const Geometry& geo_input,
                   const Coordinate& xg, const ComplexD& value = 1.0);

template <class Inverter, class T>
void set_point_src_propagator(Propagator4dT<T>& prop, const Inverter& inv,
                              const Coordinate& xg, const ComplexD& value = 1.0)
{
  TIMER_VERBOSE("set_point_src_propagator");
  const Geometry geo = geo_resize(inv.geo());
  prop.init(geo);
  FermionField4dT<T> sol, src;
  sol.init(geo);
  src.init(geo);
  for (Int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_point_src_fermion_field(src, xg, cs, value);
    set_zero(sol);
    invert(sol, src, inv);
    set_propagator_col_from_fermion_field(prop, cs, sol);
  }
}

void set_wall_src_fermion_field(FermionField4d& ff, const Int tslice,
                                const CoordinateD& lmom, const Int cs);

void set_wall_src(Propagator4d& prop, const Geometry& geo_input,
                  const Int tslice, const CoordinateD& lmom = CoordinateD());

template <class Inverter>
void set_wall_src_propagator(Propagator4d& prop, const Inverter& inv,
                                    const Int tslice,
                                    const CoordinateD& lmom = CoordinateD())
{
  TIMER_VERBOSE("set_wall_src_propagator");
  const Geometry geo = geo_resize(inv.geo());
  prop.init(geo);
  FermionField4d sol, src;
  sol.init(geo);
  src.init(geo);
  for (Int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_wall_src_fermion_field(src, tslice, lmom, cs);
    set_zero(sol);
    invert(sol, src, inv);
    set_propagator_col_from_fermion_field(prop, cs, sol);
  }
}

void set_rand_u1_src_psel(Propagator4d& prop, FieldM<ComplexD, 1>& fu1,
                          const PointsSelection& psel, const Geometry& geo_,
                          const RngState& rs);

void set_rand_u1_sol_psel(SelectedPoints<WilsonMatrix>& sp_prop,
                          const Propagator4d& prop,
                          const FieldM<ComplexD, 1>& fu1,
                          const PointsSelection& psel);

void set_rand_u1_src_fsel(Propagator4d& prop, FieldM<ComplexD, 1>& fu1,
                          const FieldSelection& fsel, const RngState& rs);

void set_rand_u1_sol_fsel(SelectedField<WilsonMatrix>& sf_prop,
                          const Propagator4d& prop,
                          const FieldM<ComplexD, 1>& fu1,
                          const FieldSelection& fsel);

void set_rand_vol_u1(Field<ComplexD>& fu1, const Geometry& geo_input,
                     const RngState& rs);

void set_rand_vol_u1_src(Propagator4d& prop, const Field<ComplexD>& fu1);

void set_mom_src_fermion_field(FermionField4d& ff, const CoordinateD& lmom,
                               const Int cs);

void set_mom_src(Propagator4d& prop, const Geometry& geo_input,
                 const CoordinateD& lmom);

template <class Inverter>
void set_mom_src_propagator(Propagator4d& prop, Inverter& inv,
                            const CoordinateD& lmom)
{
  TIMER_VERBOSE("set_mom_src_propagator");
  const Geometry& geo = geo_remult(inv.geo());
  prop.init(geo);
  qassert(prop.geo() == geo);
  FermionField4d src, sol;
  src.init(geo);
  sol.init(geo);
  for (Int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_mom_src_fermion_field(src, lmom, cs);
    set_zero(sol);
    invert(sol, src, inv);
    set_propagator_col_from_fermion_field(prop, cs, sol);
  }
}

template <class T>
void free_mom_invert(Propagator4dT<T>& sol, const Propagator4dT<T>& src,
                     const RealD mass, const RealD m5 = 1.0,
                     const CoordinateD& momtwist = CoordinateD())
// DWF infinite L_s
// M_5 <= 1.0
{
  TIMER("free_mom_invert");
  sol.init(src);
  const Geometry& geo = src.geo();
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    Coordinate kl = geo.coordinate_from_index(index);
    Coordinate kg = geo.coordinate_g_from_l(kl);
    array<double, DIMN> kk, ks;
    double p2 = 0.0;
    double wp = 1.0 - m5;
    SpinMatrixT<T> pg;
    set_zero(pg);
    for (Int i = 0; i < DIMN; ++i) {
      Coordinate total_site = geo.total_site();
      kg[i] = smod(kg[i], total_site[i]);
      kk[i] = 2.0 * PI * (kg[i] + momtwist[i]) / (double)total_site[i];
      ks[i] = sin(kk[i]);
      pg += SpinMatrixConstantsT<T>::get_cps_gammas()[i] * (ComplexT<T>)ks[i];
      p2 += sqr(ks[i]);
      wp += 2.0 * sqr(sin(kk[i] / 2.0));
    }
    const double calpha = (1.0 + sqr(wp) + p2) / 2.0 / wp;
    const double alpha = acosh(calpha);
    const double lwa = 1.0 - wp * exp(-alpha);
    SpinMatrixT<T> m;
    set_unit(m, mass * lwa);
    SpinMatrixT<T> ipgm = pg;
    ipgm *= (ComplexT<T>)(-ii);
    ipgm += m;
    ipgm *= lwa / (p2 + sqr(mass * lwa));
    WilsonMatrixT<T>& wm_sol = sol.get_elem(kl);
    if (1.0e-10 > p2 && 1.0e-10 > lwa) {
      // if (0.0 != qnorm(ipgm)) {
      //   Display(cname, fname, "kg = %s\n", show(kg).c_str());
      //   Display(cname, fname, "p2         = %13.5E\n", p2);
      //   Display(cname, fname, "wp         = %13.5E\n", wp);
      //   Display(cname, fname, "alpha      = %13.5E\n", alpha);
      //   Display(cname, fname, "lwa        = %13.5E\n", lwa);
      //   Display(cname, fname, "qnorm(ipgm) = %13.5E\n", qnorm(ipgm));
      // }
      set_zero(wm_sol);
    } else {
      const WilsonMatrixT<T>& wm_src = src.get_elem(kl);
      wm_sol = ipgm * wm_src;
    }
  }
}

template <class T>
void free_invert(Propagator4dT<T>& sol, const Propagator4dT<T>& src,
                 const RealD mass, const RealD m5 = 1.0,
                 const CoordinateD& momtwist = CoordinateD())
{
  TIMER_VERBOSE("free_invert");
  const Geometry& geo = src.geo();
  sol.init(src);
  fft_complex_field(sol, true);
  free_mom_invert(sol, sol, mass, m5, momtwist);
  fft_complex_field(sol, false);
  sol *= 1.0 / geo.total_volume();
}

void free_invert(Prop& p_sol, const Prop& p_src, const RealD mass, const RealD m5 = 1.0,
                 const CoordinateD& momtwist = CoordinateD());

void convert_wm_from_mspincolor(Propagator4d& prop_wm,
                                const Propagator4d& prop_msc);

void convert_mspincolor_from_wm(Propagator4d& prop_msc,
                                const Propagator4d& prop_wm);

void convert_wm_from_mspincolor(SelectedField<WilsonMatrix>& prop_wm,
                                const SelectedField<WilsonMatrix>& prop_msc);

void convert_mspincolor_from_wm(SelectedField<WilsonMatrix>& prop_msc,
                                const SelectedField<WilsonMatrix>& prop_wm);

void convert_wm_from_mspincolor(SelectedPoints<WilsonMatrix>& prop_wm,
                                const SelectedPoints<WilsonMatrix>& prop_msc);

void convert_mspincolor_from_wm(SelectedPoints<WilsonMatrix>& prop_msc,
                                const SelectedPoints<WilsonMatrix>& prop_wm);

void set_t_range_flip_tpbc_with_tslice(int& t_start, Int& t_stop,
                                       const Int tslice_flip_tpbc,
                                       const Int t_size);

void flip_tpbc_with_tslice(SelectedPoints<WilsonMatrix>& ps_prop,
                           const PointsSelection& psel,
                           const Int tslice_flip_tpbc, const Int t_size);

void flip_tpbc_with_tslice(SelectedField<WilsonMatrix>& s_prop,
                           const FieldSelection& fsel,
                           const Int tslice_flip_tpbc);

// -------------------------------------------------------------------------

void set_tslice_mom_src_fermion_field(FermionField4d& ff, const Int tslice,
                                      const CoordinateD& lmom, const Int cs);

void set_tslice_mom_src(Propagator4d& prop, const Geometry& geo_input,
                        const Int tslice, const CoordinateD& lmom);

template <class Inverter>
void set_tslice_mom_src_propagator(Propagator4d& prop, const Int tslice,
                                   const CoordinateD& lmom, Inverter& inverter)
{
  TIMER_VERBOSE("set_tslice_mom_src_propagator");
  const Geometry& geo = geo_remult(inverter.geo());
  prop.init(geo);
  qassert(prop.geo() == geo);
  FermionField4d src, sol;
  src.init(geo);
  sol.init(geo);
  for (Int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_tslice_mom_src_fermion_field(src, tslice, lmom, cs);
    set_zero(sol);
    invert(sol, src, inverter);
    set_propagator_col_from_fermion_field(prop, cs, sol);
  }
}

void set_volume_src_fermion_field(FermionField4d& ff, const CoordinateD& lmom,
                                  const Int cs);

void set_volume_src(Propagator4d& prop, const Geometry& geo_input,
                    const CoordinateD& lmom = CoordinateD());

template <class Inverter>
void set_volume_src_propagator(Propagator4d& prop, const Inverter& inv,
                               const CoordinateD& lmom = CoordinateD())
{
  TIMER_VERBOSE("set_volume_src_propagator");
  const Geometry geo = geo_resize(inv.geo());
  prop.init(geo);
  FermionField4d sol, src;
  sol.init(geo);
  src.init(geo);
  for (Int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_volume_src_fermion_field(src, lmom, cs);
    set_zero(sol);
    invert(sol, src, inv);
    set_propagator_col_from_fermion_field(prop, cs, sol);
  }
}

}  // namespace qlat
