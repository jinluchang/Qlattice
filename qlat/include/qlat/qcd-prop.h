#pragma once

#include <qlat/qcd.h>
#include <qlat/selected-field.h>
#include <qlat/selected-points.h>

namespace qlat
{  //

template <class T>
WilsonMatrixT<T> make_wilson_matrix_from_vectors(
    const array<ConstHandle<WilsonVectorT<T> >, 4 * NUM_COLOR>& cols)
{
  WilsonMatrixT<T> ret;
  for (int i = 0; i < 4 * NUM_COLOR; ++i) {
    for (int j = 0; j < 4 * NUM_COLOR; ++j) {
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
  for (int i = 0; i < 4 * NUM_COLOR; ++i) {
    qassert(geo == ffs[i].geo());
  }
  prop.init(geo_reform(geo));
  qassert(prop.geo() == geo);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    array<ConstHandle<WilsonVectorT<T> >, 4 * NUM_COLOR> cols;
    for (int k = 0; k < 4 * NUM_COLOR; ++k) {
      cols[k].init(ffs[k].get_elem(xl));
    }
    prop.get_elem(xl) = make_wilson_matrix_from_vectors(cols);
  }
}

template <class T>
inline void set_wilson_matrix_col_from_vector(WilsonMatrixT<T>& wm,
                                              const int idx,
                                              const WilsonVectorT<T>& col)
{
  for (int j = 0; j < 4 * NUM_COLOR; ++j) {
    wm(j, idx) = col(j);
  }
}

template <class T>
inline void set_propagator_col_from_fermion_field(Propagator4dT<T>& prop,
                                                  const int idx,
                                                  const FermionField4dT<T>& ff)
{
  TIMER("set_propagator_col_from_fermion_field");
  const Geometry& geo = ff.geo();
  qassert(geo == prop.geo());
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    set_wilson_matrix_col_from_vector(prop.get_elem(xl), idx, ff.get_elem(xl));
  }
}

template <class T>
inline void set_wilson_vector_from_matrix_col(WilsonVectorT<T>& col,
                                              const WilsonMatrixT<T>& wm,
                                              const int idx)
{
  for (int j = 0; j < 4 * NUM_COLOR; ++j) {
    col(j) = wm(j, idx);
  }
}

template <class T>
inline void set_fermion_field_from_propagator_col(FermionField4dT<T>& ff,
                                                  const Propagator4dT<T>& prop,
                                                  const int idx)
{
  TIMER("set_fermion_field_from_propagator_col");
  const Geometry& geo = prop.geo();
  ff.init(geo_reform(geo));
  qassert(geo == ff.geo());
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    set_wilson_vector_from_matrix_col(ff.get_elem(xl), prop.get_elem(xl), idx);
  }
}

template <class T>
inline void fermion_field_5d_from_4d(FermionField5dT<T>& ff5d,
                                     const FermionField4dT<T>& ff4d,
                                     const int upper, const int lower)
// ff5d need to be initialized
// upper componets are right handed
// lower componets are left handed
{
  TIMER("fermion_field_5d_from_4d");
  const Geometry& geo = ff5d.geo();
  set_zero(ff5d);
  const int sizewvh = sizeof(WilsonVectorT<T>) / 2;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate x = geo.coordinate_from_index(index);
    memcpy((char*)&(ff5d.get_elem(x, upper)), (const char*)&(ff4d.get_elem(x)),
           sizewvh);
    memcpy((char*)&(ff5d.get_elem(x, lower)) + sizewvh,
           (const char*)&(ff4d.get_elem(x)) + sizewvh, sizewvh);
  }
}

template <class T>
void fermion_field_4d_from_5d(FermionField4dT<T>& ff4d,
                              const FermionField5dT<T>& ff5d, const int upper,
                              const int lower)
// upper componets are right handed
// lower componets are left handed
{
  TIMER("fermion_field_4d_from_5d");
  const Geometry& geo = ff5d.geo();
  ff4d.init(geo_reform(geo));
  set_zero(ff4d);
  const int sizewvh = sizeof(WilsonVectorT<T>) / 2;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate x = geo.coordinate_from_index(index);
    memcpy((char*)&(ff4d.get_elem(x)), (const char*)&(ff5d.get_elem(x, upper)),
           sizewvh);
    memcpy((char*)&(ff4d.get_elem(x)) + sizewvh,
           (const char*)&(ff5d.get_elem(x, lower)) + sizewvh, sizewvh);
  }
}

template <class Inverter, class T>
inline long invert_dwf(FermionField4dT<T>& sol, const FermionField4dT<T>& src,
                       const Inverter& inv, const int ls_ = 0)
// sol do not need to be initialized
// inv.geo() must be the geometry of the fermion field
// invert(sol5d, src5d, inv) perform the inversion
{
  TIMER_VERBOSE("invert_dwf(4d,4d,inv)");
  const Geometry& geo = src.geo();
  qassert(check_matching_geo(geo, inv.geo()));
  sol.init(geo);
  const int ls = ls_ != 0 ? ls_ : inv.fa.ls;
  const Geometry geo_ls = geo_reform(inv.geo(), ls, 0);
  FermionField5dT<T> sol5d, src5d;
  sol5d.init(geo_ls);
  src5d.init(geo_ls);
  fermion_field_5d_from_4d(src5d, src, 0, ls - 1);
  fermion_field_5d_from_4d(sol5d, sol, ls - 1, 0);
  const long iter = invert(sol5d, src5d, inv);
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
  const Geometry geo = geo_reform(src.geo());
  sol.init(geo);
  FermionField4dT<T> ff_sol, ff_src;
  for (int j = 0; j < 4 * NUM_COLOR; ++j) {
    set_fermion_field_from_propagator_col(ff_src, src, j);
    invert(ff_sol, ff_src, inv);
    set_propagator_col_from_fermion_field(sol, j, ff_sol);
  }
}

template <class T>
void set_point_src_fermion_field(FermionField4dT<T>& ff, const Coordinate& xg,
                                 const int cs, const Complex& value = 1.0)
// ff need to be initialized
{
  TIMER("set_point_src_fermion_field");
  const Geometry& geo = ff.geo();
  set_zero(ff);
  const Coordinate xl = geo.coordinate_l_from_g(xg);
  if (geo.is_local(xl)) {
    ff.get_elem(xl)(cs) = value;
  }
}

template <class T>
void set_point_src(Propagator4dT<T>& prop, const Geometry& geo_input,
                   const Coordinate& xg, const Complex& value = 1.0)
{
  TIMER_VERBOSE("set_point_src");
  const Geometry geo = geo_reform(geo_input);
  prop.init(geo);
  FermionField4dT<T> src;
  src.init(geo);
  for (int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_point_src_fermion_field(src, xg, cs, value);
    set_propagator_col_from_fermion_field(prop, cs, src);
  }
}

template <class Inverter, class T>
void set_point_src_propagator(Propagator4dT<T>& prop, const Inverter& inv,
                              const Coordinate& xg, const Complex& value = 1.0)
{
  TIMER_VERBOSE("set_point_src_propagator");
  const Geometry geo = geo_reform(inv.geo());
  prop.init(geo);
  FermionField4dT<T> sol, src;
  sol.init(geo);
  src.init(geo);
  for (int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_point_src_fermion_field(src, xg, cs, value);
    set_zero(sol);
    invert(sol, src, inv);
    set_propagator_col_from_fermion_field(prop, cs, sol);
  }
}

inline void set_wall_src_fermion_field(FermionField4d& ff, const int tslice,
                                       const CoordinateD& lmom, const int cs)
// ff need to be initialized beforehand
{
  qassert(lmom[3] == 0);
  const Geometry& geo = ff.geo();
  const CoordinateD mom = lmom * lattice_mom_mult(geo);
  set_zero(ff);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    if (xg[3] == tslice) {
      double phase = 0.0;
      for (int i = 0; i < DIMN; ++i) {
        phase += mom[i] * xg[i];
      }
      ff.get_elem(xl)(cs) = std::polar(1.0, phase);
    }
  }
}

inline void set_wall_src(Propagator4d& prop, const Geometry& geo_input,
                         const int tslice,
                         const CoordinateD& lmom = CoordinateD())
{
  TIMER_VERBOSE("set_wall_src");
  const Geometry geo = geo_reform(geo_input);
  prop.init(geo);
  FermionField4d src;
  src.init(geo);
  for (int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_wall_src_fermion_field(src, tslice, lmom, cs);
    set_propagator_col_from_fermion_field(prop, cs, src);
  }
}

template <class Inverter>
inline void set_wall_src_propagator(Propagator4d& prop, const Inverter& inv,
                                    const int tslice,
                                    const CoordinateD& lmom = CoordinateD())
{
  TIMER_VERBOSE("set_wall_src_propagator");
  const Geometry geo = geo_reform(inv.geo());
  prop.init(geo);
  FermionField4d sol, src;
  sol.init(geo);
  src.init(geo);
  for (int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_wall_src_fermion_field(src, tslice, lmom, cs);
    set_zero(sol);
    invert(sol, src, inv);
    set_propagator_col_from_fermion_field(prop, cs, sol);
  }
}

inline void set_rand_u1_src_psel(Propagator4d& prop, FieldM<Complex, 1>& fu1,
                                 const PointSelection& psel,
                                 const Geometry& geo_, const RngState& rs)
{
  TIMER_VERBOSE("set_rand_u1_src_psel");
  const Geometry geo = geo_reform(geo_);
  const Coordinate total_site = geo.total_site();
  prop.init(geo);
  fu1.init(geo);
  set_zero(prop);
  set_zero(fu1);
  qthread_for(idx, (long)psel.size(), {
    const Coordinate xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      const long gindex = index_from_coordinate(xg, total_site);
      RngState rst = rs.newtype(gindex);
      const double phase = u_rand_gen(rst, PI, -PI);
      const Complex u1 = std::polar(1.0, phase);
      set_unit(prop.get_elem(xl), u1);
      fu1.get_elem(xl) = u1;
    }
  });
}

inline void set_rand_u1_sol_psel(SelectedPoints<WilsonMatrix>& sp_prop,
                                 const Propagator4d& prop,
                                 const FieldM<Complex, 1>& fu1,
                                 const PointSelection& psel)
// calculate self loop at psel locations
{
  TIMER_VERBOSE("set_rand_u1_sol_psel")
  SelectedPoints<Complex> sp_fu1;
  set_selected_points(sp_prop, prop, psel);
  set_selected_points(sp_fu1, fu1, psel);
  qthread_for(idx, (long)psel.size(), {
    const Complex& u1 = sp_fu1.get_elem(idx);
    WilsonMatrix& wm = sp_prop.get_elem(idx);
    wm *= qlat::qconj(u1);
  });
}

inline void set_rand_u1_src_fsel(Propagator4d& prop, FieldM<Complex, 1>& fu1,
                                 const FieldSelection& fsel, const RngState& rs)
{
  TIMER_VERBOSE("set_rand_u1_src_fsel");
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate total_site = geo.total_site();
  prop.init(geo);
  fu1.init(geo);
  set_zero(prop);
  set_zero(fu1);
  qthread_for(idx, fsel.n_elems, {
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    qassert(geo.is_local(xl));
    const long gindex = index_from_coordinate(xg, total_site);
    RngState rst = rs.newtype(gindex);
    const double phase = u_rand_gen(rst, PI, -PI);
    const Complex u1 = std::polar(1.0, phase);
    set_unit(prop.get_elem(xl), u1);
    fu1.get_elem(xl) = u1;
  });
}

inline void set_rand_u1_sol_fsel(SelectedField<WilsonMatrix>& sf_prop,
                                 const Propagator4d& prop,
                                 const FieldM<Complex, 1>& fu1,
                                 const FieldSelection& fsel)
// calculate self loop at fsel locations
{
  TIMER_VERBOSE("set_rand_u1_sol_fsel")
  SelectedField<Complex> sf_fu1;
  set_selected_field(sf_prop, prop, fsel);
  set_selected_field(sf_fu1, fu1, fsel);
  qthread_for(idx, fsel.n_elems, {
    const Complex& u1 = sf_fu1.get_elem(idx);
    WilsonMatrix& wm = sf_prop.get_elem(idx);
    wm *= qlat::qconj(u1);
  });
}

inline void set_mom_src_fermion_field(FermionField4d& ff,
                                      const CoordinateD& lmom, const int cs)
// ff need to be initialized beforehand
{
  const Geometry& geo = ff.geo();
  const CoordinateD mom = lmom * lattice_mom_mult(geo);
  set_zero(ff);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    double phase = 0.0;
    for (int i = 0; i < DIMN; ++i) {
      phase += mom[i] * xg[i];
    }
    ff.get_elem(xl)(cs) = std::polar(1.0, phase);
  }
}

inline void set_mom_src(Propagator4d& prop, const Geometry& geo_input,
                        const CoordinateD& lmom)
{
  TIMER_VERBOSE("set_mom_src");
  const Geometry& geo = geo_reform(geo_input);
  prop.init(geo);
  qassert(prop.geo() == geo);
  FermionField4d src;
  src.init(geo);
  for (int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_mom_src_fermion_field(src, lmom, cs);
    set_propagator_col_from_fermion_field(prop, cs, src);
  }
}

template <class Inverter>
inline void set_mom_src_propagator(Propagator4d& prop, Inverter& inv,
                                   const CoordinateD& lmom)
{
  TIMER_VERBOSE("set_mom_src_propagator");
  const Geometry& geo = geo_remult(inv.geo());
  prop.init(geo);
  qassert(prop.geo() == geo);
  FermionField4d src, sol;
  src.init(geo);
  sol.init(geo);
  for (int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_mom_src_fermion_field(src, lmom, cs);
    set_zero(sol);
    invert(sol, src, inv);
    set_propagator_col_from_fermion_field(prop, cs, sol);
  }
}

template <class T>
void free_mom_invert(Propagator4dT<T>& sol, const Propagator4dT<T>& src,
                     const double mass, const double m5 = 1.0,
                     const CoordinateD& momtwist = CoordinateD())
// DWF infinite L_s
// M_5 <= 1.0
{
  TIMER("free_mom_invert");
  sol.init(src);
  const Geometry& geo = src.geo();
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    Coordinate kl = geo.coordinate_from_index(index);
    Coordinate kg = geo.coordinate_g_from_l(kl);
    array<double, DIMN> kk, ks;
    double p2 = 0.0;
    double wp = 1.0 - m5;
    SpinMatrixT<T> pg;
    set_zero(pg);
    for (int i = 0; i < DIMN; ++i) {
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
                 const double mass, const double m5 = 1.0,
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

inline void convert_wm_from_mspincolor(Propagator4d& prop_wm, const Propagator4d& prop_msc)
{
  TIMER("convert_wm_from_mspincolor");
  const Geometry& geo = prop_msc.geo();
  prop_wm.init(geo);
  qassert(geo.is_only_local);
  qassert(prop_wm.geo().is_only_local);
  qacc_for(index, geo.local_volume(), {
    WilsonMatrix& wm = prop_wm.get_elem(index);
    const WilsonMatrix& msc = prop_msc.get_elem(index);
    convert_wm_from_mspincolor(wm, msc);
  });
}

inline void convert_mspincolor_from_wm(Propagator4d& prop_msc, const Propagator4d& prop_wm)
{
  TIMER("convert_mspincolor_from_wm");
  const Geometry& geo = prop_wm.geo();
  prop_msc.init(geo);
  qassert(geo.is_only_local);
  qassert(prop_wm.geo().is_only_local);
  qacc_for(index, geo.local_volume(), {
    WilsonMatrix& msc = prop_msc.get_elem(index);
    const WilsonMatrix& wm = prop_wm.get_elem(index);
    convert_mspincolor_from_wm(msc, wm);
  });
}

inline void convert_wm_from_mspincolor(SelectedField<WilsonMatrix>& prop_wm, const SelectedField<WilsonMatrix>& prop_msc)
{
  TIMER("convert_wm_from_mspincolor(s_prop)");
  qassert(prop_msc.geo().multiplicity == 1);
  prop_wm.init(prop_msc.geo(), prop_msc.n_elems, 1);
  qacc_for(idx, prop_msc.n_elems, {
    WilsonMatrix& wm = prop_wm.get_elem(idx);
    const WilsonMatrix& msc = prop_msc.get_elem(idx);
    convert_wm_from_mspincolor(wm, msc);
  });
}

inline void convert_mspincolor_from_wm(SelectedField<WilsonMatrix>& prop_msc, const SelectedField<WilsonMatrix>& prop_wm)
{
  TIMER("convert_mspincolor_from_wm(s_prop)");
  qassert(prop_wm.geo().multiplicity == 1);
  prop_msc.init(prop_wm.geo(), prop_wm.n_elems, 1);
  qacc_for(idx, prop_wm.n_elems, {
    WilsonMatrix& msc = prop_msc.get_elem(idx);
    const WilsonMatrix& wm = prop_wm.get_elem(idx);
    convert_mspincolor_from_wm(msc, wm);
  });
}

inline void convert_wm_from_mspincolor(SelectedPoints<WilsonMatrix>& prop_wm, const SelectedPoints<WilsonMatrix>& prop_msc)
{
  TIMER("convert_wm_from_mspincolor(sp_prop)");
  qassert(prop_msc.multiplicity == 1);
  prop_wm.init(prop_msc.n_points, 1);
  qacc_for(idx, prop_msc.n_points, {
    WilsonMatrix& wm = prop_wm.get_elem(idx);
    const WilsonMatrix& msc = prop_msc.get_elem(idx);
    convert_wm_from_mspincolor(wm, msc);
  });
}

inline void convert_mspincolor_from_wm(SelectedPoints<WilsonMatrix>& prop_msc, const SelectedPoints<WilsonMatrix>& prop_wm)
{
  TIMER("convert_mspincolor_from_wm(sp_prop)");
  qassert(prop_wm.multiplicity == 1);
  prop_msc.init(prop_wm.n_points, 1);
  qacc_for(idx, prop_wm.n_points, {
    WilsonMatrix& msc = prop_msc.get_elem(idx);
    const WilsonMatrix& wm = prop_wm.get_elem(idx);
    convert_mspincolor_from_wm(msc, wm);
  });
}

inline void set_t_range_flip_tpbc_with_tslice(int& t_start, int& t_stop,
                                              const int tslice_flip_tpbc,
                                              const int t_size)
// flip range t_start <= t < t_stop
// tslice_flip_tpbc - t_size_half ... tslice_flip_tpbc ... tslice_flip_tpbc +
// t_size_half - 1
{
  const int t_size_half = t_size / 2;
  qassert(t_size == t_size_half * 2);
  if (tslice_flip_tpbc + t_size_half < t_size) {
    t_start = tslice_flip_tpbc + t_size_half;
    t_stop = t_size;
  } else if (tslice_flip_tpbc - t_size_half >= 0) {
    t_start = 0;
    t_stop = tslice_flip_tpbc - t_size_half;
  } else {
    qassert(false);
  }
}

inline void flip_tpbc_with_tslice(SelectedPoints<WilsonMatrix>& ps_prop,
                                  const PointSelection& psel,
                                  const int tslice_flip_tpbc, const int t_size)
{
  if (tslice_flip_tpbc < 0) {
    return;
  }
  TIMER_VERBOSE("flip_tpbc_with_tslice(psel)");
  qassert(t_size > 0);
  qassert(t_size > tslice_flip_tpbc);
  qassert(ps_prop.n_points == (long)psel.size());
  int t_start, t_stop;
  set_t_range_flip_tpbc_with_tslice(t_start, t_stop, tslice_flip_tpbc, t_size);
  qthread_for(idx, ps_prop.n_points, {
    const Coordinate xg = psel[idx];
    const int t = xg[3];
    qassert(t_size > t);
    if (t_start <= t and t < t_stop) {
      ps_prop.get_elem(idx) *= -1;
    }
  });
}

inline void flip_tpbc_with_tslice(SelectedField<WilsonMatrix>& s_prop,
                                  const FieldSelection& fsel,
                                  const int tslice_flip_tpbc)
{
  if (tslice_flip_tpbc < 0) {
    return;
  }
  TIMER_VERBOSE("flip_tpbc_with_tslice(fsel)");
  qassert(s_prop.n_elems == (long)fsel.n_elems);
  const Geometry& geo = fsel.f_rank.geo();
  const int t_size = geo.total_site()[3];
  int t_start, t_stop;
  set_t_range_flip_tpbc_with_tslice(t_start, t_stop, tslice_flip_tpbc, t_size);
  qacc_for(idx, fsel.n_elems, {
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const int t = xg[3];
    if (t_start <= t and t < t_stop) {
      s_prop.get_elem(idx) *= -1;
    }
  });
}

// -------------------------------------------------------------------------

inline void set_tslice_mom_src_fermion_field(FermionField4d& ff,
                                             const int tslice,
                                             const CoordinateD& lmom,
                                             const int cs)
// ff need to be initialized beforehand
{
  qassert(lmom[3] == 0);
  const Geometry& geo = ff.geo();
  const CoordinateD mom = lmom * lattice_mom_mult(geo);
  set_zero(ff);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    if (xg[3] == tslice) {
      double phase = 0.0;
      for (int i = 0; i < DIMN; ++i) {
        phase += mom[i] * xg[i];
      }
      ff.get_elem(xl)(cs) = std::polar(1.0, phase);
    }
  }
}

inline void set_tslice_mom_src(Propagator4d& prop, const Geometry& geo_input,
                               const int tslice, const CoordinateD& lmom)
{
  TIMER_VERBOSE("set_tslice_mom_src");
  const Geometry geo = geo_reform(geo_input);
  prop.init(geo);
  FermionField4d src;
  src.init(geo);
  for (int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_tslice_mom_src_fermion_field(src, tslice, lmom, cs);
    set_propagator_col_from_fermion_field(prop, cs, src);
  }
}

template <class Inverter>
inline void set_tslice_mom_src_propagator(Propagator4d& prop, const int tslice,
                                          const CoordinateD& lmom,
                                          Inverter& inverter)
{
  TIMER_VERBOSE("set_tslice_mom_src_propagator");
  const Geometry& geo = geo_remult(inverter.geo());
  prop.init(geo);
  qassert(prop.geo() == geo);
  FermionField4d src, sol;
  src.init(geo);
  sol.init(geo);
  for (int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_tslice_mom_src_fermion_field(src, tslice, lmom, cs);
    set_zero(sol);
    invert(sol, src, inverter);
    set_propagator_col_from_fermion_field(prop, cs, sol);
  }
}

inline void set_volume_src_fermion_field(FermionField4d& ff,
                                         const CoordinateD& lmom, const int cs)
// ff need to be initialized beforehand
{
  qassert(lmom[3] == 0);
  const Geometry& geo = ff.geo();
  const CoordinateD mom = lmom * lattice_mom_mult(geo);
  set_zero(ff);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    double phase = 0.0;
    for (int i = 0; i < DIMN; ++i) {
      phase += mom[i] * xg[i];
    }
    ff.get_elem(xl)(cs) = std::polar(1.0, phase);
  }
}

inline void set_volume_src(Propagator4d& prop, const Geometry& geo_input,
                           const CoordinateD& lmom = CoordinateD())
{
  TIMER_VERBOSE("set_volume_src");
  const Geometry geo = geo_reform(geo_input);
  prop.init(geo);
  FermionField4d src;
  src.init(geo);
  for (int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_volume_src_fermion_field(src, lmom, cs);
    set_propagator_col_from_fermion_field(prop, cs, src);
  }
}

template <class Inverter>
inline void set_volume_src_propagator(Propagator4d& prop, const Inverter& inv,
                                      const CoordinateD& lmom = CoordinateD())
{
  TIMER_VERBOSE("set_volume_src_propagator");
  const Geometry geo = geo_reform(inv.geo());
  prop.init(geo);
  FermionField4d sol, src;
  sol.init(geo);
  src.init(geo);
  for (int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_volume_src_fermion_field(src, lmom, cs);
    set_zero(sol);
    invert(sol, src, inv);
    set_propagator_col_from_fermion_field(prop, cs, sol);
  }
}

}  // namespace qlat
