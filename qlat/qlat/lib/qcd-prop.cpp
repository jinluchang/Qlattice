#include <qlat/qcd-prop.h>

namespace qlat
{  //

void set_wall_src_fermion_field(FermionField4d& ff, const int tslice,
                                const CoordinateD& lmom, const int cs)
// ff need to be initialized beforehand
{
  qassert(lmom[3] == 0);
  const Geometry& geo = ff.geo();
  const CoordinateD mom = lmom * lattice_mom_mult(geo);
  set_zero(ff);
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    if (xg[3] == tslice) {
      double phase = 0.0;
      for (int i = 0; i < DIMN; ++i) {
        phase += mom[i] * xg[i];
      }
      ff.get_elem(xl)(cs) = qpolar(1.0, phase);
    }
  }
}

void set_wall_src(Propagator4d& prop, const Geometry& geo_input,
                  const int tslice, const CoordinateD& lmom)
{
  TIMER_VERBOSE("set_wall_src");
  const Geometry geo = geo_resize(geo_input);
  prop.init(geo);
  FermionField4d src;
  src.init(geo);
  for (int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_wall_src_fermion_field(src, tslice, lmom, cs);
    set_propagator_col_from_fermion_field(prop, cs, src);
  }
}

void set_rand_u1_src_psel(Propagator4d& prop, FieldM<ComplexD, 1>& fu1,
                          const PointsSelection& psel, const Geometry& geo_,
                          const RngState& rs)
{
  TIMER_VERBOSE("set_rand_u1_src_psel");
  const Geometry geo = geo_resize(geo_);
  const Coordinate total_site = geo.total_site();
  prop.init(geo);
  fu1.init(geo);
  set_zero(prop);
  set_zero(fu1);
  qthread_for(idx, (Long)psel.size(), {
    const Coordinate xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      const Long gindex = index_from_coordinate(xg, total_site);
      RngState rst = rs.newtype(gindex);
      const double phase = u_rand_gen(rst, PI, -PI);
      const ComplexD u1 = qpolar(1.0, phase);
      set_unit(prop.get_elem(xl), u1);
      fu1.get_elem(xl) = u1;
    }
  });
}

void set_rand_u1_sol_psel(SelectedPoints<WilsonMatrix>& sp_prop,
                          const Propagator4d& prop,
                          const FieldM<ComplexD, 1>& fu1,
                          const PointsSelection& psel)
// calculate self loop at psel locations
{
  TIMER_VERBOSE("set_rand_u1_sol_psel");
  SelectedPoints<ComplexD> sp_fu1;
  set_selected_points(sp_prop, prop, psel);
  set_selected_points(sp_fu1, fu1, psel);
  qthread_for(idx, (Long)psel.size(), {
    const ComplexD& u1 = sp_fu1.get_elem(idx);
    WilsonMatrix& wm = sp_prop.get_elem(idx);
    wm *= qlat::qconj(u1);
  });
}

void set_rand_u1_src_fsel(Propagator4d& prop, FieldM<ComplexD, 1>& fu1,
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
    const Long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    qassert(geo.is_local(xl));
    const Long gindex = index_from_coordinate(xg, total_site);
    RngState rst = rs.newtype(gindex);
    const double phase = u_rand_gen(rst, PI, -PI);
    const ComplexD u1 = qpolar(1.0, phase);
    set_unit(prop.get_elem(xl), u1);
    fu1.get_elem(xl) = u1;
  });
}

void set_rand_u1_sol_fsel(SelectedField<WilsonMatrix>& sf_prop,
                          const Propagator4d& prop,
                          const FieldM<ComplexD, 1>& fu1,
                          const FieldSelection& fsel)
// calculate self loop at fsel locations
{
  TIMER_VERBOSE("set_rand_u1_sol_fsel");
  SelectedField<ComplexD> sf_fu1;
  set_selected_field(sf_prop, prop, fsel);
  set_selected_field(sf_fu1, fu1, fsel);
  qthread_for(idx, fsel.n_elems, {
    const ComplexD& u1 = sf_fu1.get_elem(idx);
    WilsonMatrix& wm = sf_prop.get_elem(idx);
    wm *= qlat::qconj(u1);
  });
}

void set_rand_vol_u1(Field<ComplexD>& fu1, const Geometry& geo_input,
                     const RngState& rs)
// fu1.multiplicity == 1
{
  TIMER_VERBOSE("set_rand_vol_u1");
  const Geometry geo = geo_resize(geo_input);
  const Coordinate total_site = geo.total_site();
  fu1.init(geo, 1);
  set_zero(fu1);
  qthread_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    qassert(geo.is_local(xl));
    const Long gindex = index_from_coordinate(xg, total_site);
    RngState rst = rs.newtype(gindex);
    const RealD phase = u_rand_gen(rst, PI, -PI);
    const ComplexD u1 = qpolar(1.0, phase);
    fu1.get_elem(xl) = u1;
  });
}

void set_rand_vol_u1_src(Propagator4d& prop, const Field<ComplexD>& fu1)
// fu1.multiplicity == 1
{
  TIMER_VERBOSE("set_rand_vol_u1_src");
  const Geometry geo = geo_resize(fu1.geo());
  prop.init(geo);
  set_zero(prop);
  qthread_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const ComplexD u1 = fu1.get_elem(xl);
    set_unit(prop.get_elem(xl), u1);
  });
}

void set_mom_src_fermion_field(FermionField4d& ff, const CoordinateD& lmom,
                               const int cs)
// ff need to be initialized beforehand
{
  const Geometry& geo = ff.geo();
  const CoordinateD mom = lmom * lattice_mom_mult(geo);
  set_zero(ff);
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    double phase = 0.0;
    for (int i = 0; i < DIMN; ++i) {
      phase += mom[i] * xg[i];
    }
    ff.get_elem(xl)(cs) = qpolar(1.0, phase);
  }
}

void set_mom_src(Propagator4d& prop, const Geometry& geo_input,
                 const CoordinateD& lmom)
{
  TIMER_VERBOSE("set_mom_src");
  const Geometry geo = geo_resize(geo_input);
  prop.init(geo);
  qassert(prop.geo() == geo);
  FermionField4d src;
  src.init(geo);
  for (int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_mom_src_fermion_field(src, lmom, cs);
    set_propagator_col_from_fermion_field(prop, cs, src);
  }
}

void convert_wm_from_mspincolor(Propagator4d& prop_wm,
                                const Propagator4d& prop_msc)
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

void convert_mspincolor_from_wm(Propagator4d& prop_msc,
                                const Propagator4d& prop_wm)
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

void convert_wm_from_mspincolor(SelectedField<WilsonMatrix>& prop_wm,
                                const SelectedField<WilsonMatrix>& prop_msc)
{
  TIMER("convert_wm_from_mspincolor(s_prop)");
  qassert(prop_msc.multiplicity == 1);
  prop_wm.init(prop_msc.geo(), prop_msc.n_elems, 1);
  qacc_for(idx, prop_msc.n_elems, {
    WilsonMatrix& wm = prop_wm.get_elem(idx);
    const WilsonMatrix& msc = prop_msc.get_elem(idx);
    convert_wm_from_mspincolor(wm, msc);
  });
}

void convert_mspincolor_from_wm(SelectedField<WilsonMatrix>& prop_msc,
                                const SelectedField<WilsonMatrix>& prop_wm)
{
  TIMER("convert_mspincolor_from_wm(s_prop)");
  qassert(prop_wm.multiplicity == 1);
  prop_msc.init(prop_wm.geo(), prop_wm.n_elems, 1);
  qacc_for(idx, prop_wm.n_elems, {
    WilsonMatrix& msc = prop_msc.get_elem(idx);
    const WilsonMatrix& wm = prop_wm.get_elem(idx);
    convert_mspincolor_from_wm(msc, wm);
  });
}

void convert_wm_from_mspincolor(SelectedPoints<WilsonMatrix>& prop_wm,
                                const SelectedPoints<WilsonMatrix>& prop_msc)
{
  TIMER("convert_wm_from_mspincolor(sp_prop)");
  qassert(prop_msc.multiplicity == 1);
  prop_wm.init(prop_msc.n_points, 1, prop_msc.points_dist_type);
  qacc_for(idx, prop_msc.n_points, {
    WilsonMatrix& wm = prop_wm.get_elem(idx);
    const WilsonMatrix& msc = prop_msc.get_elem(idx);
    convert_wm_from_mspincolor(wm, msc);
  });
}

void convert_mspincolor_from_wm(SelectedPoints<WilsonMatrix>& prop_msc,
                                const SelectedPoints<WilsonMatrix>& prop_wm)
{
  TIMER("convert_mspincolor_from_wm(sp_prop)");
  qassert(prop_wm.multiplicity == 1);
  prop_msc.init(prop_wm.n_points, 1, prop_wm.points_dist_type);
  qacc_for(idx, prop_wm.n_points, {
    WilsonMatrix& msc = prop_msc.get_elem(idx);
    const WilsonMatrix& wm = prop_wm.get_elem(idx);
    convert_mspincolor_from_wm(msc, wm);
  });
}

void set_t_range_flip_tpbc_with_tslice(int& t_start, int& t_stop,
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

void flip_tpbc_with_tslice(SelectedPoints<WilsonMatrix>& ps_prop,
                           const PointsSelection& psel,
                           const int tslice_flip_tpbc, const int t_size)
{
  if (tslice_flip_tpbc < 0) {
    return;
  }
  TIMER_VERBOSE("flip_tpbc_with_tslice(psel)");
  qassert(t_size > 0);
  qassert(t_size > tslice_flip_tpbc);
  qassert(ps_prop.n_points == (Long)psel.size());
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

void flip_tpbc_with_tslice(SelectedField<WilsonMatrix>& s_prop,
                           const FieldSelection& fsel,
                           const int tslice_flip_tpbc)
{
  if (tslice_flip_tpbc < 0) {
    return;
  }
  TIMER_VERBOSE("flip_tpbc_with_tslice(fsel)");
  qassert(s_prop.n_elems == (Long)fsel.n_elems);
  const Geometry& geo = fsel.f_rank.geo();
  const int t_size = geo.total_site()[3];
  int t_start, t_stop;
  set_t_range_flip_tpbc_with_tslice(t_start, t_stop, tslice_flip_tpbc, t_size);
  qacc_for(idx, fsel.n_elems, {
    const Long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const int t = xg[3];
    if (t_start <= t and t < t_stop) {
      s_prop.get_elem(idx) *= -1;
    }
  });
}

void set_tslice_mom_src_fermion_field(FermionField4d& ff, const int tslice,
                                      const CoordinateD& lmom, const int cs)
// ff need to be initialized beforehand
{
  qassert(lmom[3] == 0);
  const Geometry& geo = ff.geo();
  const CoordinateD mom = lmom * lattice_mom_mult(geo);
  set_zero(ff);
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    if (xg[3] == tslice) {
      double phase = 0.0;
      for (int i = 0; i < DIMN; ++i) {
        phase += mom[i] * xg[i];
      }
      ff.get_elem(xl)(cs) = qpolar(1.0, phase);
    }
  }
}

void set_tslice_mom_src(Propagator4d& prop, const Geometry& geo_input,
                        const int tslice, const CoordinateD& lmom)
{
  TIMER_VERBOSE("set_tslice_mom_src");
  const Geometry geo = geo_resize(geo_input);
  prop.init(geo);
  FermionField4d src;
  src.init(geo);
  for (int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_tslice_mom_src_fermion_field(src, tslice, lmom, cs);
    set_propagator_col_from_fermion_field(prop, cs, src);
  }
}

void set_volume_src_fermion_field(FermionField4d& ff, const CoordinateD& lmom,
                                  const int cs)
// ff need to be initialized beforehand
{
  qassert(lmom[3] == 0);
  const Geometry& geo = ff.geo();
  const CoordinateD mom = lmom * lattice_mom_mult(geo);
  set_zero(ff);
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    double phase = 0.0;
    for (int i = 0; i < DIMN; ++i) {
      phase += mom[i] * xg[i];
    }
    ff.get_elem(xl)(cs) = qpolar(1.0, phase);
  }
}

void set_volume_src(Propagator4d& prop, const Geometry& geo_input,
                    const CoordinateD& lmom)
{
  TIMER_VERBOSE("set_volume_src");
  const Geometry geo = geo_resize(geo_input);
  prop.init(geo);
  FermionField4d src;
  src.init(geo);
  for (int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_volume_src_fermion_field(src, lmom, cs);
    set_propagator_col_from_fermion_field(prop, cs, src);
  }
}

}  // namespace qlat
