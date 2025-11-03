#include <qlat/qcd-prop.h>

namespace qlat
{  //

void set_ff_vec_from_prop(std::vector<FermionField4d>& ff_vec,
                          const Propagator4d& prop)
{
  TIMER_FLOPS("set_ff_vec_from_prop(ff_vec,prop)");
  timer.flops += get_data(prop).data_size();
  const Int num_field = 12;
  ff_vec.resize(num_field);
  const Geometry geo = prop.geo.get();
  qassert(geo.is_only_local);
  vector<FermionField4d> ffv_vec(num_field, MemType::Cpu);
  set_zero(ffv_vec);
  qfor(id_field, num_field, {
    ff_vec[id_field].init();
    ff_vec[id_field].set_mem_type(prop.field.mem_type);
    ff_vec[id_field].init(geo);
    ffv_vec[id_field].set_view(ff_vec[id_field]);
  });
  ffv_vec.set_mem_type(MemType::Acc);
  qacc_for(index, geo.local_volume(), {
    const WilsonMatrix& wm = prop.get_elem(index);
    for (Int id_field = 0; id_field < num_field; ++id_field) {
      WilsonVector& wv = ffv_vec[id_field].get_elem(index);
      set_wilson_vector_from_matrix_col(wv, wm, id_field);
    }
  });
}

void set_prop_from_ff_vec(Propagator4d& prop,
                          const std::vector<FermionField4d>& ff_vec)
{
  TIMER_FLOPS("set_prop_from_ff_vec(prop,ff_vec)");
  const Int num_field = 12;
  qassert(ff_vec.size() == num_field);
  timer.flops += get_data(ff_vec[0]).data_size() * num_field;
  const Geometry geo = ff_vec[0].geo.get();
  qassert(geo.is_only_local);
  prop.init();
  prop.set_mem_type(ff_vec[0].field.mem_type);
  prop.init(geo);
  vector<FermionField4d> ffv_vec(num_field, MemType::Cpu);
  set_zero(ffv_vec);
  qfor(id_field, num_field, {
    qassert(ff_vec[id_field].geo.get() == geo);
    ffv_vec[id_field].set_view(ff_vec[id_field]);
  });
  ffv_vec.set_mem_type(MemType::Acc);
  qacc_for(index, geo.local_volume(), {
    WilsonMatrix& wm = prop.get_elem(index);
    for (Int id_field = 0; id_field < num_field; ++id_field) {
      const WilsonVector& wv = ffv_vec[id_field].get_elem(index);
      set_wilson_matrix_col_from_vector(wm, id_field, wv);
    }
  });
}

template <class T>
static void set_point_src_fermion_field(FermionField4dT<T>& ff,
                                        const Coordinate& xg, const Int cs,
                                        const ComplexD& value = 1.0)
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
static void set_point_src(Propagator4dT<T>& prop, const Geometry& geo_input,
                          const Coordinate& xg, const ComplexD& value = 1.0)
{
  TIMER_VERBOSE("set_point_src");
  const Geometry geo = geo_resize(geo_input);
  prop.init(geo);
  FermionField4dT<T> src;
  src.init(geo);
  for (Int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_point_src_fermion_field(src, xg, cs, value);
    set_propagator_col_from_fermion_field(prop, cs, src);
  }
}

void set_point_src_fermion_field(FermionField4dT<RealD>& ff,
                                 const Coordinate& xg, const Int cs,
                                 const ComplexD& value)
{
  set_point_src_fermion_field<RealD>(ff, xg, cs, value);
}

void set_point_src_fermion_field(FermionField4dT<RealF>& ff,
                                 const Coordinate& xg, const Int cs,
                                 const ComplexD& value)
{
  set_point_src_fermion_field<RealF>(ff, xg, cs, value);
}

void set_point_src(Propagator4dT<RealD>& prop, const Geometry& geo_input,
                   const Coordinate& xg, const ComplexD& value)
{
  set_point_src<RealD>(prop, geo_input, xg, value);
}

void set_point_src(Propagator4dT<RealF>& prop, const Geometry& geo_input,
                   const Coordinate& xg, const ComplexD& value)
{
  set_point_src<RealF>(prop, geo_input, xg, value);
}

void set_wall_src_fermion_field(FermionField4d& ff, const Int tslice,
                                const CoordinateD& lmom, const Int cs)
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
      for (Int i = 0; i < DIMN; ++i) {
        phase += mom[i] * xg[i];
      }
      ff.get_elem(xl)(cs) = qpolar(1.0, phase);
    }
  }
}

void set_wall_src(Propagator4d& prop, const Geometry& geo_input,
                  const Int tslice, const CoordinateD& lmom)
{
  TIMER_VERBOSE("set_wall_src");
  const Geometry geo = geo_resize(geo_input);
  prop.init(geo);
  FermionField4d src;
  src.init(geo);
  for (Int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
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
                               const Int cs)
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
    for (Int i = 0; i < DIMN; ++i) {
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
  for (Int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_mom_src_fermion_field(src, lmom, cs);
    set_propagator_col_from_fermion_field(prop, cs, src);
  }
}

void free_invert(Prop& p_sol, const Prop& p_src, const RealD mass, const RealD m5,
                 const CoordinateD& momtwist)
{
  TIMER("free_invert");
  free_invert<RealD>(p_sol, p_src, mass, m5, momtwist);
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
                                       const Int tslice_flip_tpbc,
                                       const Int t_size)
// flip range t_start <= t < t_stop
// tslice_flip_tpbc - t_size_half ... tslice_flip_tpbc ... tslice_flip_tpbc +
// t_size_half - 1
{
  const Int t_size_half = t_size / 2;
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
                           const Int tslice_flip_tpbc, const Int t_size)
{
  if (tslice_flip_tpbc < 0) {
    return;
  }
  TIMER_VERBOSE("flip_tpbc_with_tslice(psel)");
  qassert(t_size > 0);
  qassert(t_size > tslice_flip_tpbc);
  qassert(ps_prop.n_points == (Long)psel.size());
  Int t_start, t_stop;
  set_t_range_flip_tpbc_with_tslice(t_start, t_stop, tslice_flip_tpbc, t_size);
  qthread_for(idx, ps_prop.n_points, {
    const Coordinate xg = psel[idx];
    const Int t = xg[3];
    qassert(t_size > t);
    if (t_start <= t and t < t_stop) {
      ps_prop.get_elem(idx) *= -1;
    }
  });
}

void flip_tpbc_with_tslice(SelectedField<WilsonMatrix>& s_prop,
                           const FieldSelection& fsel,
                           const Int tslice_flip_tpbc)
{
  if (tslice_flip_tpbc < 0) {
    return;
  }
  TIMER_VERBOSE("flip_tpbc_with_tslice(fsel)");
  qassert(s_prop.n_elems == (Long)fsel.n_elems);
  const Geometry& geo = fsel.f_rank.geo();
  const Int t_size = geo.total_site()[3];
  Int t_start, t_stop;
  set_t_range_flip_tpbc_with_tslice(t_start, t_stop, tslice_flip_tpbc, t_size);
  qacc_for(idx, fsel.n_elems, {
    const Long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Int t = xg[3];
    if (t_start <= t and t < t_stop) {
      s_prop.get_elem(idx) *= -1;
    }
  });
}

void set_tslice_mom_src_fermion_field(FermionField4d& ff, const Int tslice,
                                      const CoordinateD& lmom, const Int cs)
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
      for (Int i = 0; i < DIMN; ++i) {
        phase += mom[i] * xg[i];
      }
      ff.get_elem(xl)(cs) = qpolar(1.0, phase);
    }
  }
}

void set_tslice_mom_src(Propagator4d& prop, const Geometry& geo_input,
                        const Int tslice, const CoordinateD& lmom)
{
  TIMER_VERBOSE("set_tslice_mom_src");
  const Geometry geo = geo_resize(geo_input);
  prop.init(geo);
  FermionField4d src;
  src.init(geo);
  for (Int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_tslice_mom_src_fermion_field(src, tslice, lmom, cs);
    set_propagator_col_from_fermion_field(prop, cs, src);
  }
}

void set_volume_src_fermion_field(FermionField4d& ff, const CoordinateD& lmom,
                                  const Int cs)
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
    for (Int i = 0; i < DIMN; ++i) {
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
  for (Int cs = 0; cs < 4 * NUM_COLOR; ++cs) {
    set_volume_src_fermion_field(src, lmom, cs);
    set_propagator_col_from_fermion_field(prop, cs, src);
  }
}

}  // namespace qlat
