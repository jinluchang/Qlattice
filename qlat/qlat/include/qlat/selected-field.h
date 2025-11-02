#pragma once

#include <qlat/field.h>
#include <qlat/selected-points.h>

namespace qlat
{  //

void mk_field_selection(FieldRank& f_rank, const Geometry& geo,
                        const int64_t val = 0);

void mk_field_selection(FieldRank& f_rank, const Coordinate& total_site,
                        const int64_t val = 0);

void mk_field_selection(FieldRank& f_rank, const Coordinate& total_site,
                        const PointsSelection& psel,
                        const Long rank_xgs = 1024L * 1024L * 1024L * 1024L *
                                              1024L);

void set_selected_gindex(SelectedField<Long>& sfgi, const FieldSelection& fsel);

void mk_grid_field_selection(FieldRank& f_rank, const Coordinate& total_site,
                             const Long n_per_tslice_, const RngState& rs);

void mk_field_selection(FieldRank& f_rank, const Coordinate& total_site,
                        const Long n_per_tslice, const RngState& rs);

void add_field_selection(FieldRank& f_rank, const PointsSelection& psel,
                         const Long rank_psel = 1024L * 1024L * 1024L * 1024L *
                                                1024L);

void add_field_selection(FieldRank& f_rank, const FieldSelection& fsel);

void select_rank_range(FieldRank& f_rank, const Long rank_start = 0,
                       const Long rank_stop = -1);

void select_t_range(FieldRank& f_rank, const Long t_start = 0,
                    const Long t_stop = -1);

void set_n_per_tslice(FieldRank& f_rank, const Long n_per_tslice);

void update_field_selection(FieldSelection& fsel);

void set_grid_field_selection(FieldSelection& fsel,
                              const Coordinate& total_site,
                              const Long n_per_tslice, const RngState& rs);

void set_field_selection(FieldSelection& fsel, const FieldRank& f_rank);

void set_field_selection(FieldSelection& fsel, const Coordinate& total_site);

void set_field_selection(FieldSelection& fsel, const Coordinate& total_site,
                         const Long n_per_tslice, const RngState& rs);

void set_field_selection(FieldSelection& fsel, const Coordinate& total_site,
                         const Long n_per_tslice, const RngState& rs,
                         const PointsSelection& psel);

bool is_matching_fsel(const FieldSelection& fsel1, const FieldSelection& fsel2);

bool is_containing(const PointsSelection& psel, const FieldSelection& fsel_small);

bool is_containing(const PointsSelection& psel, const PointsSelection& psel_small);

bool is_containing(const FieldSelection& fsel, const FieldSelection& fsel_small);

bool is_containing(const FieldSelection& fsel, const PointsSelection& psel_small);

void intersect_with(FieldSelection& fsel, const FieldSelection& fsel1);

PointsSelection intersect(const FieldSelection& fsel,
                          const PointsSelection& psel);

PointsSelection psel_from_fsel(const FieldSelection& fsel);

PointsSelection psel_from_fsel_local(const FieldSelection& fsel);

Long write_field_selection(const FieldSelection& fsel, const std::string& path);

Long read_field_selection(FieldSelection& fsel, const std::string& path);

Long idx_from_xg(const Coordinate& xg, const FieldSelection& fsel);

// ---------------------------------------

template <class M>
bool is_initialized(const SelectedField<M>& sf)
{
  return sf.initialized;
}

template <class M>
bool is_consistent(const SelectedField<M>& sf, const FieldSelection& fsel)
{
  return sf.initialized and sf.n_elems == fsel.n_elems and
         sf.geo() == fsel.f_local_idx.geo() and
         fsel.f_rank.geo() == fsel.f_local_idx.geo() and
         (Long) sf.field.size() == sf.n_elems * (Long)sf.multiplicity and
         fsel.f_local_idx.geo().is_only_local;
}

template <class M>
SelectedField<M>& operator+=(SelectedField<M>& f, const SelectedField<M>& f1)
{
  TIMER("sel_field_operator+=");
  if (not f.initialized) {
    f = f1;
  } else {
    qassert(f1.initialized);
    qassert(is_matching_geo(f.geo(), f1.geo()));
    qassert(f.multiplicity == f1.multiplicity);
    qassert(f.field.size() == f1.field.size());
#pragma omp parallel for
    for (Long k = 0; k < (Long)f.field.size(); ++k) {
      f.field[k] += f1.field[k];
    }
  }
  return f;
}

template <class M>
SelectedField<M>& operator-=(SelectedField<M>& f, const SelectedField<M>& f1)
{
  TIMER("sel_field_operator-=");
  if (not f.initialized) {
    f.init(f1.geo(), f1.n_elems, f1.multiplicity);
    set_zero(f);
    f -= f1;
  } else {
    qassert(f1.initialized);
    qassert(is_matching_geo(f.geo(), f1.geo()));
    qassert(f.multiplicity == f1.multiplicity);
    qassert(f.field.size() == f1.field.size());
#pragma omp parallel for
    for (Long k = 0; k < (Long)f.field.size(); ++k) {
      f.field[k] -= f1.field[k];
    }
  }
  return f;
}

template <class M>
SelectedField<M>& operator*=(SelectedField<M>& f, const double factor)
{
  TIMER("sel_field_operator*=(F,D)");
  qassert(f.initialized);
#pragma omp parallel for
  for (Long k = 0; k < (Long)f.field.size(); ++k) {
    f.field[k] *= factor;
  }
  return f;
}

template <class M>
SelectedField<M>& operator*=(SelectedField<M>& f, const ComplexD factor)
{
  TIMER("sel_field_operator*=(F,C)");
  qassert(f.initialized);
#pragma omp parallel for
  for (Long k = 0; k < (Long)f.field.size(); ++k) {
    f.field[k] *= factor;
  }
  return f;
}

template <class M>
void only_keep_selected_points(Field<M>& f, const FieldSelection& fsel)
{
  TIMER("only_keep_selected_points");
  qassert(f.geo().is_only_local);
  qassert(fsel.f_local_idx.geo().is_only_local);
  qassert(f.geo() == fsel.f_local_idx.geo());
  const Geometry& geo = f.geo();
  const FieldM<Long, 1>& f_local_idx = fsel.f_local_idx;
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Long idx = f_local_idx.get_elems_const(index)[0];
    if (idx < 0) {
      qassert(idx == -1);
      Vector<M> fv = f.get_elems(index);
      set_zero(fv);
    }
  }
}

template <class M>
RealD qnorm(const SelectedField<M>& sf)
{
  TIMER("qnorm(sf)");
  RealD s = qnorm(sf.field);
  glb_sum(s);
  return s;
}

template <class M>
void qnorm_field(SelectedField<RealD>& f, const SelectedField<M>& f1)
{
  TIMER("qnorm_field(f,f1)");
  const Geometry& geo = f1.geo();
  f.init(geo, f1.n_elems, 1);
  qacc_for(idx, f.n_elems, {
    const Vector<M> f1v = f1.get_elems_const(idx);
    f.get_elem(idx) = qnorm(f1v);
  });
}

void set_sqrt_field(SelectedField<RealD>& f, const SelectedField<RealD>& f1);

// -------------------------------------------

template <class M>
void set_selected_field(SelectedField<M>& sf, const Field<M>& f,
                        const FieldSelection& fsel)
{
  TIMER("set_selected_field(sf,f,fsel)");
  qassert(f.geo().is_only_local);
  qassert(fsel.f_local_idx.geo().is_only_local);
  qassert(f.geo() == fsel.f_local_idx.geo());
  // const Geometry& geo = f.geo();
  const Int multiplicity = f.multiplicity;
  sf.init(fsel, multiplicity);
  qacc_for(idx, fsel.n_elems, {
    const Long index = fsel.indices[idx];
    const Vector<M> fv = f.get_elems_const(index);
    Vector<M> sfv = sf.get_elems(idx);
    for (Int m = 0; m < multiplicity; ++m) {
      sfv[m] = fv[m];
    }
  });
}

template <class M>
void set_selected_field(SelectedField<M>& sf, const SelectedField<M>& sf0,
                        const FieldSelection& fsel, const FieldSelection& fsel0,
                        const bool is_keeping_data = true)
{
  if (&sf == &sf0) {
    return;
  }
  TIMER("set_selected_field(sf,sf0,fsel,fsel0)");
  if (&fsel == &fsel0) {
    sf = sf0;
    return;
  }
  qassert(not sf0.geo.null());
  qassert(not fsel.f_local_idx.geo.null());
  qassert(not fsel0.f_local_idx.geo.null());
  qassert(sf0.geo().is_only_local);
  qassert(fsel.f_local_idx.geo().is_only_local);
  qassert(fsel0.f_local_idx.geo().is_only_local);
  qassert(sf0.geo() == fsel0.f_local_idx.geo());
  qassert(sf0.geo() == fsel.f_local_idx.geo());
  // const Geometry& geo = sf0.geo();
  const Int multiplicity = sf0.multiplicity;
  if (is_keeping_data) {
    sf.init_zero(fsel, multiplicity);
  } else {
    sf.init(fsel, multiplicity);
    set_zero(sf);
  }
  qacc_for(idx, fsel.n_elems, {
    const Long index = fsel.indices[idx];
    const Long idx0 = fsel0.f_local_idx.get_elem(index);
    if (idx0 >= 0) {
      Vector<M> sfv = sf.get_elems(idx);
      const Vector<M> fv = sf0.get_elems_const(idx0);
      for (Int m = 0; m < multiplicity; ++m) {
        sfv[m] = fv[m];
      }
    }
  });
}

template <class M>
void set_selected_field(SelectedField<M>& sf, const SelectedPoints<M>& sp,
                        const FieldSelection& fsel, const PointsSelection& psel,
                        const bool is_keeping_data = true)
{
  TIMER("set_selected_field(sf,sp,fsel,psel)");
  qassert(fsel.f_local_idx.geo().is_only_local);
  qassert(sf.geo() == fsel.f_local_idx.geo());
  const Long n_points = sp.n_points;
  qassert(n_points == (Long)psel.size());
  const Geometry& geo = fsel.f_rank.geo();
  const Int multiplicity = sp.multiplicity;
  if (is_keeping_data) {
    sf.init_zero(fsel, multiplicity);
  } else {
    sf.init(fsel, multiplicity);
    set_zero(sf);
  }
  qthread_for(idx, n_points, {
    const Coordinate& xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      const Long sf_idx = fsel.f_local_idx.get_elem(xl);
      if (sf_idx >= 0) {
        qassert(sf_idx < sf.n_elems);
        const Vector<M> spv = sp.get_elems_const(idx);
        Vector<M> fv = sf.get_elems(sf_idx);
        for (Int m = 0; m < multiplicity; ++m) {
          fv[m] = spv[m];
        }
      }
    }
  });
}

template <class M>
void set_selected_points(SelectedPoints<M>& sp, const SelectedField<M>& sf,
                         const PointsSelection& psel,
                         const FieldSelection& fsel,
                         const bool is_keeping_data = true)
// only assign available points
{
  TIMER("set_selected_points(sp,sf,psel,fsel)");
  const Geometry& geo = sf.geo();
  qassert(is_consistent(sf, fsel));
  qassert(psel.points_dist_type == PointsDistType::Global);
  const Long n_points = psel.size();
  SelectedPoints<M> sp_tmp;
  sp_tmp.init(psel, sf.multiplicity);
  set_zero(sp_tmp);
  SelectedPoints<int8_t> sp_count;
  sp_count.init(psel, 1);
  set_zero(sp_count);
  qthread_for(idx, n_points, {
    const Coordinate& xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      const Long sf_idx = fsel.f_local_idx.get_elem(xl);
      if (sf_idx >= 0) {
        qassert(sf_idx < sf.n_elems);
        sp_count.get_elem(idx) += 1;
        Vector<M> spv = sp_tmp.get_elems(idx);
        const Vector<M> fv = sf.get_elems_const(sf_idx);
        for (Int m = 0; m < sf.multiplicity; ++m) {
          spv[m] = fv[m];
        }
      }
    }
  });
  glb_sum(get_data_char(sp_tmp.points));
  if (is_keeping_data) {
    glb_sum(get_data_char(sp_count.points));
    sp.init_zero(psel, sf.multiplicity);
    qthread_for(idx, n_points, {
      if (sp_count.get_elem(idx) > 0) {
        Vector<M> spv = sp.get_elems(idx);
        const Vector<M> spv_tmp = sp_tmp.get_elems_const(idx);
        for (Int m = 0; m < sf.multiplicity; ++m) {
          spv[m] = spv_tmp[m];
        }
      }
    });
  } else {
    sp = sp_tmp;
  }
}

template <class M>
void set_field_selected(Field<M>& f, const SelectedField<M>& sf,
                        const FieldSelection& fsel,
                        const bool is_keeping_data = true)
{
  TIMER("set_field_selected(f,sf,fsel)");
  qassert(sf.geo().is_only_local);
  qassert(fsel.f_local_idx.geo().is_only_local);
  qassert(sf.geo() == fsel.f_local_idx.geo());
  const Geometry& geo = sf.geo();
  const Int multiplicity = sf.multiplicity;
  if (is_keeping_data) {
    f.init_zero(geo, multiplicity);
  } else {
    f.init(geo, multiplicity);
    set_zero(f);
  }
  qacc_for(idx, fsel.n_elems, {
    const Long index = fsel.indices[idx];
    Vector<M> fv = f.get_elems(index);
    const Vector<M> sfv = sf.get_elems_const(idx);
    for (Int m = 0; m < multiplicity; ++m) {
      fv[m] = sfv[m];
    }
  });
}

// -------------------------------------------

template <class M>
bool is_consistent(const SelectedPoints<M>& sp, const SelectedField<M>& sf,
                   const PointsSelection& psel, const FieldSelection& fsel)
{
  TIMER("is_consistent(sp,sf)");
  qassert(is_consistent(sp, psel));
  qassert(is_consistent(sf, fsel));
  const Geometry& geo = sf.geo();
  const Long n_points = psel.size();
  double qnorm_diff = 0.0;
  qfor(idx, n_points, {
    const Coordinate& xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      const Long sf_idx = fsel.f_local_idx.get_elem(xl);
      if (sf_idx >= 0) {
        const Vector<M> fv = sf.get_elems_const(sf_idx);
        const Vector<M> spv = sp.get_elems_const(idx);
        for (Int m = 0; m < sf.multiplicity; ++m) {
          qnorm_diff += qnorm(spv[m] - fv[m]);
        }
      }
    }
  });
  glb_sum(qnorm_diff);
  return qnorm_diff == 0.0;
}

template <class M>
void acc_field(Field<M>& f, const SelectedField<M>& sf,
               const FieldSelection& fsel)
// f can be empty
{
  TIMER("acc_field(f,sf,fsel)");
  const Geometry& geo = fsel.f_rank.geo();
  const Int multiplicity = sf.multiplicity;
  if (not is_initialized(f)) {
    f.init(geo, multiplicity);
    set_zero(f);
  }
  qassert(multiplicity == f.multiplicity);
  qassert(sf.n_elems == fsel.n_elems);
  qacc_for(idx, fsel.n_elems, {
    const Long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> fv = f.get_elems(xl);
    const Vector<M> sfv = sf.get_elems_const(idx);
    for (Int m = 0; m < multiplicity; ++m) {
      fv[m] += sfv[m];
    }
  });
}

template <class M>
std::vector<M> field_sum_tslice(const SelectedField<M>& sf,
                                const FieldSelection& fsel, const Int t_dir = 3)
// length = t_size * multiplicity
{
  TIMER("field_sum_tslice");
  qassert(sf.geo().is_only_local);
  qassert(fsel.f_local_idx.geo().is_only_local);
  qassert(sf.geo() == fsel.f_local_idx.geo());
  const Geometry& geo = sf.geo();
  const Int t_size = geo.total_site()[t_dir];
  const Int multiplicity = sf.multiplicity;
  std::vector<M> vec(t_size * multiplicity);
  set_zero(vec);
  for (Long idx = 0; idx < fsel.n_elems; ++idx) {
    const Long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Vector<M> sfv = sf.get_elems_const(idx);
    for (Int m = 0; m < multiplicity; ++m) {
      vec[xg[t_dir] * multiplicity + m] += sfv[m];
    }
  }
  return vec;
}

template <class M>
void field_glb_sum_tslice(SelectedPoints<M>& sp, const SelectedField<M>& sf,
                          const FieldSelection& fsel, const Int t_dir = 3)
{
  TIMER("field_glb_sum_tslice(sp,sf,fsel)");
  sp.init();
  const Geometry& geo = sf.geo();
  const Int t_size = geo.total_site()[t_dir];
  const Int multiplicity = sf.multiplicity;
  std::vector<M> vec = field_sum_tslice(sf, fsel, t_dir);
  glb_sum(vec);
  sp.init(t_size, multiplicity, PointsDistType::Global);
  sp.points = vec;
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
void set_u_rand(SelectedField<M>& sf, const FieldSelection& fsel,
                       const RngState& rs, const RealD upper = 1.0,
                       const RealD lower = -1.0)
{
  TIMER("set_u_rand(sf,fsel,rs,upper,lower)");
  if (not is_composed_of_real<M>()) {
    qassert(is_composed_of_real<M>());
    return;
  }
  using Real = typename IsDataValueType<M>::ElementaryType;
  const Geometry& geo = sf.geo();
  qassert(geo.is_only_local);
  qassert(fsel.f_local_idx.geo().is_only_local);
  qassert(geo == fsel.f_local_idx.geo());
  qthread_for(idx, fsel.n_elems, {
    const Long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Long gindex = geo.g_index_from_g_coordinate(xg);
    RngState rsi = rs.newtype(gindex);
    Vector<M> v = sf.get_elems(idx);
    Vector<Real> dv((Real*)v.data(), v.data_size() / sizeof(Real));
    for (Int m = 0; m < dv.size(); ++m) {
      dv[m] = u_rand_gen(rsi, upper, lower);
    }
  });
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
void set_g_rand(SelectedField<M>& sf, const FieldSelection& fsel,
                       const RngState& rs, const RealD center = 0.0,
                       const RealD sigma = 1.0)
{
  TIMER("set_g_rand(sf,fsel,rs,center,sigma)");
  if (not is_composed_of_real<M>()) {
    qassert(is_composed_of_real<M>());
    return;
  }
  using Real = typename IsDataValueType<M>::ElementaryType;
  const Geometry& geo = sf.geo();
  qassert(geo.is_only_local);
  qassert(fsel.f_local_idx.geo().is_only_local);
  qassert(geo == fsel.f_local_idx.geo());
  qthread_for(idx, fsel.n_elems, {
    const Long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Long gindex = geo.g_index_from_g_coordinate(xg);
    RngState rsi = rs.newtype(gindex);
    Vector<M> v = sf.get_elems(idx);
    Vector<Real> dv((Real*)v.data(), v.data_size() / sizeof(Real));
    for (Int m = 0; m < dv.size(); ++m) {
      dv[m] = g_rand_gen(rsi, center, sigma);
    }
  });
}

template <class M, class N>
void convert_field_float_from_double(SelectedField<N>& ff,
                                     const SelectedField<M>& f)
// interface_function
{
  TIMER("convert_field_float_from_double(sf)");
  qassert(f.geo().is_only_local);
  qassert(sizeof(M) % sizeof(RealD) == 0);
  qassert(sizeof(N) % sizeof(RealF) == 0);
  qassert(f.multiplicity * sizeof(M) / 2 % sizeof(N) == 0);
  const Int multiplicity = f.multiplicity * sizeof(M) / 2 / sizeof(N);
  const Geometry& geo = f.geo();
  const Long n_elems = f.n_elems;
  ff.init(geo, n_elems, multiplicity);
  const Vector<M> fdata = get_data(f);
  const Vector<double> fd((double*)fdata.data(),
                          fdata.data_size() / sizeof(RealD));
  Vector<N> ffdata = get_data(ff);
  Vector<RealF> ffd((RealF*)ffdata.data(), ffdata.data_size() / sizeof(RealF));
  qassert(ffd.size() == fd.size());
  qacc_for(i, ffd.size(), { ffd[i] = fd[i]; });
}

template <class M, class N>
void convert_field_double_from_float(SelectedField<N>& ff,
                                     const SelectedField<M>& f)
// interface_function
{
  TIMER("convert_field_double_from_float(sf)");
  qassert(f.geo().is_only_local);
  qassert(sizeof(M) % sizeof(RealF) == 0);
  qassert(sizeof(N) % sizeof(RealD) == 0);
  qassert(f.multiplicity * sizeof(M) * 2 % sizeof(N) == 0);
  const Int multiplicity = f.multiplicity * sizeof(M) * 2 / sizeof(N);
  const Geometry& geo = f.geo();
  const Long n_elems = f.n_elems;
  ff.init(geo, n_elems, multiplicity);
  const Vector<M> fdata = get_data(f);
  const Vector<RealF> fd((RealF*)fdata.data(),
                         fdata.data_size() / sizeof(RealF));
  Vector<N> ffdata = get_data(ff);
  Vector<double> ffd((double*)ffdata.data(),
                     ffdata.data_size() / sizeof(RealD));
  qassert(ffd.size() == fd.size());
  qacc_for(i, ffd.size(), { ffd[i] = fd[i]; });
}

}  // namespace qlat
