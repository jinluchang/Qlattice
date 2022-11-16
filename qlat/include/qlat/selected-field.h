#pragma once

#include <qlat/field.h>
#include <qlat/selected-points.h>

namespace qlat
{  //

inline void add_field_selection(FieldM<int64_t, 1>& f_rank,
                                const PointSelection& psel,
                                const long rank_psel = 1024L * 1024L * 1024L *
                                                       1024L * 1024L)
// interface function
// add psel points to f_rank. (only lower rank if already selected)
{
  TIMER_VERBOSE("add_field_selection(psel)");
  const Geometry& geo = f_rank.geo();
#pragma omp parallel for
  for (long i = 0; i < (long)psel.size(); ++i) {
    const Coordinate xl = geo.coordinate_l_from_g(psel[i]);
    if (geo.is_local(xl)) {
      int64_t& rank = f_rank.get_elem(xl);
      if (rank < 0 or rank > rank_psel) {
        rank = rank_psel;
      }
    }
  }
}

inline void mk_field_selection(FieldM<int64_t, 1>& f_rank,
                               const Coordinate& total_site,
                               const int64_t val = 0)
// interface function
// select everything with val
// default val = 0 ; means selection everything
// val = -1 deselection everything
{
  TIMER_VERBOSE("mk_field_selection(f_rank,total_site)");
  Geometry geo;
  geo.init(total_site, 1);
  f_rank.init();
  f_rank.init(geo);
  qassert(f_rank.geo().is_only_local);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    int64_t& rank = f_rank.get_elem(index);
    rank = val;
  }
}

inline void mk_field_selection(FieldM<int64_t, 1>& f_rank,
                               const Coordinate& total_site,
                               const std::vector<Coordinate>& xgs,
                               const long rank_xgs = 1024L * 1024L * 1024L *
                                                     1024L * 1024L)
{
  TIMER_VERBOSE("mk_field_selection(xgs)");
  mk_field_selection(f_rank, total_site, -1);
  add_field_selection(f_rank, xgs, rank_xgs);
}

inline void select_rank_range(FieldM<int64_t, 1>& f_rank,
                              const long rank_start = 0,
                              const long rank_stop = -1)
// keep rank info if rank_start <= rank and (rank < rank_stop or rank_stop == -1)
// otherwise rank = -1
// default parameter does not change selection
// but will erase the rank information for points not selected (rank = -1)
{
  TIMER_VERBOSE("select_rank_range");
  const Geometry& geo = f_rank.geo();
  qassert(geo.is_only_local);
  //const Coordinate total_site = geo.total_site();
  qacc_for(index, geo.local_volume(), {
    int64_t& rank = f_rank.get_elem(index);
    if (not(rank_start <= rank and (rank < rank_stop or rank_stop == -1))) {
      rank = -1;
    }
  });
}

inline void select_t_range(FieldM<int64_t, 1>& f_rank, const long t_start = 0,
                           const long t_stop = -1)
// keep rank info if t_start <= t and (t < t_stop or t_stop == -1)
// otherwise rank = -1
// default parameter does not change selection
// but will erase the rank information for points not selected (rank = -1)
{
  TIMER_VERBOSE("select_t_range");
  const Geometry& geo = f_rank.geo();
  qassert(geo.is_only_local);
  //const Coordinate total_site = geo.total_site();
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const int t = xg[3];
    if (not(t_start <= t and (t < t_stop or t_stop == -1))) {
      int64_t& rank = f_rank.get_elem(index);
      rank = -1;
    }
  });
}

inline void set_n_per_tslice(FieldM<int64_t, 1>& f_rank,
                             const long n_per_tslice)
// will erase the rank information for points not selected (rank = -1)
//
// if n_per_tslice == -1 then all points are selected regardless of rank
// (if a point was not selected before (rank < 0), then rank will be set to be
// rank = spatial_vol).
//
// otherwise: n_per_tslice is not enforced but only serve as an limit for f_rank
//
// 0 <= rank < n_per_tslice
//
// if point is not selected, rank = -1
{
  TIMER_VERBOSE("set_n_per_tslice");
  const Geometry& geo = f_rank.geo();
  qassert(geo.is_only_local);
  const Coordinate total_site = geo.total_site();
  const long spatial_vol = total_site[0] * total_site[1] * total_site[2];
  qassert(n_per_tslice == -1 or
          (0 <= n_per_tslice and n_per_tslice <= spatial_vol));
  qacc_for(index, geo.local_volume(), {
    int64_t& rank = f_rank.get_elem(index);
    if (n_per_tslice == -1 and rank < 0) {
      rank = spatial_vol;
    } else if (not(0 <= rank and rank < n_per_tslice)) {
      rank = -1;
    }
  });
}

inline void update_field_selection(FieldSelection& fsel)
// interface function
// update fsel based only on f_rank
// do not touch n_per_tslice and prob at all
{
  TIMER_VERBOSE("update_field_selection");
  const Geometry& geo = fsel.f_rank.geo();
  qassert(geo.is_only_local);
  fsel.f_local_idx.init();
  fsel.f_local_idx.init(geo);
  long n_elems = 0;
  for (long index = 0; index < geo.local_volume(); ++index) {
    const int64_t& rank = fsel.f_rank.get_elem(index);
    long& idx = fsel.f_local_idx.get_elem(index);
    if (0 <= rank) {
      idx = n_elems;
      n_elems += 1;
    } else {
      idx = -1;
    }
  }
  fsel.n_elems = n_elems;
  fsel.ranks.resize(fsel.n_elems);
  fsel.indices.resize(fsel.n_elems);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const long idx = fsel.f_local_idx.get_elems_const(index)[0];
    if (idx >= 0) {
      const long rank = fsel.f_rank.get_elem(index);
      fsel.ranks[idx] = rank;
      fsel.indices[idx] = index;
    }
  }
}

inline void update_field_selection(FieldSelection& fsel,
                                   const long n_per_tslice_)
// interface function
// only adjust parameter, do not change contents
{
  const Geometry& geo = fsel.f_rank.geo();
  qassert(geo.is_only_local);
  const Coordinate total_site = geo.total_site();
  const long spatial_vol = total_site[0] * total_site[1] * total_site[2];
  qassert(n_per_tslice_ == -1 or
          (0 <= n_per_tslice_ and n_per_tslice_ <= spatial_vol));
  fsel.n_per_tslice = n_per_tslice_ == -1 ? spatial_vol : n_per_tslice_;
  fsel.prob = (double)fsel.n_per_tslice / (double)spatial_vol;
}

inline void set_field_selection(FieldSelection& fsel,
                                const FieldM<int64_t, 1>& f_rank,
                                const long n_per_tslice_ = 0,
                                const bool is_limit_on_rank = false)
// call set_n_per_tslice if is_limit_on_rank = true
// otherwise will strictly follow f_rank without constraint of n_per_tslice
{
  TIMER_VERBOSE("set_field_selection(fsel,f_rank,n_per_tslice)");
  fsel.init();
  fsel.f_rank = f_rank;
  if (is_limit_on_rank) {
    set_n_per_tslice(fsel.f_rank, n_per_tslice_);
  }
  update_field_selection(fsel);
  update_field_selection(fsel, n_per_tslice_);
}

inline void set_field_selection(FieldSelection& fsel,
                                const Coordinate& total_site)
// select everything with rank = 0
{
  TIMER_VERBOSE("set_field_selection(fsel,total_site)");
  fsel.init();
  mk_field_selection(fsel.f_rank, total_site);
  update_field_selection(fsel);
  update_field_selection(fsel, -1);  // select all points
}

inline bool is_matching_fsel(const FieldSelection& fsel1,
                             const FieldSelection& fsel2)
// only check selection, does not check rank or parameter
{
  const long n_elems = fsel1.n_elems;
  if (n_elems != fsel2.n_elems) {
    return false;
  }
  bool is_same = true;
  const vector_acc<long>& indices1 = fsel1.indices;
  const vector_acc<long>& indices2 = fsel2.indices;
  qassert(indices1.size() == n_elems);
  qassert(indices2.size() == n_elems);
  qthread_for(idx, n_elems, {
    if (indices1[idx] != indices2[idx]) {
      is_same = false;
    }
  });
  return is_same;
}

inline PointSelection psel_from_fsel(const FieldSelection& fsel)
{
  TIMER("psel_from_fsel")
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate total_site = geo.total_site();
  long n_elems = fsel.n_elems;
  long total_n_elems = n_elems;
  glb_sum(total_n_elems);
  //const int num_node = geo.geon.num_node;
  const int id_node = geo.geon.id_node;
  vector<long> vec(geo.geon.num_node, 0);
  all_gather(get_data(vec), get_data_one_elem(n_elems));
  long idx_offset = 0;
  for (int i = 0; i < id_node; ++i) {
    idx_offset += vec[i];
  }
  qassert(idx_offset <= total_n_elems);
  vector<long> vec_gindex(total_n_elems, 0);
  qthread_for(idx, fsel.n_elems, {
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const long gindex = index_from_coordinate(xg, total_site);
    vec_gindex[idx_offset + idx] = gindex;
  });
  glb_sum(get_data(vec_gindex));
  PointSelection psel(total_n_elems);
  qthread_for(idx, (long)psel.size(), {
    long gindex = vec_gindex[idx];
    psel[idx] = coordinate_from_index(gindex, total_site);
  });
  return psel;
}

inline PointSelection psel_from_fsel_local(const FieldSelection& fsel)
{
  TIMER("psel_from_fsel_local")
  const Geometry& geo = fsel.f_rank.geo();
  //const Coordinate total_site = geo.total_site();
  long n_elems = fsel.n_elems;
  PointSelection psel(n_elems);
  qthread_for(idx, (long)psel.size(), {
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    psel[idx] = xg;
  });
  return psel;
}

template <class M, class N>
SelectedField<M>& qcast(SelectedField<N>& x)
// IMPORTANT: will modify the multiplicity of x, need to cast back after finish.
{
  if (x.initialized) {
    const int size = x.geo().multiplicity * sizeof(N);
    x.geo().multiplicity = size / sizeof(M);
    qassert(x.geo().multiplicity * (int)sizeof(M) == size);
  }
  return (SelectedField<M>&)x;
}

template <class M, class N>
const SelectedField<M>& qcast_const(const SelectedField<N>& x)
// IMPORTANT: will modify the multiplicity of x, need to cast back after finish.
{
  return qcast((SelectedField<N>&)x);
}

template <class M>
bool is_initialized(const SelectedField<M>& sf)
{
  return sf.initialized;
}

template <class M>
bool is_consistent(const SelectedField<M>& sf, const FieldSelection& fsel)
{
  return sf.initialized and sf.n_elems == fsel.n_elems and
         geo_remult(sf.geo(), 1) == fsel.f_local_idx.geo() and
         fsel.f_rank.geo() == fsel.f_local_idx.geo() and
         (long) sf.field.size() == sf.n_elems * (long)sf.geo().multiplicity and
         fsel.f_local_idx.geo().is_only_local;
}

template <class M>
void qswap(SelectedField<M>& f1, SelectedField<M>& f2)
{
  std::swap(f1.initialized, f2.initialized);
  std::swap(f1.n_elems, f2.n_elems);
  qswap(f1.geo, f2.geo);
  qswap(f1.field, f2.field);
}

template <class M>
Vector<M> get_data(const SelectedField<M>& sf)
{
  return get_data(sf.field);
}

template <class M>
void set_zero(SelectedField<M>& sf)
{
  TIMER("set_zero(SelectedField)");
  set_zero(get_data(sf));
}

template <class M>
double qnorm(const SelectedField<M>& sf)
{
  double s = qnorm(sf.field);
  glb_sum(s);
  return s;
}

template <class M>
const SelectedField<M>& operator+=(SelectedField<M>& f,
                                   const SelectedField<M>& f1)
{
  TIMER("sel_field_operator+=");
  if (not f.initialized) {
    f = f1;
  } else {
    qassert(f1.initialized);
    qassert(is_matching_geo_mult(f.geo(), f1.geo()));
    qassert(f.field.size() == f1.field.size());
#pragma omp parallel for
    for (long k = 0; k < (long)f.field.size(); ++k) {
      f.field[k] += f1.field[k];
    }
  }
  return f;
}

template <class M>
const SelectedField<M>& operator-=(SelectedField<M>& f,
                                   const SelectedField<M>& f1)
{
  TIMER("sel_field_operator-=");
  if (not f.initialized) {
    f.init(f1.geo(), f1.n_elems, f1.geo().multiplicity);
    set_zero(f);
    f -= f1;
  } else {
    qassert(f1.initialized);
    qassert(is_matching_geo_mult(f.geo(), f1.geo()));
    qassert(f.field.size() == f1.field.size());
#pragma omp parallel for
    for (long k = 0; k < (long)f.field.size(); ++k) {
      f.field[k] -= f1.field[k];
    }
  }
  return f;
}

template <class M>
const SelectedField<M>& operator*=(SelectedField<M>& f, const double factor)
{
  TIMER("sel_field_operator*=(F,D)");
  qassert(f.initialized);
#pragma omp parallel for
  for (long k = 0; k < (long)f.field.size(); ++k) {
    f.field[k] *= factor;
  }
  return f;
}

template <class M>
const SelectedField<M>& operator*=(SelectedField<M>& f, const Complex factor)
{
  TIMER("sel_field_operator*=(F,C)");
  qassert(f.initialized);
#pragma omp parallel for
  for (long k = 0; k < (long)f.field.size(); ++k) {
    f.field[k] *= factor;
  }
  return f;
}

inline void set_selected_gindex(SelectedField<long>& sfgi,
                                const FieldSelection& fsel)
{
  TIMER_VERBOSE("set_selected_gindex");
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate total_site = geo.total_site();
  sfgi.init(fsel, 1);
  qthread_for(idx, fsel.n_elems, {
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const long gindex = index_from_coordinate(xg, total_site);
    sfgi.get_elem(idx) = gindex;
  });
}

template <class M>
void only_keep_selected_points(Field<M>& f, const FieldSelection& fsel)
{
  TIMER("only_keep_selected_points");
  qassert(f.geo().is_only_local);
  qassert(fsel.f_local_idx.geo().is_only_local);
  qassert(geo_remult(f.geo()) == geo_remult(fsel.f_local_idx.geo()));
  const Geometry& geo = f.geo();
  const FieldM<long, 1>& f_local_idx = fsel.f_local_idx;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const long idx = f_local_idx.get_elems_const(index)[0];
    if (idx < 0) {
      qassert(idx == -1);
      Vector<M> fv = f.get_elems(index);
      set_zero(fv);
    }
  }
}

template <class M>
void set_selected_field(SelectedField<M>& sf, const Field<M>& f,
                        const FieldSelection& fsel)
{
  TIMER("set_selected_field(sf,f,fsel)");
  qassert(f.geo().is_only_local);
  qassert(fsel.f_local_idx.geo().is_only_local);
  qassert(geo_remult(f.geo()) == fsel.f_local_idx.geo());
  const Geometry& geo = f.geo();
  const int multiplicity = geo.multiplicity;
  sf.init(fsel, multiplicity);
  qacc_for(idx, fsel.n_elems, {
    const long index = fsel.indices[idx];
    const Vector<M> fv = f.get_elems_const(index);
    Vector<M> sfv = sf.get_elems(idx);
    for (int m = 0; m < multiplicity; ++m) {
      sfv[m] = fv[m];
    }
  });
}

template <class M>
void set_selected_field(SelectedField<M>& sf, const SelectedField<M>& sf0,
                        const FieldSelection& fsel, const FieldSelection& fsel0)
// Does not clear sf's original value if not assigned
{
  TIMER("set_selected_field(sf,sf0,fsel,fsel0)");
  qassert(sf0.geo().is_only_local);
  qassert(fsel.f_local_idx.geo().is_only_local);
  qassert(fsel0.f_local_idx.geo().is_only_local);
  qassert(geo_remult(sf0.geo()) == fsel0.f_local_idx.geo());
  qassert(geo_remult(sf0.geo()) == fsel.f_local_idx.geo());
  const Geometry& geo = sf0.geo();
  const int multiplicity = geo.multiplicity;
  sf.init(fsel, multiplicity);
  qacc_for(idx, fsel.n_elems, {
    const long index = fsel.indices[idx];
    const long idx0 = fsel0.f_local_idx.get_elem(index);
    if (idx0 >= 0) {
      Vector<M> sfv = sf.get_elems(idx);
      const Vector<M> fv = sf0.get_elems_const(idx0);
      for (int m = 0; m < multiplicity; ++m) {
        sfv[m] = fv[m];
      }
    }
  });
}

template <class M>
void set_selected_field(SelectedField<M>& sf, const SelectedPoints<M>& sp,
                        const FieldSelection& fsel, const PointSelection& psel)
// Does not clear sf's original value if not assigned
{
  TIMER("set_selected_field(sf,sp,fsel,psel)");
  qassert(fsel.f_local_idx.geo().is_only_local);
  qassert(geo_remult(sf.geo()) == fsel.f_local_idx.geo());
  const long n_points = sp.n_points;
  qassert(n_points == (long)psel.size());
  const Geometry& geo = fsel.f_rank.geo();
  const int multiplicity = sp.multiplicity;
  sf.init(fsel, multiplicity);
  qthread_for(idx, n_points, {
    const Coordinate& xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      const long sf_idx = fsel.f_local_idx.get_elem(xl);
      if (sf_idx >= 0) {
        qassert(sf_idx < sf.n_elems);
        const Vector<M> spv = sp.get_elems_const(idx);
        Vector<M> fv = sf.get_elems(sf_idx);
        for (int m = 0; m < geo.multiplicity; ++m) {
          fv[m] = spv[m];
        }
      }
    }
  });
}

template <class M>
void set_selected_points(SelectedPoints<M>& sp, const SelectedField<M>& sf,
                         const PointSelection& psel, const FieldSelection& fsel)
// only assign available points
{
  TIMER("set_selected_points(sp,sf,psel,fsel)");
  const Geometry& geo = sf.geo();
  qassert(is_consistent(sf, fsel));
  const long n_points = psel.size();
  sp.init(psel, geo.multiplicity);
  set_zero(sp.points);
  qthread_for(idx, n_points, {
    const Coordinate& xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      const long sf_idx = fsel.f_local_idx.get_elem(xl);
      if (sf_idx >= 0) {
        qassert(sf_idx < sf.n_elems);
        const Vector<M> fv = sf.get_elems_const(sf_idx);
        Vector<M> spv = sp.get_elems(idx);
        for (int m = 0; m < geo.multiplicity; ++m) {
          spv[m] = fv[m];
        }
      }
    }
  });
  glb_sum_byte_vec(get_data(sp.points));
}

template <class M>
void set_field_selected(Field<M>& f, const SelectedField<M>& sf,
                        const FieldSelection& fsel,
                        const bool is_keeping_data = false)
{
  TIMER("set_field_selected(f,sf,fsel)");
  qassert(sf.geo().is_only_local);
  qassert(fsel.f_local_idx.geo().is_only_local);
  qassert(geo_remult(sf.geo()) == fsel.f_local_idx.geo());
  const Geometry& geo = sf.geo();
  if (not is_keeping_data) {
    f.init();
    f.init(sf.geo());
    set_zero(f);
  }
  const int multiplicity = geo.multiplicity;
  qacc_for(idx, fsel.n_elems, {
    const long index = fsel.indices[idx];
    Vector<M> fv = f.get_elems(index);
    const Vector<M> sfv = sf.get_elems_const(idx);
    for (int m = 0; m < multiplicity; ++m) {
      fv[m] = sfv[m];
    }
  });
}

template <class M>
bool is_consistent(const SelectedPoints<M>& sp, const SelectedField<M>& sf,
                   const PointSelection& psel, const FieldSelection& fsel)
{
  TIMER("is_consistent(sp,sf)");
  qassert(is_consistent(sp, psel));
  qassert(is_consistent(sf, fsel));
  const Geometry& geo = sf.geo();
  const long n_points = psel.size();
  double qnorm_diff = 0.0;
  qfor(idx, n_points, {
    const Coordinate& xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      const long sf_idx = fsel.f_local_idx.get_elem(xl);
      if (sf_idx >= 0) {
        const Vector<M> fv = sf.get_elems_const(sf_idx);
        const Vector<M> spv = sp.get_elems_const(idx);
        for (int m = 0; m < geo.multiplicity; ++m) {
          qnorm_diff += qnorm(spv[m] - fv[m]);
        }
      }
    }
  });
  glb_sum(qnorm_diff);
  return qnorm_diff == 0.0;
}

template <class M>
void acc_field(Field<M>& f, const Complex& coef, const SelectedField<M>& sf,
               const FieldSelection& fsel)
// f can be empty
{
  TIMER("acc_field(f,coef,sf,fsel)");
  const Geometry& geo = fsel.f_rank.geo();
  const int multiplicity = sf.geo().multiplicity;
  if (not is_initialized(f)) {
    f.init(geo_remult(geo, multiplicity));
    set_zero(f);
  }
  qassert(multiplicity == f.geo().multiplicity);
  qassert(sf.n_elems == fsel.n_elems);
  qacc_for(idx, fsel.n_elems, {
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> fv = f.get_elems(xl);
    const Vector<M> sfv = sf.get_elems_const(idx);
    for (int m = 0; m < multiplicity; ++m) {
      M x = sfv[m];
      x *= coef;
      fv[m] += x;
    }
  });
}

template <class M>
void acc_field(Field<M>& f, const SelectedField<M>& sf,
               const FieldSelection& fsel)
// f can be empty
{
  TIMER("acc_field(f,sf,fsel)");
  const Geometry& geo = fsel.f_rank.geo();
  const int multiplicity = sf.geo().multiplicity;
  if (not is_initialized(f)) {
    f.init(geo_remult(geo, multiplicity));
    set_zero(f);
  }
  qassert(multiplicity == f.geo().multiplicity);
  qassert(sf.n_elems == fsel.n_elems);
  qacc_for(idx, fsel.n_elems, {
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<M> fv = f.get_elems(xl);
    const Vector<M> sfv = sf.get_elems_const(idx);
    for (int m = 0; m < multiplicity; ++m) {
      fv[m] += sfv[m];
    }
  });
}

template <class M>
std::vector<M> field_sum_tslice(const SelectedField<M>& sf,
                                const FieldSelection& fsel,
                                const int t_dir = 3)
// length = t_size * multiplicity
{
  TIMER("field_sum_tslice");
  qassert(sf.geo().is_only_local);
  qassert(fsel.f_local_idx.geo().is_only_local);
  qassert(geo_remult(sf.geo()) == fsel.f_local_idx.geo());
  const Geometry& geo = sf.geo();
  const int t_size = geo.total_site()[t_dir];
  const int multiplicity = geo.multiplicity;
  std::vector<M> vec(t_size * multiplicity);
  set_zero(vec);
  for (long idx = 0; idx < fsel.n_elems; ++idx) {
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Vector<M> sfv = sf.get_elems_const(idx);
    for (int m = 0; m < multiplicity; ++m) {
      vec[xg[t_dir] * multiplicity + m] += sfv[m];
    }
  }
  return vec;
}

template <class M>
void field_glb_sum_tslice_double(SelectedPoints<M>& sp,
                                 const SelectedField<M>& sf,
                                 const FieldSelection& fsel,
                                 const int t_dir = 3)
{
  TIMER("field_glb_sum_tslice_double(sp,sf,fsel)");
  sp.init();
  const Geometry& geo = sf.geo();
  const int t_size = geo.total_site()[t_dir];
  const int multiplicity = geo.multiplicity;
  std::vector<M> vec = field_sum_tslice(sf, fsel, t_dir);
  glb_sum_double_vec(get_data(vec));
  sp.init(t_size, multiplicity);
  sp.points = vec;
}

template <class M>
void field_glb_sum_tslice_long(SelectedPoints<M>& sp,
                               const SelectedField<M>& sf,
                               const FieldSelection& fsel, const int t_dir = 3)
{
  TIMER("field_glb_sum_tslice_long(sp,sf,fsel)");
  sp.init();
  const Geometry& geo = sf.geo();
  const int t_size = geo.total_site()[t_dir];
  const int multiplicity = geo.multiplicity;
  std::vector<M> vec = field_sum_tslice(sf, fsel, t_dir);
  glb_sum_long_vec(get_data(vec));
  sp.init(t_size, multiplicity);
  sp.points = vec;
}

template <class M>
void qnorm_field(SelectedField<double>& f, const SelectedField<M>& f1)
{
  TIMER("qnorm_field");
  const Geometry& geo = f1.geo();
  f.init();
  f.init(geo, f1.n_elems, 1);
  qacc_for(idx, f.n_elems, {
    const Vector<M> f1v = f1.get_elems_const(idx);
    f.get_elem(idx) = qnorm(f1v);
  });
}

template <class M>
void set_u_rand_double(SelectedField<M>& sf, const FieldSelection& fsel,
                       const RngState& rs, const double upper = 1.0,
                       const double lower = -1.0)
{
  TIMER("set_u_rand_double(sf,fsel,rs)");
  const Geometry& geo = sf.geo();
  qassert(geo.is_only_local);
  qassert(fsel.f_local_idx.geo().is_only_local);
  qassert(geo_remult(geo) == fsel.f_local_idx.geo());
  qacc_for(idx, fsel.n_elems, {
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const long gindex = geo.g_index_from_g_coordinate(xg);
    RngState rsi = rs.newtype(gindex);
    Vector<M> v = sf.get_elems(idx);
    Vector<double> dv((double*)v.data(), v.data_size() / sizeof(double));
    for (int m = 0; m < dv.size(); ++m) {
      dv[m] = u_rand_gen(rsi, upper, lower);
    }
  });
}

// old code

template <class M>
void set_selected_field_slow(SelectedField<M>& sf, const Field<M>& f,
                             const FieldSelection& fsel)
{
  TIMER("set_selected_field_slow");
  qassert(f.geo().is_only_local);
  qassert(fsel.f_local_idx.geo().is_only_local);
  qassert(geo_remult(f.geo()) == fsel.f_local_idx.geo());
  const Geometry& geo = f.geo();
  const int multiplicity = geo.multiplicity;
  sf.init(fsel, multiplicity);
  const FieldM<long, 1>& f_local_idx = fsel.f_local_idx;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const long idx = f_local_idx.get_elems_const(index)[0];
    if (idx >= 0) {
      qassert(idx < fsel.n_elems);
      const long offset = idx * multiplicity;
      const Vector<M> fv = f.get_elems_const(index);
      for (int m = 0; m < multiplicity; ++m) {
        sf.field[offset + m] = fv[m];
      }
    }
  }
}

template <class M>
void set_field_selected_slow(Field<M>& f, const SelectedField<M>& sf,
                             const FieldSelection& fsel)
{
  TIMER("set_field_selected_slow");
  qassert(sf.geo().is_only_local);
  qassert(fsel.f_local_idx.geo().is_only_local);
  qassert(geo_remult(sf.geo()) == fsel.f_local_idx.geo());
  const Geometry& geo = sf.geo();
  f.init();
  f.init(sf.geo());
  set_zero(f);
  const int multiplicity = geo.multiplicity;
  const FieldM<long, 1>& f_local_idx = fsel.f_local_idx;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const long idx = f_local_idx.get_elems_const(index)[0];
    if (idx >= 0) {
      qassert(idx < fsel.n_elems);
      const long offset = idx * multiplicity;
      Vector<M> fv = f.get_elems(index);
      for (int m = 0; m < multiplicity; ++m) {
        fv[m] = sf.field[offset + m];
      }
    }
  }
}

}  // namespace qlat

