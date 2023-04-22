#pragma once

#include <qlat/field.h>
#include <qlat/selected-points.h>

namespace qlat
{  //

void add_field_selection(FieldM<int64_t, 1>& f_rank, const PointSelection& psel,
                         const long rank_psel = 1024L * 1024L * 1024L * 1024L *
                                                1024L);

void mk_field_selection(FieldM<int64_t, 1>& f_rank,
                        const Coordinate& total_site, const int64_t val = 0);

void mk_field_selection(FieldM<int64_t, 1>& f_rank,
                        const Coordinate& total_site,
                        const std::vector<Coordinate>& xgs,
                        const long rank_xgs = 1024L * 1024L * 1024L * 1024L *
                                              1024L);

void select_rank_range(FieldM<int64_t, 1>& f_rank, const long rank_start = 0,
                       const long rank_stop = -1);

void select_t_range(FieldM<int64_t, 1>& f_rank, const long t_start = 0,
                    const long t_stop = -1);

void set_n_per_tslice(FieldM<int64_t, 1>& f_rank, const long n_per_tslice);

void update_field_selection(FieldSelection& fsel);

void update_field_selection(FieldSelection& fsel, const long n_per_tslice_);

void set_grid_field_selection(FieldSelection& fsel,
                              const Coordinate& total_site,
                              const long n_per_tslice, const RngState& rs);

void set_field_selection(FieldSelection& fsel, const FieldM<int64_t, 1>& f_rank,
                         const long n_per_tslice_ = 0,
                         const bool is_limit_on_rank = false);

void set_field_selection(FieldSelection& fsel, const Coordinate& total_site);

void set_field_selection(FieldSelection& fsel, const Coordinate& total_site,
                         const long n_per_tslice, const RngState& rs);

void set_field_selection(FieldSelection& fsel, const Coordinate& total_site,
                         const long n_per_tslice, const RngState& rs,
                         const PointSelection& psel);

bool is_matching_fsel(const FieldSelection& fsel1, const FieldSelection& fsel2);

PointSelection psel_from_fsel(const FieldSelection& fsel);

PointSelection psel_from_fsel_local(const FieldSelection& fsel);

void set_selected_gindex(SelectedField<long>& sfgi, const FieldSelection& fsel);

void mk_grid_field_selection(FieldM<int64_t, 1>& f_rank,
                             const Coordinate& total_site,
                             const long n_per_tslice_, const RngState& rs);

void mk_field_selection(FieldM<int64_t, 1>& f_rank,
                        const Coordinate& total_site, const long n_per_tslice,
                        const RngState& rs);

long write_field_selection(const FieldSelection& fsel, const std::string& path);

long read_field_selection(FieldSelection& fsel, const std::string& path,
                          const long n_per_tslice);

// ---------------------------------------

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
  return qcast<M, N>((SelectedField<N>&)x);
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
                                const FieldSelection& fsel, const int t_dir = 3)
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
  qthread_for(idx, fsel.n_elems, {
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

template <class M, class N>
void convert_field_float_from_double(SelectedField<N>& ff,
                                     const SelectedField<M>& f)
// interface_function
{
  TIMER("convert_field_float_from_double(sf)");
  qassert(f.geo().is_only_local);
  qassert(sizeof(M) % sizeof(double) == 0);
  qassert(sizeof(N) % sizeof(float) == 0);
  qassert(f.geo().multiplicity * sizeof(M) / 2 % sizeof(N) == 0);
  const int multiplicity = f.geo().multiplicity * sizeof(M) / 2 / sizeof(N);
  const Geometry geo = geo_remult(f.geo(), multiplicity);
  const long n_elems = f.n_elems;
  ff.init(geo, n_elems, multiplicity);
  const Vector<M> fdata = get_data(f);
  const Vector<double> fd((double*)fdata.data(),
                          fdata.data_size() / sizeof(double));
  Vector<N> ffdata = get_data(ff);
  Vector<float> ffd((float*)ffdata.data(), ffdata.data_size() / sizeof(float));
  qassert(ffd.size() == fd.size());
  qacc_for(i, ffd.size(), {
    ffd[i] = fd[i];
  });
}

template <class M, class N>
void convert_field_double_from_float(SelectedField<N>& ff,
                                     const SelectedField<M>& f)
// interface_function
{
  TIMER("convert_field_double_from_float(sf)");
  qassert(f.geo().is_only_local);
  qassert(sizeof(M) % sizeof(float) == 0);
  qassert(sizeof(N) % sizeof(double) == 0);
  qassert(f.geo().multiplicity * sizeof(M) * 2 % sizeof(N) == 0);
  const int multiplicity = f.geo().multiplicity * sizeof(M) * 2 / sizeof(N);
  const Geometry geo = geo_remult(f.geo(), multiplicity);
  const long n_elems = f.n_elems;
  ff.init(geo, n_elems, multiplicity);
  const Vector<M> fdata = get_data(f);
  const Vector<float> fd((float*)fdata.data(),
                         fdata.data_size() / sizeof(float));
  Vector<N> ffdata = get_data(ff);
  Vector<double> ffd((double*)ffdata.data(),
                     ffdata.data_size() / sizeof(double));
  qassert(ffd.size() == fd.size());
  qacc_for(i, ffd.size(), {
    ffd[i] = fd[i];
  });
}

}  // namespace qlat
