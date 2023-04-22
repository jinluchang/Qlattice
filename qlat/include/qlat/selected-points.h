#pragma once

#include <qlat/field.h>

namespace qlat
{  //

PointSelection mk_tslice_point_selection(const int t_size, const int t_dir = 3);

PointSelection mk_tslice_point_selection(const Coordinate& total_site,
                                         const int t_dir = 3);

PointSelection mk_random_point_selection(const Coordinate& total_site,
                                         const long num, const RngState& rs,
                                         const long pool_factor = 2);

void save_point_selection(const PointSelection& psel, const std::string& path);

void save_point_selection_info(const PointSelection& psel,
                               const std::string& path);

PointSelection load_point_selection(const std::string& path);

PointSelection load_point_selection_info(const std::string& path);

// -----------------------

template <class M, class N>
SelectedPoints<M>& qcast(SelectedPoints<N>& x)
// IMPORTANT: will modify the multiplicity of x, need to cast back after finish.
{
  if (x.initialized) {
    const int size = x.multiplicity * sizeof(N);
    x.multiplicity = size / sizeof(M);
    qassert(x.multiplicity * sizeof(M) == size);
  }
  return (SelectedPoints<M>&)x;
}

template <class M, class N>
const SelectedPoints<M>& qcast_const(const SelectedPoints<N>& x)
// IMPORTANT: will modify the multiplicity of x, need to cast back after finish.
{
  return qcast<M, N>((SelectedPoints<N>&)x);
}

template <class M>
bool is_initialized(const SelectedPoints<M>& sp)
{
  return sp.initialized;
}

template <class M>
bool is_consistent(const SelectedPoints<M>& sp, const PointSelection& psel)
{
  return sp.initialized and sp.n_points == (long)psel.size();
}

template <class M>
void qswap(SelectedPoints<M>& f1, SelectedPoints<M>& f2)
{
  std::swap(f1.initialized, f2.initialized);
  std::swap(f1.n_points, f2.n_points);
  std::swap(f1.multiplicity, f2.multiplicity);
  qswap(f1.points, f2.points);
}

template <class M>
const SelectedPoints<M>& operator+=(SelectedPoints<M>& f,
                                    const SelectedPoints<M>& f1)
{
  TIMER("sel_points_operator+=");
  if (not f.initialized) {
    f = f1;
  } else {
    qassert(f1.initialized);
    qassert(f.multiplicity == f1.multiplicity);
    qassert(f.n_points == f1.n_points);
    qassert(f.points.size() == f1.points.size());
#pragma omp parallel for
    for (long k = 0; k < (long)f.points.size(); ++k) {
      f.points[k] += f1.points[k];
    }
  }
  return f;
}

template <class M>
const SelectedPoints<M>& operator-=(SelectedPoints<M>& f,
                                    const SelectedPoints<M>& f1)
{
  TIMER("sel_points_operator-=");
  if (not f.initialized) {
    f.init(f1.n_points, f1.multiplicity);
    set_zero(f);
    f -= f1;
  } else {
    qassert(f1.initialized);
    qassert(f.multiplicity == f1.multiplicity);
    qassert(f.n_points == f1.n_points);
    qassert(f.points.size() == f1.points.size());
#pragma omp parallel for
    for (long k = 0; k < (long)f.points.size(); ++k) {
      f.points[k] -= f1.points[k];
    }
  }
  return f;
}

template <class M>
const SelectedPoints<M>& operator*=(SelectedPoints<M>& f, const double factor)
{
  TIMER("sel_points_operator*=(F,D)");
  qassert(f.initialized);
#pragma omp parallel for
  for (long k = 0; k < (long)f.points.size(); ++k) {
    f.points[k] *= factor;
  }
  return f;
}

template <class M>
const SelectedPoints<M>& operator*=(SelectedPoints<M>& f, const Complex factor)
{
  TIMER("sel_points_operator*=(F,C)");
  qassert(f.initialized);
#pragma omp parallel for
  for (long k = 0; k < (long)f.points.size(); ++k) {
    f.points[k] *= factor;
  }
  return f;
}

template <class M>
void only_keep_selected_points(Field<M>& f, const PointSelection& psel)
{
  TIMER("only_keep_selected_points");
  const Geometry& geo = f.geo();
  qassert(geo.is_only_local);
  Field<M> f1;
  f1.init(geo);
  set_zero(f1);
  const long n_points = psel.size();
#pragma omp parallel for
  for (long idx = 0; idx < n_points; ++idx) {
    const Coordinate& xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      const Vector<M> fv = f.get_elems_const(xl);
      Vector<M> f1v = f1.get_elems(xl);
      for (int m = 0; m < geo.multiplicity; ++m) {
        f1v[m] = fv[m];
      }
    }
  }
  qswap(f, f1);
}

template <class M>
double qnorm(const SelectedPoints<M>& sp)
{
  TIMER("qnorm");
  return qnorm(sp.points);
}

template <class M>
void qnorm_field(SelectedPoints<double>& sp, const SelectedPoints<M>& sp1)
{
  TIMER("qnorm_field");
  sp.init();
  sp.init(sp1.n_points, 1);
  qacc_for(idx, sp.n_points,
           { sp.get_elem(idx) = qnorm(sp1.get_elems_const(idx)); });
}

template <class M>
void set_selected_points(SelectedPoints<M>& sp, const Field<M>& f,
                         const PointSelection& psel)
{
  TIMER("set_selected_points(sp,f,psel)");
  const Geometry& geo = f.geo();
  qassert(geo.is_only_local);
  const long n_points = psel.size();
  sp.init(psel, geo.multiplicity);
  set_zero(sp.points);
  qthread_for(idx, n_points, {
    const Coordinate& xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      const Vector<M> fv = f.get_elems_const(xl);
      Vector<M> spv = sp.get_elems(idx);
      for (int m = 0; m < geo.multiplicity; ++m) {
        spv[m] = fv[m];
      }
    }
  });
  glb_sum_byte_vec(get_data(sp.points));
}

template <class M>
void set_field_selected(Field<M>& f, const SelectedPoints<M>& sp,
                        const Geometry& geo_, const PointSelection& psel)
{
  TIMER("set_field_selected");
  const Geometry geo = geo_reform(geo_, sp.multiplicity, 0);
  qassert(geo.is_only_local);
  const long n_points = sp.n_points;
  qassert(n_points == (long)psel.size());
  f.init();
  f.init(geo);
  set_zero(f);
#pragma omp parallel for
  for (long idx = 0; idx < n_points; ++idx) {
    const Coordinate& xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      const Vector<M> spv = sp.get_elems_const(idx);
      Vector<M> fv = f.get_elems(xl);
      for (int m = 0; m < geo.multiplicity; ++m) {
        fv[m] = spv[m];
      }
    }
  }
}

template <class M>
void set_field_selected(Field<M>& f, const SelectedPoints<M>& sp,
                        const PointSelection& psel)
// deprecated
{
  TIMER("set_field_selected");
  const Geometry& geo = f.geo();
  qassert(geo.multiplicity == sp.multiplicity);
  const long n_points = sp.n_points;
  qassert(n_points == (long)psel.size());
#pragma omp parallel for
  for (long idx = 0; idx < n_points; ++idx) {
    const Coordinate& xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      const Vector<M> spv = sp.get_elems_const(idx);
      Vector<M> fv = f.get_elems(xl);
      for (int m = 0; m < geo.multiplicity; ++m) {
        fv[m] = spv[m];
      }
    }
  }
}

template <class M>
void acc_field(Field<M>& f, const SelectedPoints<M>& sp, const Geometry& geo_,
               const PointSelection& psel)
{
  TIMER("acc_field");
  const Geometry geo = geo_reform(geo_, sp.multiplicity, 0);
  qassert(geo.is_only_local);
  const long n_points = sp.n_points;
  qassert(n_points == (long)psel.size());
  f.init(geo);
#pragma omp parallel for
  for (long idx = 0; idx < n_points; ++idx) {
    const Coordinate& xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      const Vector<M> spv = sp.get_elems_const(idx);
      Vector<M> fv = f.get_elems(xl);
      for (int m = 0; m < geo.multiplicity; ++m) {
        fv[m] += spv[m];
      }
    }
  }
}

template <class M>
void field_glb_sum_tslice_double(SelectedPoints<M>& sp, const Field<M>& f,
                                 const int t_dir = 3)
{
  TIMER("field_glb_sum_tslice_double(sp,f)");
  sp.init();
  const Geometry& geo = f.geo();
  const int t_size = geo.total_site()[t_dir];
  const int multiplicity = geo.multiplicity;
  std::vector<M> vec = field_sum_tslice(f, t_dir);
  glb_sum_double_vec(get_data(vec));
  sp.init(t_size, multiplicity);
  sp.points = vec;
}

template <class M>
void field_glb_sum_tslice_long(SelectedPoints<M>& sp, const Field<M>& f,
                               const int t_dir = 3)
{
  TIMER("field_glb_sum_tslice_long(sp,f)");
  sp.init();
  const Geometry& geo = f.geo();
  const int t_size = geo.total_site()[t_dir];
  const int multiplicity = geo.multiplicity;
  std::vector<M> vec = field_sum_tslice(f, t_dir);
  glb_sum_long_vec(get_data(vec));
  sp.init(t_size, multiplicity);
  sp.points = vec;
}

template <class M>
LatData lat_data_from_selected_points_complex(const SelectedPoints<M>& sp)
{
  TIMER("lat_data_from_selected_points_complex");
  LatData ld;
  ld.info.push_back(lat_dim_number("idx", 0, sp.n_points - 1));
  ld.info.push_back(lat_dim_number("m", 0, sp.multiplicity - 1));
  qassert(sizeof(M) >= sizeof(Complex));
  ld.info.push_back(lat_dim_number("v", 0, (long)(sizeof(M) / sizeof(Complex)) - 1));
  ld.info.push_back(lat_dim_re_im());
  lat_data_alloc(ld);
  assign(get_data(ld.res), get_data(sp.points));
  return ld;
}

template <class M>
void selected_points_from_lat_data_complex(SelectedPoints<M>& sp,
                                           const LatData& ld)
{
  TIMER("selected_points_from_lat_data_complex");
  qassert(ld.info.size() == 4);
  qassert(ld.info[0].name == "idx");
  qassert(ld.info[1].name == "m");
  qassert(ld.info[2].name == "v");
  qassert(ld.info[3].name == "re-im");
  const long n_points = ld.info[0].size;
  const long multiplicity = ld.info[1].size;
  const long sizof_M_vs_sizeof_complex = ld.info[2].size;
  qassert(sizeof(M) == sizof_M_vs_sizeof_complex * sizeof(Complex));
  qassert(ld.info[3].size == 2);
  sp.init(n_points, multiplicity);
  assign(get_data(sp.points), get_data(ld.res));
}

template <class M>
void save_selected_points_complex(const SelectedPoints<M>& sp,
                                  const std::string& path)
{
  TIMER_VERBOSE("save_selected_points_complex");
  if (get_id_node() == 0) {
    const LatData ld = lat_data_from_selected_points_complex(sp);
    ld.save(path);
  }
}

template <class M>
void load_selected_points_complex(SelectedPoints<M>& sp,
                                  const std::string& path)
{
  TIMER_VERBOSE("load_selected_points_complex");
  long n_points = 0;
  long multiplicity = 0;
  if (get_id_node() == 0) {
    LatData ld;
    ld.load(path);
    selected_points_from_lat_data_complex(sp, ld);
    n_points = sp.n_points;
    multiplicity = sp.multiplicity;
  }
  bcast(get_data_one_elem(n_points));
  bcast(get_data_one_elem(multiplicity));
  if (get_id_node() != 0) {
    sp.init(n_points, multiplicity);
  }
  vector<M> buffer(sp.points.size());
  assign(get_data(buffer), get_data(sp.points));
  bcast(get_data(buffer));
  assign(get_data(sp.points), get_data(buffer));
}

// --------------------

#ifdef QLAT_INSTANTIATE_SELECTED_POINTS
#define QLAT_EXTERN
#else
#define QLAT_EXTERN extern
#endif

#define QLAT_EXTERN_TEMPLATE(TYPENAME)                                       \
                                                                             \
  QLAT_EXTERN template void qswap<TYPENAME>(SelectedPoints<TYPENAME> & f1,   \
                                            SelectedPoints<TYPENAME> & f2);  \
                                                                             \
  QLAT_EXTERN template const SelectedPoints<TYPENAME>& operator+=<TYPENAME>( \
      SelectedPoints<TYPENAME>& f, const SelectedPoints<TYPENAME>& f1);      \
                                                                             \
  QLAT_EXTERN template const SelectedPoints<TYPENAME>& operator-=<TYPENAME>( \
      SelectedPoints<TYPENAME>& f, const SelectedPoints<TYPENAME>& f1);      \
                                                                             \
  QLAT_EXTERN template const SelectedPoints<TYPENAME>& operator*=            \
      <TYPENAME>(SelectedPoints<TYPENAME>& f, const double factor);          \
                                                                             \
  QLAT_EXTERN template const SelectedPoints<TYPENAME>& operator*=            \
      <TYPENAME>(SelectedPoints<TYPENAME>& f, const Complex factor);         \
                                                                             \
  QLAT_EXTERN template void only_keep_selected_points<TYPENAME>(             \
      Field<TYPENAME> & f, const PointSelection& psel);                      \
                                                                             \
  QLAT_EXTERN template double qnorm<TYPENAME>(                               \
      const SelectedPoints<TYPENAME>& sp);                                   \
                                                                             \
  QLAT_EXTERN template void qnorm_field<TYPENAME>(                           \
      SelectedPoints<double> & sp, const SelectedPoints<TYPENAME>& sp1);     \
                                                                             \
  QLAT_EXTERN template void set_selected_points<TYPENAME>(                   \
      SelectedPoints<TYPENAME> & sp, const Field<TYPENAME>& f,               \
      const PointSelection& psel);                                           \
                                                                             \
  QLAT_EXTERN template void set_field_selected<TYPENAME>(                    \
      Field<TYPENAME> & f, const SelectedPoints<TYPENAME>& sp,               \
      const Geometry& geo_, const PointSelection& psel);                     \
                                                                             \
  QLAT_EXTERN template void set_field_selected<TYPENAME>(                    \
      Field<TYPENAME> & f, const SelectedPoints<TYPENAME>& sp,               \
      const PointSelection& psel);                                           \
                                                                             \
  QLAT_EXTERN template void acc_field<TYPENAME>(                             \
      Field<TYPENAME> & f, const SelectedPoints<TYPENAME>& sp,               \
      const Geometry& geo_, const PointSelection& psel);                     \
                                                                             \
  QLAT_EXTERN template void field_glb_sum_tslice_double<TYPENAME>(           \
      SelectedPoints<TYPENAME> & sp, const Field<TYPENAME>& f,               \
      const int t_dir);                                                      \
                                                                             \
  QLAT_EXTERN template void field_glb_sum_tslice_long<TYPENAME>(             \
      SelectedPoints<TYPENAME> & sp, const Field<TYPENAME>& f,               \
      const int t_dir);                                                      \
                                                                             \
  QLAT_EXTERN template LatData                                               \
  lat_data_from_selected_points_complex<TYPENAME>(                           \
      const SelectedPoints<TYPENAME>& sp);                                   \
                                                                             \
  QLAT_EXTERN template void selected_points_from_lat_data_complex<TYPENAME>( \
      SelectedPoints<TYPENAME> & sp, const LatData& ld);                     \
                                                                             \
  QLAT_EXTERN template void save_selected_points_complex<TYPENAME>(          \
      const SelectedPoints<TYPENAME>& sp, const std::string& path);          \
                                                                             \
  QLAT_EXTERN template void load_selected_points_complex<TYPENAME>(          \
      SelectedPoints<TYPENAME> & sp, const std::string& path)

QLAT_CALL_WITH_TYPES(QLAT_EXTERN_TEMPLATE);

#undef QLAT_EXTERN
#undef QLAT_EXTERN_TEMPLATE
#undef QLAT_EXTERN_CLASS

}  // namespace qlat
