#pragma once

#include <qlat/field.h>

namespace qlat
{  //

PointsSelection mk_tslice_point_selection(const int t_size, const int t_dir = 3);

PointsSelection mk_tslice_point_selection(const Coordinate& total_site,
                                          const int t_dir = 3);

PointsSelection mk_random_point_selection(const Coordinate& total_site,
                                         const Long num, const RngState& rs,
                                         const Long pool_factor = 2);

void save_point_selection(const PointsSelection& psel, const std::string& path);

void save_point_selection_info(const PointsSelection& psel,
                               const std::string& path);

PointsSelection load_point_selection(const std::string& path);

PointsSelection load_point_selection_info(const std::string& path);

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
bool is_consistent(const SelectedPoints<M>& sp, const PointsSelection& psel)
{
  return sp.initialized and sp.n_points == (Long)psel.size();
}

template <class M>
SelectedPoints<M>& operator+=(SelectedPoints<M>& f, const SelectedPoints<M>& f1)
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
    for (Long k = 0; k < (Long)f.points.size(); ++k) {
      f.points[k] += f1.points[k];
    }
  }
  return f;
}

template <class M>
SelectedPoints<M>& operator-=(SelectedPoints<M>& f, const SelectedPoints<M>& f1)
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
    for (Long k = 0; k < (Long)f.points.size(); ++k) {
      f.points[k] -= f1.points[k];
    }
  }
  return f;
}

template <class M>
SelectedPoints<M>& operator*=(SelectedPoints<M>& f, const double factor)
{
  TIMER("sel_points_operator*=(F,D)");
  qassert(f.initialized);
#pragma omp parallel for
  for (Long k = 0; k < (Long)f.points.size(); ++k) {
    f.points[k] *= factor;
  }
  return f;
}

template <class M>
SelectedPoints<M>& operator*=(SelectedPoints<M>& f, const ComplexD factor)
{
  TIMER("sel_points_operator*=(F,C)");
  qassert(f.initialized);
#pragma omp parallel for
  for (Long k = 0; k < (Long)f.points.size(); ++k) {
    f.points[k] *= factor;
  }
  return f;
}

template <class M>
void only_keep_selected_points(Field<M>& f, const PointsSelection& psel)
{
  TIMER("only_keep_selected_points");
  const Geometry& geo = f.geo();
  qassert(geo.is_only_local);
  Field<M> f1;
  f1.init(geo);
  set_zero(f1);
  const Long n_points = psel.size();
#pragma omp parallel for
  for (Long idx = 0; idx < n_points; ++idx) {
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
RealD qnorm(const SelectedPoints<M>& sp)
{
  TIMER("qnorm(sp)");
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

void set_sqrt_field(SelectedPoints<RealD>& sp,
                    const SelectedPoints<RealD>& sp1);

// -------------------------------------------

template <class M>
void set_selected_points(SelectedPoints<M>& sp, const Field<M>& f,
                         const PointsSelection& psel)
{
  TIMER("set_selected_points(sp,f,psel)");
  const Geometry& geo = f.geo();
  qassert(geo.is_only_local);
  const Long n_points = psel.size();
  sp.init(psel, geo.multiplicity);
  set_zero(sp);  // has to set_zero for glb_sum_byte_vec
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
void set_selected_points(SelectedPoints<M>& sp, const Field<M>& f,
                         const PointsSelection& psel, const int m)
{
  TIMER("set_selected_points(sp,f,psel,m)");
  const Geometry& geo = f.geo();
  qassert(geo.is_only_local);
  const Long n_points = psel.size();
  sp.init(psel, 1);
  set_zero(sp);  // has to set_zero for glb_sum_byte_vec
  qthread_for(idx, n_points, {
    const Coordinate& xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      const Vector<M> fv = f.get_elems_const(xl);
      sp.get_elem(idx) = fv[m];
    }
  });
  glb_sum_byte_vec(get_data(sp.points));
}

template <class M>
void set_field_selected(Field<M>& f, const SelectedPoints<M>& sp,
                        const Geometry& geo_, const PointsSelection& psel,
                        const bool is_keeping_data = true)
{
  TIMER("set_field_selected(f,sp,geo,psel)");
  const Geometry geo = geo_reform(geo_, sp.multiplicity, 0);
  qassert(geo.is_only_local);
  const Long n_points = sp.n_points;
  qassert(n_points == (Long)psel.size());
  if (is_keeping_data) {
    f.init_zero(geo);
  } else {
    f.init(geo);
    set_zero(f);
  }
  qthread_for(idx, n_points, {
    const Coordinate& xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      const Vector<M> spv = sp.get_elems_const(idx);
      Vector<M> fv = f.get_elems(xl);
      for (int m = 0; m < geo.multiplicity; ++m) {
        fv[m] = spv[m];
      }
    }
  });
}

template <class M>
void set_field_selected(Field<M>& f, const SelectedPoints<M>& sp,
                        const Geometry& geo_, const PointsSelection& psel,
                        const int m, const bool is_keeping_data = true)
{
  TIMER("set_field_selected(f,sp,geo,psel,m)");
  if (is_keeping_data) {
    f.init_zero(geo_);
  } else {
    f.init(geo_);
    set_zero(f);
  }
  const Geometry& geo = f.geo();
  const Long n_points = sp.n_points;
  qassert(n_points == (Long)psel.size());
  qassert(sp.multiplicity == 1);
  qthread_for(idx, n_points, {
    const Coordinate& xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      Vector<M> fv = f.get_elems(xl);
      fv[m] = sp.get_elem(idx);
    }
  });
}

// -------------------------------------------

template <class M>
void acc_field(Field<M>& f, const SelectedPoints<M>& sp, const Geometry& geo_,
               const PointsSelection& psel)
{
  TIMER("acc_field(f,sp,geo,psel)");
  const Geometry geo = geo_reform(geo_, sp.multiplicity, 0);
  qassert(geo.is_only_local);
  const Long n_points = sp.n_points;
  qassert(n_points == (Long)psel.size());
  f.init(geo);
#pragma omp parallel for
  for (Long idx = 0; idx < n_points; ++idx) {
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
void field_glb_sum(SelectedPoints<M>& sp, const Field<M>& f)
{
  TIMER("field_glb_sum(sp,f)");
  sp.init();
  const Geometry& geo = f.geo();
  const int multiplicity = geo.multiplicity;
  std::vector<M> vec = field_glb_sum(f);
  sp.init(1, multiplicity);
  sp.points = vec;
}

template <class M>
void field_glb_sum_tslice(SelectedPoints<M>& sp, const Field<M>& f,
                          const int t_dir = 3)
{
  TIMER("field_glb_sum_tslice(sp,f)");
  sp.init();
  const Geometry& geo = f.geo();
  const int t_size = geo.total_site()[t_dir];
  const int multiplicity = geo.multiplicity;
  std::vector<M> vec = field_glb_sum_tslice(f, t_dir);
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
  qassert(sizeof(M) >= sizeof(ComplexD));
  ld.info.push_back(lat_dim_number("v", 0, (Long)(sizeof(M) / sizeof(ComplexD)) - 1));
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
  const Long n_points = ld.info[0].size;
  const Long multiplicity = ld.info[1].size;
  const Long sizof_M_vs_sizeof_complex = ld.info[2].size;
  qassert(sizeof(M) == sizof_M_vs_sizeof_complex * sizeof(ComplexD));
  qassert(ld.info[3].size == 2);
  sp.init(n_points, multiplicity);
  assign(get_data(sp.points), get_data(ld.res));
}

template <class M>
LatData lat_data_from_selected_points_double(const SelectedPoints<M>& sp)
{
  TIMER("lat_data_from_selected_points_double");
  LatData ld;
  ld.info.push_back(lat_dim_number("idx", 0, sp.n_points - 1));
  ld.info.push_back(lat_dim_number("m", 0, sp.multiplicity - 1));
  qassert(sizeof(M) >= sizeof(double));
  ld.info.push_back(
      lat_dim_number("v", 0, (Long)(sizeof(M) / sizeof(double)) - 1));
  lat_data_alloc(ld);
  assign(get_data(ld.res), get_data(sp.points));
  return ld;
}

template <class M>
void selected_points_from_lat_data_double(SelectedPoints<M>& sp,
                                          const LatData& ld)
{
  TIMER("selected_points_from_lat_data_double");
  qassert(ld.info.size() == 4);
  qassert(ld.info[0].name == "idx");
  qassert(ld.info[1].name == "m");
  qassert(ld.info[2].name == "v");
  const Long n_points = ld.info[0].size;
  const Long multiplicity = ld.info[1].size;
  const Long sizof_M_vs_sizeof_double = ld.info[2].size;
  qassert(sizeof(M) == sizof_M_vs_sizeof_double * sizeof(double));
  sp.init(n_points, multiplicity);
  assign(get_data(sp.points), get_data(ld.res));
}

template <class M>
LatData lat_data_from_selected_points(const SelectedPoints<M>& sp)
{
  TIMER("lat_data_from_selected_points");
  if (is_composed_of_complex_d<M>()) {
    return lat_data_from_selected_points_complex(sp);
  } else if (is_composed_of_real_d<M>()) {
    return lat_data_from_selected_points_double(sp);
  } else {
    qerr(fname + ssprintf(": get_type_name(M)=%s", get_type_name<M>().c_str()));
    return LatData();
  }
}

template <class M>
void selected_points_from_lat_data(SelectedPoints<M>& sp, const LatData& ld)
{
  TIMER("selected_points_from_lat_data");
  if (is_composed_of_complex_d<M>()) {
    return selected_points_from_lat_data_complex(sp, ld);
  } else if (is_composed_of_real_d<M>()) {
    return selected_points_from_lat_data_double(sp, ld);
  } else {
    qerr(fname + ssprintf(": get_type_name(M)=%s", get_type_name<M>().c_str()));
  }
}

template <class M>
void save_selected_points(const SelectedPoints<M>& sp, const std::string& path)
{
  TIMER_VERBOSE("save_selected_points");
  if (get_id_node() == 0) {
    const LatData ld = lat_data_from_selected_points(sp);
    ld.save(path);
  }
}

template <class M>
void load_selected_points(SelectedPoints<M>& sp, const std::string& path)
{
  TIMER_VERBOSE("load_selected_points");
  Long n_points = 0;
  Long multiplicity = 0;
  if (get_id_node() == 0) {
    LatData ld;
    ld.load(path);
    selected_points_from_lat_data(sp, ld);
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

#define QLAT_EXTERN_TEMPLATE(TYPENAME)                                   \
                                                                         \
  QLAT_EXTERN template SelectedPoints<TYPENAME>& operator+=<TYPENAME>(   \
      SelectedPoints<TYPENAME>& f, const SelectedPoints<TYPENAME>& f1);  \
                                                                         \
  QLAT_EXTERN template SelectedPoints<TYPENAME>& operator-=<TYPENAME>(   \
      SelectedPoints<TYPENAME>& f, const SelectedPoints<TYPENAME>& f1);  \
                                                                         \
  QLAT_EXTERN template SelectedPoints<TYPENAME>& operator*=              \
      <TYPENAME>(SelectedPoints<TYPENAME>& f, const double factor);      \
                                                                         \
  QLAT_EXTERN template SelectedPoints<TYPENAME>& operator*=              \
      <TYPENAME>(SelectedPoints<TYPENAME>& f, const ComplexD factor);    \
                                                                         \
  QLAT_EXTERN template void only_keep_selected_points<TYPENAME>(         \
      Field<TYPENAME> & f, const PointsSelection& psel);                 \
                                                                         \
  QLAT_EXTERN template double qnorm<TYPENAME>(                           \
      const SelectedPoints<TYPENAME>& sp);                               \
                                                                         \
  QLAT_EXTERN template void qnorm_field<TYPENAME>(                       \
      SelectedPoints<double> & sp, const SelectedPoints<TYPENAME>& sp1); \
                                                                         \
  QLAT_EXTERN template void set_selected_points<TYPENAME>(               \
      SelectedPoints<TYPENAME> & sp, const Field<TYPENAME>& f,           \
      const PointsSelection& psel);                                      \
                                                                         \
  QLAT_EXTERN template void set_field_selected<TYPENAME>(                \
      Field<TYPENAME> & f, const SelectedPoints<TYPENAME>& sp,           \
      const Geometry& geo_, const PointsSelection& psel,                 \
      const bool is_keeping_data);                                       \
                                                                         \
  QLAT_EXTERN template void set_field_selected<TYPENAME>(                \
      Field<TYPENAME> & f, const SelectedPoints<TYPENAME>& sp,           \
      const Geometry& geo_, const PointsSelection& psel, const int m,    \
      const bool is_keeping_data);                                       \
                                                                         \
  QLAT_EXTERN template void acc_field<TYPENAME>(                         \
      Field<TYPENAME> & f, const SelectedPoints<TYPENAME>& sp,           \
      const Geometry& geo_, const PointsSelection& psel);                \
                                                                         \
  QLAT_EXTERN template void field_glb_sum_tslice<TYPENAME>(              \
      SelectedPoints<TYPENAME> & sp, const Field<TYPENAME>& f,           \
      const int t_dir);                                                  \
                                                                         \
  QLAT_EXTERN template LatData lat_data_from_selected_points<TYPENAME>(  \
      const SelectedPoints<TYPENAME>& sp);                               \
                                                                         \
  QLAT_EXTERN template void selected_points_from_lat_data<TYPENAME>(     \
      SelectedPoints<TYPENAME> & sp, const LatData& ld);                 \
                                                                         \
  QLAT_EXTERN template void save_selected_points<TYPENAME>(              \
      const SelectedPoints<TYPENAME>& sp, const std::string& path);      \
                                                                         \
  QLAT_EXTERN template void load_selected_points<TYPENAME>(              \
      SelectedPoints<TYPENAME> & sp, const std::string& path)

QLAT_CALL_WITH_TYPES(QLAT_EXTERN_TEMPLATE);
#undef QLAT_EXTERN_TEMPLATE

#undef QLAT_EXTERN

}  // namespace qlat
