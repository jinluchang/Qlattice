#pragma once

#include <qlat/field.h>

namespace qlat
{  //

PointsSelection mk_tslice_points_selection(const Coordinate& total_site,
                                          const Int t_dir = 3);

PointsSelection mk_random_points_selection(const Coordinate& total_site,
                                           const Long num, const RngState& rs,
                                           const Long pool_factor = 2);

void set_psel_full(PointsSelection& psel, const Geometry& geo);

void lat_data_from_points_selection(LatDataInt& ld, const PointsSelection& psel);

void points_selection_from_lat_data(PointsSelection& psel, const LatDataInt& ld, const PointsDistType points_dist_type = PointsDistType::Global);

void save_points_selection(const PointsSelection& psel, const std::string& path);

void save_points_selection_info(const PointsSelection& psel,
                               const std::string& path);

PointsSelection load_points_selection(const std::string& path);

PointsSelection load_points_selection_info(const std::string& path);

crc32_t crc32_par(const PointsSelection& psel);

// -----------------------

template <class M>
bool is_initialized(const SelectedPoints<M>& sp)
{
  return sp.initialized;
}

template <class M>
bool is_consistent(const SelectedPoints<M>& sp, const PointsSelection& psel)
{
  return sp.initialized and sp.n_points == (Long)psel.size() and sp.points_dist_type == psel.points_dist_type;
}

template <class M>
SelectedPoints<M>& operator+=(SelectedPoints<M>& f, const SelectedPoints<M>& f1)
{
  TIMER("sel_points_operator+=");
  if (not f.initialized) {
    f = f1;
  } else {
    Qassert(f1.initialized);
    Qassert(f.multiplicity == f1.multiplicity);
    Qassert(f.n_points == f1.n_points);
    Qassert(f.points.size() == f1.points.size());
    qacc_for(k, f.points.size(), { f.points[k] += f1.points[k]; });
  }
  return f;
}

template <class M>
SelectedPoints<M>& operator-=(SelectedPoints<M>& f, const SelectedPoints<M>& f1)
{
  TIMER("sel_points_operator-=");
  if (not f.initialized) {
    f.init(f1.n_points, f1.multiplicity, f1.points_dist_type);
    set_zero(f);
    f -= f1;
  } else {
    Qassert(f1.initialized);
    Qassert(f.multiplicity == f1.multiplicity);
    Qassert(f.n_points == f1.n_points);
    Qassert(f.points.size() == f1.points.size());
    qacc_for(k, f.points.size(), { f.points[k] -= f1.points[k]; });
  }
  return f;
}

template <class M>
SelectedPoints<M>& operator*=(SelectedPoints<M>& f, const double factor)
{
  TIMER("sel_points_operator*=(F,D)");
  Qassert(f.initialized);
  qacc_for(k, f.points.size(), { f.points[k] *= factor; });
  return f;
}

template <class M>
SelectedPoints<M>& operator*=(SelectedPoints<M>& f, const ComplexD factor)
{
  TIMER("sel_points_operator*=(F,C)");
  Qassert(f.initialized);
  qacc_for(k, f.points.size(), { f.points[k] *= factor; });
  return f;
}

template <class M>
void only_keep_selected_points(Field<M>& f, const PointsSelection& psel)
{
  TIMER("only_keep_selected_points");
  const Geometry& geo = f.geo();
  Qassert(geo.is_only_local);
  Field<M> f1;
  f1.init(geo, f.multiplicity);
  set_zero(f1);
  const Long n_points = psel.size();
  qacc_for(idx, n_points, {
    const Coordinate& xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      const Vector<M> fv = f.get_elems_const(xl);
      Vector<M> f1v = f1.get_elems(xl);
      for (Int m = 0; m < f.multiplicity; ++m) {
        f1v[m] = fv[m];
      }
    }
  });
  qswap(f, f1);
}

template <class M>
RealD qnorm(const SelectedPoints<M>& sp)
{
  TIMER("qnorm(sp)");
  return qnorm(sp.points);
}

template <class M>
void qnorm_field(SelectedPoints<RealD>& sp, const SelectedPoints<M>& sp1)
{
  TIMER("qnorm_field(sp,sp1)");
  sp.init(sp1.n_points, 1, sp1.points_dist_type);
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
  Qassert(geo.is_only_local);
  Qassert(psel.points_dist_type == PointsDistType::Global);
  const Long n_points = psel.size();
  sp.init(psel, f.multiplicity);
  set_zero(sp);  // has to set_zero for glb_sum
  qacc_for(idx, n_points, {
    const Coordinate& xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      const Vector<M> fv = f.get_elems_const(xl);
      Vector<M> spv = sp.get_elems(idx);
      for (Int m = 0; m < f.multiplicity; ++m) {
        spv[m] = fv[m];
      }
    }
  });
  glb_sum(get_data_char(sp.points));
}

template <class M>
void set_selected_points(SelectedPoints<M>& sp, const Field<M>& f,
                         const PointsSelection& psel, const Int m)
{
  TIMER("set_selected_points(sp,f,psel,m)");
  const Geometry& geo = f.geo();
  Qassert(geo.is_only_local);
  Qassert(psel.points_dist_type == PointsDistType::Global);
  const Long n_points = psel.size();
  sp.init(psel, 1);
  set_zero(sp);  // has to set_zero for glb_sum
  qacc_for(idx, n_points, {
    const Coordinate& xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      const Vector<M> fv = f.get_elems_const(xl);
      sp.get_elem(idx) = fv[m];
    }
  });
  glb_sum(get_data_char(sp.points));
}

template <class M>
void set_selected_points(SelectedPoints<M>& sp, const SelectedPoints<M>& sp0,
                         const PointsSelection& psel,
                         const PointsSelection& psel0,
                         const bool is_keeping_data = true)
// Most efficient if psel and psel0 is the same.
// If not, more efficient if psel and psel0 has the same order.
{
  if (&sp == &sp0) {
    return;
  }
  TIMER("set_selected_points(sp,sp0,psel,psel0)");
  if (&psel == &psel0) {
    sp = sp0;
    return;
  }
  Qassert(psel.points_dist_type == psel0.points_dist_type);
  const Long n_points = psel.size();
  const Long n_points0 = psel0.size();
  const Int multiplicity = sp0.multiplicity;
  Qassert(sp0.n_points == n_points0);
  if (is_keeping_data) {
    sp.init_zero(psel, multiplicity);
  } else {
    sp.init(psel, multiplicity);
    set_zero(sp);
  }
  Qassert(sp.n_points == n_points);
  SelectedPoints<Long> sp_idx;
  sp_idx.init(psel, 1);
  Qassert(sp_idx.n_points == n_points);
  Long idx_last = -1;
  qfor(idx, n_points, {
    const Coordinate& xg = psel[idx];
    Long idx0 = -1;
    for (Long i = 0; i < n_points0; ++i) {
      idx_last += 1;
      if (idx_last >= n_points0) {
        idx_last = idx_last % n_points0;
      }
      const Coordinate& xg0 = psel0[idx_last];
      if (xg0 == xg) {
        idx0 = idx_last;
        break;
      }
    }
    qassert(idx0 < n_points0);
    sp_idx.get_elem(idx) = idx0;
  });
  qacc_for(idx, n_points, {
    const Long idx0 = sp_idx.get_elem(idx);
    if (idx0 >= 0) {
      qassert(idx0 < n_points0);
      const Vector<M> spv0 = sp0.get_elems_const(idx0);
      Vector<M> spv = sp.get_elems(idx);
      for (Int m = 0; m < multiplicity; ++m) {
        spv[m] = spv0[m];
      }
    }
  });
}

template <class M>
void set_field_selected(Field<M>& f, const SelectedPoints<M>& sp,
                        const Geometry& geo, const PointsSelection& psel,
                        const bool is_keeping_data = true)
{
  TIMER("set_field_selected(f,sp,geo,psel)");
  Qassert(geo.is_only_local);
  const Long n_points = sp.n_points;
  Qassert(n_points == (Long)psel.size());
  if (is_keeping_data) {
    f.init_zero(geo, sp.multiplicity);
  } else {
    f.init(geo, sp.multiplicity);
    set_zero(f);
  }
  qacc_for(idx, n_points, {
    const Coordinate& xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      const Vector<M> spv = sp.get_elems_const(idx);
      Vector<M> fv = f.get_elems(xl);
      for (Int m = 0; m < sp.multiplicity; ++m) {
        fv[m] = spv[m];
      }
    }
  });
}

template <class M>
void set_field_selected(Field<M>& f, const SelectedPoints<M>& sp,
                        const Geometry& geo, const PointsSelection& psel,
                        const Int multiplicity, const Int m,
                        const bool is_keeping_data = true)
{
  TIMER("set_field_selected(f,sp,geo,psel,mult,m)");
  if (is_keeping_data) {
    f.init_zero(geo, multiplicity);
  } else {
    f.init(geo, multiplicity);
    set_zero(f);
  }
  Qassert(0 <= m and m < multiplicity);
  const Long n_points = sp.n_points;
  Qassert(n_points == (Long)psel.size());
  Qassert(sp.multiplicity == 1);
  qacc_for(idx, n_points, {
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
void acc_field(Field<M>& f, const SelectedPoints<M>& sp, const Geometry& geo,
               const PointsSelection& psel)
{
  TIMER("acc_field(f,sp,geo,psel)");
  Qassert(geo.is_only_local);
  const Long n_points = sp.n_points;
  Qassert(n_points == (Long)psel.size());
  f.init(geo, sp.multiplicity);
  qacc_for(idx, n_points, {
    const Coordinate& xg = psel[idx];
    const Coordinate xl = geo.coordinate_l_from_g(xg);
    if (geo.is_local(xl)) {
      const Vector<M> spv = sp.get_elems_const(idx);
      Vector<M> fv = f.get_elems(xl);
      for (Int m = 0; m < sp.multiplicity; ++m) {
        fv[m] += spv[m];
      }
    }
  });
}

template <class M>
void field_glb_sum(SelectedPoints<M>& sp, const Field<M>& f)
{
  TIMER("field_glb_sum(sp,f)");
  sp.init();
  const Int multiplicity = f.multiplicity;
  std::vector<M> vec = field_glb_sum(f);
  sp.init(1, multiplicity, PointsDistType::Global);
  sp.points = vec;
}

template <class M>
void field_glb_sum_tslice(SelectedPoints<M>& sp, const Field<M>& f,
                          const Int t_dir = 3)
{
  TIMER("field_glb_sum_tslice(sp,f)");
  sp.init();
  const Geometry geo = f.get_geo();
  const Int t_size = geo.total_site()[t_dir];
  const Int multiplicity = f.multiplicity;
  std::vector<M> vec = field_glb_sum_tslice(f, t_dir);
  sp.init(t_size, multiplicity, PointsDistType::Global);
  sp.points = vec;
}

template <class M>
void field_glb_max(SelectedPoints<M>& sp, const Field<M>& f)
{
  TIMER("field_glb_max(sp,f)");
  sp.init();
  const Int multiplicity = f.multiplicity;
  std::vector<M> vec = field_glb_max(f);
  sp.init(1, multiplicity, PointsDistType::Global);
  sp.points = vec;
}

template <class M>
void field_glb_min(SelectedPoints<M>& sp, const Field<M>& f)
{
  TIMER("field_glb_min(sp,f)");
  sp.init();
  const Int multiplicity = f.multiplicity;
  std::vector<M> vec = field_glb_min(f);
  sp.init(1, multiplicity, PointsDistType::Global);
  sp.points = vec;
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
void set_u_rand(SelectedPoints<M>& sp, const PointsSelection& psel,
                const RngState& rs, const RealD upper = 1.0,
                const RealD lower = -1.0)
{
  TIMER("set_u_rand(sp,psel,rs)");
  if (not is_composed_of_real<M>()) {
    Qassert(is_composed_of_real<M>());
    return;
  }
  using Real = typename IsDataValueType<M>::ElementaryType;
  const Long n_points = sp.n_points;
  const Int multiplicity = sp.multiplicity;
  const Coordinate total_site = psel.total_site;
  Qassert(n_points == psel.size());
  Qassert(sp.points.size() == n_points * multiplicity);
  qthread_for(idx, n_points, {
    const Coordinate xg = psel[idx];
    const Long gindex = index_from_coordinate(xg, total_site);
    RngState rsi = rs.newtype(gindex);
    Vector<M> v = sp.get_elems(idx);
    Vector<Real> dv((Real*)v.data(), v.data_size() / sizeof(Real));
    for (Int m = 0; m < dv.size(); ++m) {
      dv[m] = u_rand_gen(rsi, upper, lower);
    }
  });
}

template <class M, QLAT_ENABLE_IF(is_data_value_type<M>())>
void set_g_rand(SelectedPoints<M>& sp, const PointsSelection& psel,
                const RngState& rs, const RealD center = 0.0,
                const RealD sigma = 1.0)
{
  TIMER("set_g_rand(sp,psel,rs,center,sigma)");
  if (not is_composed_of_real<M>()) {
    Qassert(is_composed_of_real<M>());
    return;
  }
  using Real = typename IsDataValueType<M>::ElementaryType;
  const Long n_points = sp.n_points;
  const Int multiplicity = sp.multiplicity;
  const Coordinate total_site = psel.total_site;
  Qassert(n_points == psel.size());
  Qassert(sp.points.size() == n_points * multiplicity);
  qthread_for(idx, n_points, {
    const Coordinate xg = psel[idx];
    const Long gindex = index_from_coordinate(xg, total_site);
    RngState rsi = rs.newtype(gindex);
    Vector<M> v = sp.get_elems(idx);
    Vector<Real> dv((Real*)v.data(), v.data_size() / sizeof(Real));
    for (Int m = 0; m < dv.size(); ++m) {
      dv[m] = g_rand_gen(rsi, center, sigma);
    }
  });
}

// -------------------------------------------

template <class M>
void lat_data_from_selected_points(LatData& ld, const SelectedPoints<M>& sp)
{
  TIMER("lat_data_from_selected_points(sp)");
  Qassert(is_composed_of_real_d<M>());
  ld.init();
  ld.info.push_back(lat_dim_number("idx", 0, sp.n_points - 1));
  ld.info.push_back(lat_dim_number("m", 0, sp.multiplicity - 1));
  if (is_composed_of_complex_d<M>()) {
    const Long n_v = sizeof(M) / sizeof(ComplexD);
    Qassert(n_v * (Long)sizeof(ComplexD) == (Long)sizeof(M));
    ld.info.push_back(lat_dim_number("v", 0, n_v - 1));
    ld.info.push_back(lat_dim_re_im());
  } else if (is_composed_of_real_d<M>()) {
    const Long n_v = sizeof(M) / sizeof(RealD);
    Qassert(n_v * (Long)sizeof(RealD) == (Long)sizeof(M));
    ld.info.push_back(lat_dim_number("v", 0, n_v - 1));
  } else {
    qerr(fname + ssprintf(": get_type_name(M)=%s", get_type_name<M>().c_str()));
  }
  lat_data_alloc(ld);
  assign(get_data(ld.res), get_data(sp.points));
}

template <class M>
void selected_points_from_lat_data(SelectedPoints<M>& sp, const LatData& ld)
{
  TIMER("selected_points_from_lat_data(sp,ld)");
  Qassert(is_composed_of_real_d<M>());
  if (is_composed_of_complex_d<M>()) {
    Qassert(ld.info.size() == 4);
  } else if (is_composed_of_real_d<M>()) {
    Qassert(ld.info.size() == 3);
  } else {
    qerr(fname + ssprintf(": get_type_name(M)=%s", get_type_name<M>().c_str()));
  }
  Qassert(ld.info[0].name == "idx");
  Qassert(ld.info[1].name == "m");
  Qassert(ld.info[2].name == "v");
  const Long n_points = ld.info[0].size;
  const Long multiplicity = ld.info[1].size;
  const Long sizof_M_vs_sizeof_v = ld.info[2].size;
  if (is_composed_of_complex_d<M>()) {
    Qassert(sizeof(M) == sizof_M_vs_sizeof_v * sizeof(ComplexD));
    Qassert(ld.info[3].name == "re-im");
    Qassert(ld.info[3].size == 2);
  } else {
    Qassert(sizeof(M) == sizof_M_vs_sizeof_v * sizeof(RealD));
  }
  sp.init(n_points, multiplicity, PointsDistType::Global);
  assign(get_data(sp.points), get_data(ld.res));
}

template <class M>
void lat_data_from_selected_points(LatDataRealF& ld,
                                   const SelectedPoints<M>& sp)
{
  TIMER("lat_data_from_selected_points(sp)");
  Qassert(is_composed_of_real_f<M>());
  ld.init();
  ld.info.push_back(lat_dim_number("idx", 0, sp.n_points - 1));
  ld.info.push_back(lat_dim_number("m", 0, sp.multiplicity - 1));
  if (is_composed_of_complex_f<M>()) {
    const Long n_v = sizeof(M) / sizeof(ComplexF);
    Qassert(n_v * (Long)sizeof(ComplexF) == (Long)sizeof(M));
    ld.info.push_back(lat_dim_number("v", 0, n_v - 1));
    ld.info.push_back(lat_dim_re_im());
  } else if (is_composed_of_real_f<M>()) {
    const Long n_v = sizeof(M) / sizeof(RealF);
    Qassert(n_v * (Long)sizeof(RealF) == (Long)sizeof(M));
    ld.info.push_back(lat_dim_number("v", 0, n_v - 1));
  } else {
    qerr(fname + ssprintf(": get_type_name(M)=%s", get_type_name<M>().c_str()));
  }
  lat_data_alloc(ld);
  assign(get_data(ld.res), get_data(sp.points));
}

template <class M>
void selected_points_from_lat_data(SelectedPoints<M>& sp,
                                   const LatDataRealF& ld)
{
  TIMER("selected_points_from_lat_data(sp,ld)");
  Qassert(is_composed_of_real_f<M>());
  if (is_composed_of_complex_f<M>()) {
    Qassert(ld.info.size() == 4);
  } else if (is_composed_of_real_f<M>()) {
    Qassert(ld.info.size() == 3);
  } else {
    qerr(fname + ssprintf(": get_type_name(M)=%s", get_type_name<M>().c_str()));
  }
  Qassert(ld.info[0].name == "idx");
  Qassert(ld.info[1].name == "m");
  Qassert(ld.info[2].name == "v");
  const Long n_points = ld.info[0].size;
  const Long multiplicity = ld.info[1].size;
  const Long sizof_M_vs_sizeof_v = ld.info[2].size;
  if (is_composed_of_complex_f<M>()) {
    Qassert(sizeof(M) == sizof_M_vs_sizeof_v * sizeof(ComplexF));
    Qassert(ld.info[3].name == "re-im");
    Qassert(ld.info[3].size == 2);
  } else {
    Qassert(sizeof(M) == sizof_M_vs_sizeof_v * sizeof(RealF));
  }
  sp.init(n_points, multiplicity, PointsDistType::Global);
  assign(get_data(sp.points), get_data(ld.res));
}

template <class M>
void lat_data_from_selected_points(LatDataLong& ld, const SelectedPoints<M>& sp)
{
  TIMER("lat_data_from_selected_points(sp)");
  Qassert(is_composed_of_long<M>());
  ld.init();
  ld.info.push_back(lat_dim_number("idx", 0, sp.n_points - 1));
  ld.info.push_back(lat_dim_number("m", 0, sp.multiplicity - 1));
  const Long n_v = sizeof(M) / sizeof(Long);
  Qassert(n_v * (Long)sizeof(Long) == (Long)sizeof(M));
  ld.info.push_back(lat_dim_number("v", 0, n_v - 1));
  lat_data_alloc(ld);
  assign(get_data(ld.res), get_data(sp.points));
}

template <class M>
void selected_points_from_lat_data(SelectedPoints<M>& sp, const LatDataLong& ld)
{
  TIMER("selected_points_from_lat_data(sp,ld)");
  Qassert(is_composed_of_long<M>());
  Qassert(ld.info.size() == 3);
  Qassert(ld.info[0].name == "idx");
  Qassert(ld.info[1].name == "m");
  Qassert(ld.info[2].name == "v");
  const Long n_points = ld.info[0].size;
  const Long multiplicity = ld.info[1].size;
  const Long sizof_M_vs_sizeof_v = ld.info[2].size;
  Qassert(sizeof(M) == sizof_M_vs_sizeof_v * sizeof(Long));
  sp.init(n_points, multiplicity, PointsDistType::Global);
  assign(get_data(sp.points), get_data(ld.res));
}

template <class M>
void lat_data_from_selected_points(LatDataInt& ld, const SelectedPoints<M>& sp)
{
  TIMER("lat_data_from_selected_points(sp)");
  Qassert(is_composed_of_int<M>());
  ld.init();
  ld.info.push_back(lat_dim_number("idx", 0, sp.n_points - 1));
  ld.info.push_back(lat_dim_number("m", 0, sp.multiplicity - 1));
  const Long n_v = sizeof(M) / sizeof(Int);
  Qassert(n_v * (Long)sizeof(Int) == (Long)sizeof(M));
  ld.info.push_back(lat_dim_number("v", 0, n_v - 1));
  lat_data_alloc(ld);
  assign(get_data(ld.res), get_data(sp.points));
}

template <class M>
void selected_points_from_lat_data(SelectedPoints<M>& sp, const LatDataInt& ld)
{
  TIMER("selected_points_from_lat_data(sp,ld)");
  Qassert(is_composed_of_int<M>());
  Qassert(ld.info.size() == 3);
  Qassert(ld.info[0].name == "idx");
  Qassert(ld.info[1].name == "m");
  Qassert(ld.info[2].name == "v");
  const Long n_points = ld.info[0].size;
  const Long multiplicity = ld.info[1].size;
  const Long sizof_M_vs_sizeof_v = ld.info[2].size;
  Qassert(sizeof(M) == sizof_M_vs_sizeof_v * sizeof(Int));
  sp.init(n_points, multiplicity, PointsDistType::Global);
  assign(get_data(sp.points), get_data(ld.res));
}

// -------------------------------------------

template <class M>
void save_selected_points(const SelectedPoints<M>& sp, QFile& qfile)
{
  TIMER_VERBOSE("save_selected_points(sp,qfile)");
  Qassert(sp.points_dist_type == PointsDistType::Global);
  if (get_id_node() == 0) {
    if (is_composed_of_real_d<M>()) {
      LatData ld;
      lat_data_from_selected_points(ld, sp);
      ld.save(qfile);
    } else if (is_composed_of_real_f<M>()) {
      LatDataRealF ld;
      lat_data_from_selected_points(ld, sp);
      ld.save(qfile);
    } else if (is_composed_of_long<M>()) {
      LatDataLong ld;
      lat_data_from_selected_points(ld, sp);
      ld.save(qfile);
    } else if (is_composed_of_int<M>()) {
      LatDataInt ld;
      lat_data_from_selected_points(ld, sp);
      ld.save(qfile);
    } else {
      Qassert(false);
    }
  }
}

template <class M>
void load_selected_points(SelectedPoints<M>& sp, QFile& qfile)
{
  TIMER_VERBOSE("load_selected_points(sp,qfile)");
  Long n_points = 0;
  Long multiplicity = 0;
  if (get_id_node() == 0) {
    if (is_composed_of_real_d<M>()) {
      LatData ld;
      ld.load(qfile);
      selected_points_from_lat_data(sp, ld);
    } else if (is_composed_of_real_f<M>()) {
      LatDataRealF ld;
      ld.load(qfile);
      selected_points_from_lat_data(sp, ld);
    } else if (is_composed_of_long<M>()) {
      LatDataLong ld;
      ld.load(qfile);
      selected_points_from_lat_data(sp, ld);
    } else if (is_composed_of_int<M>()) {
      LatDataInt ld;
      ld.load(qfile);
      selected_points_from_lat_data(sp, ld);
    } else {
      Qassert(false);
    }
    n_points = sp.n_points;
    multiplicity = sp.multiplicity;
  }
  bcast(get_data_one_elem(n_points));
  bcast(get_data_one_elem(multiplicity));
  if (get_id_node() != 0) {
    sp.init(n_points, multiplicity, PointsDistType::Global);
  }
  vector<M> buffer(sp.points.size());
  assign(get_data(buffer), get_data(sp.points));
  bcast(get_data(buffer));
  assign(get_data(sp.points), get_data(buffer));
}

template <class M>
void save_selected_points(const SelectedPoints<M>& sp, const std::string& fn)
{
  QFile qfile;
  if (get_id_node() == 0) {
    qfile = qfopen(fn + ".partial", "w");
  }
  save_selected_points(sp, qfile);
  qfclose(qfile);
  qrename_info(fn + ".partial", fn);
}

template <class M>
void load_selected_points(SelectedPoints<M>& sp, const std::string& fn)
{
  QFile qfile;
  if (get_id_node() == 0) {
    qfile = qfopen(fn, "r");
  }
  load_selected_points(sp, qfile);
  qfclose(qfile);
}

template <class M>
std::string save_selected_points_str(const SelectedPoints<M>& sp)
// Right now only return content at node 0.
{
  QFile qfile;
  if (get_id_node() == 0) {
    qfile =
        qfopen(QFileType::String, "/ save SelectedPoints /", QFileMode::Write);
  }
  save_selected_points(sp, qfile);
  std::string ret;
  if (get_id_node() == 0) {
    ret = qfile.content();
  }
  qfclose(qfile);
  return ret;
}

template <class M>
void load_selected_points_str(SelectedPoints<M>& sp, std::string& content)
// Allow to destroy `content` to be more efficient.
// Only need to set the content at node 0.
{
  QFile qfile;
  if (get_id_node() == 0) {
    qfile = qfopen(QFileType::String, "/ load SelectedPoints /",
                   QFileMode::Read, content);
  }
  load_selected_points(sp, qfile);
  qfclose(qfile);
}

// --------------------

#ifdef QLAT_INSTANTIATE_SELECTED_POINTS
#define QLAT_EXTERN
#else
#define QLAT_EXTERN extern
#endif

#define QLAT_EXTERN_TEMPLATE(TYPENAME)                                    \
                                                                          \
  QLAT_EXTERN template SelectedPoints<TYPENAME>& operator+= <TYPENAME>(   \
      SelectedPoints<TYPENAME> & f, const SelectedPoints<TYPENAME>& f1);  \
                                                                          \
  QLAT_EXTERN template SelectedPoints<TYPENAME>& operator-= <TYPENAME>(   \
      SelectedPoints<TYPENAME> & f, const SelectedPoints<TYPENAME>& f1);  \
                                                                          \
  QLAT_EXTERN template SelectedPoints<TYPENAME>& operator*=               \
      <TYPENAME>(SelectedPoints<TYPENAME> & f, const double factor);      \
                                                                          \
  QLAT_EXTERN template SelectedPoints<TYPENAME>& operator*=               \
      <TYPENAME>(SelectedPoints<TYPENAME> & f, const ComplexD factor);    \
                                                                          \
  QLAT_EXTERN template void only_keep_selected_points<TYPENAME>(          \
      Field<TYPENAME> & f, const PointsSelection& psel);                  \
                                                                          \
  QLAT_EXTERN template double qnorm<TYPENAME>(                            \
      const SelectedPoints<TYPENAME>& sp);                                \
                                                                          \
  QLAT_EXTERN template void qnorm_field<TYPENAME>(                        \
      SelectedPoints<double> & sp, const SelectedPoints<TYPENAME>& sp1);  \
                                                                          \
  QLAT_EXTERN template void set_selected_points<TYPENAME>(                \
      SelectedPoints<TYPENAME> & sp, const Field<TYPENAME>& f,            \
      const PointsSelection& psel);                                       \
                                                                          \
  QLAT_EXTERN template void set_selected_points<TYPENAME>(                \
      SelectedPoints<TYPENAME> & sp, const Field<TYPENAME>& f,            \
      const PointsSelection& psel, const Int m);                          \
                                                                          \
  QLAT_EXTERN template void set_selected_points<TYPENAME>(                \
      SelectedPoints<TYPENAME> & sp, const SelectedPoints<TYPENAME>& sp0, \
      const PointsSelection& psel, const PointsSelection& psel0,          \
      const bool is_keeping_data);                                        \
                                                                          \
  QLAT_EXTERN template void set_field_selected<TYPENAME>(                 \
      Field<TYPENAME> & f, const SelectedPoints<TYPENAME>& sp,            \
      const Geometry& geo_, const PointsSelection& psel,                  \
      const bool is_keeping_data);                                        \
                                                                          \
  QLAT_EXTERN template void set_field_selected<TYPENAME>(                 \
      Field<TYPENAME> & f, const SelectedPoints<TYPENAME>& sp,            \
      const Geometry& geo_, const PointsSelection& psel,                  \
      const Int multiplicity, const Int m, const bool is_keeping_data);   \
                                                                          \
  QLAT_EXTERN template void acc_field<TYPENAME>(                          \
      Field<TYPENAME> & f, const SelectedPoints<TYPENAME>& sp,            \
      const Geometry& geo_, const PointsSelection& psel);                 \
                                                                          \
  QLAT_EXTERN template void field_glb_sum_tslice<TYPENAME>(               \
      SelectedPoints<TYPENAME> & sp, const Field<TYPENAME>& f,            \
      const Int t_dir);                                                   \
                                                                          \
  QLAT_EXTERN template void lat_data_from_selected_points<TYPENAME>(      \
      LatData & ld, const SelectedPoints<TYPENAME>& sp);                  \
                                                                          \
  QLAT_EXTERN template void selected_points_from_lat_data<TYPENAME>(      \
      SelectedPoints<TYPENAME> & sp, const LatData& ld);                  \
                                                                          \
  QLAT_EXTERN template void lat_data_from_selected_points<TYPENAME>(      \
      LatDataRealF & ld, const SelectedPoints<TYPENAME>& sp);             \
                                                                          \
  QLAT_EXTERN template void selected_points_from_lat_data<TYPENAME>(      \
      SelectedPoints<TYPENAME> & sp, const LatDataRealF& ld);             \
                                                                          \
  QLAT_EXTERN template void lat_data_from_selected_points<TYPENAME>(      \
      LatDataLong & ld, const SelectedPoints<TYPENAME>& sp);              \
                                                                          \
  QLAT_EXTERN template void selected_points_from_lat_data<TYPENAME>(      \
      SelectedPoints<TYPENAME> & sp, const LatDataLong& ld);              \
                                                                          \
  QLAT_EXTERN template void lat_data_from_selected_points<TYPENAME>(      \
      LatDataInt & ld, const SelectedPoints<TYPENAME>& sp);               \
                                                                          \
  QLAT_EXTERN template void selected_points_from_lat_data<TYPENAME>(      \
      SelectedPoints<TYPENAME> & sp, const LatDataInt& ld);               \
                                                                          \
  QLAT_EXTERN template void save_selected_points<TYPENAME>(               \
      const SelectedPoints<TYPENAME>& sp, QFile& qfile);                  \
                                                                          \
  QLAT_EXTERN template void load_selected_points<TYPENAME>(               \
      SelectedPoints<TYPENAME> & sp, QFile & qfile);                      \
                                                                          \
  QLAT_EXTERN template void save_selected_points<TYPENAME>(               \
      const SelectedPoints<TYPENAME>& sp, const std::string& fn);         \
                                                                          \
  QLAT_EXTERN template void load_selected_points<TYPENAME>(               \
      SelectedPoints<TYPENAME> & sp, const std::string& fn);              \
                                                                          \
  QLAT_EXTERN template std::string save_selected_points_str<TYPENAME>(    \
      const SelectedPoints<TYPENAME>& sp);                                \
                                                                          \
  QLAT_EXTERN template void load_selected_points_str<TYPENAME>(           \
      SelectedPoints<TYPENAME> & sp, std::string & content);              \
                                                                          \
  QLAT_EXTERN template void set_u_rand<TYPENAME>(                         \
      SelectedPoints<TYPENAME> & sp, const PointsSelection& psel,         \
      const RngState& rs, const RealD upper, const RealD lower);          \
                                                                          \
  QLAT_EXTERN template void set_g_rand<TYPENAME>(                         \
      SelectedPoints<TYPENAME> & sp, const PointsSelection& psel,         \
      const RngState& rs, const RealD center, const RealD sigma)

QLAT_CALL_WITH_TYPES(QLAT_EXTERN_TEMPLATE);
#undef QLAT_EXTERN_TEMPLATE

#undef QLAT_EXTERN

}  // namespace qlat
