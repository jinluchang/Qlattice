#pragma once

#include <qlat/field.h>

namespace qlat
{  //

typedef std::vector<Coordinate> PointSelection;

inline PointSelection mk_tslice_point_selection(const int t_size,
                                                const int t_dir = 3)
{
  PointSelection psel;
  psel.resize(t_size);
  qassert(0 <= t_dir and t_dir < 4);
  qthread_for(idx, t_size, {
    psel[idx] = Coordinate();
    psel[idx][t_dir] = idx;
  });
  return psel;
}

inline PointSelection mk_tslice_point_selection(const Coordinate& total_site,
                                                const int t_dir = 3)
{
  return mk_tslice_point_selection(total_site[t_dir], t_dir);
}

inline PointSelection mk_random_point_selection(const Coordinate& total_site,
                                                const long num,
                                                const RngState& rs,
                                                const long pool_factor = 2)
// same rs for all node for uniform result
{
  TIMER_VERBOSE("mk_random_point_selection");
  if (num == 0) {
    PointSelection psel;
    return psel;
  }
  qassert(num > 0);
  PointSelection psel_pool(pool_factor * num);
#pragma omp parallel for
  for (long i = 0; i < (long)psel_pool.size(); ++i) {
    RngState rsi = rs.split(i);
    Coordinate xg;
    for (int m = 0; m < 4; ++m) {
      xg[m] = modl(rand_gen(rsi), total_site[m]);
    }
    psel_pool[i] = xg;
  }
  PointSelection psel(num, Coordinate(-1, -1, -1, -1));
  long idx = 0;
  for (long i = 0; i < (long)psel.size(); ++i) {
    while (idx < (long)psel_pool.size()) {
      const Coordinate xg = psel_pool[idx];
      idx += 1;
      bool is_repeat = false;
      for (long j = 0; j < i; ++j) {
        if (xg == psel[j]) {
          is_repeat = true;
          break;
        }
      }
      if (not is_repeat) {
        psel[i] = xg;
        break;
      }
    }
  }
  if (psel.back() != Coordinate(-1, -1, -1, -1)) {
    return psel;
  } else {
    displayln_info(
        fname +
        ssprintf(": pool_factor=%d is too small, rerun with larger factor.",
                 pool_factor));
    return mk_random_point_selection(total_site, num, rs, pool_factor + 2);
  }
}

inline void save_point_selection(const PointSelection& psel,
                                 const std::string& path)
{
  TIMER_VERBOSE("save_point_selection");
  QFile qfile = qfopen(path + ".partial", "w");
  qfprintf(qfile, "%ld\n", (long)psel.size());
  for (long i = 0; i < (long)psel.size(); ++i) {
    const Coordinate& c = psel[i];
    qfprintf(qfile, "%5ld    %3d %3d %3d %3d\n", i, c[0], c[1], c[2], c[3]);
  }
  qfclose(qfile);
  qrename(path + ".partial", path);
}

inline void save_point_selection_info(const PointSelection& psel,
                                      const std::string& path)
{
  TIMER_VERBOSE("save_point_selection_info");
  if (0 == get_id_node()) {
    save_point_selection(psel, path);
  }
}

inline PointSelection load_point_selection(const std::string& path)
{
  TIMER_VERBOSE("load_point_selection");
  const std::vector<std::string> lines = qgetlines(path);
  qassert(lines.size() > 0);
  const long len = read_long(lines[0]);
  qassert(len + 1 <= (long)lines.size());
  PointSelection psel;
  for (long k = 1; k < len + 1; ++k) {
    const std::vector<std::string> strs = split_line_with_spaces(lines[k]);
    if (strs.size() >= 5) {
      qassert(k - 1 == read_long(strs[0]));
      const Coordinate xg(read_long(strs[1]), read_long(strs[2]),
                          read_long(strs[3]), read_long(strs[4]));
      psel.push_back(xg);
    } else {
      displayln(fname + ssprintf(": line is '%s'.", lines[k].c_str()));
      qassert(false);
    }
  }
  return psel;
}

inline PointSelection load_point_selection_info(const std::string& path)
{
  TIMER_VERBOSE("load_point_selection_info");
  PointSelection psel;
  if (0 == get_id_node()) {
    psel = load_point_selection(path);
  }
  bcast(psel);
  return psel;
}

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
  return qcast((SelectedPoints<N>&)x);
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
Vector<M> get_data(const SelectedPoints<M>& sp)
{
  return get_data(sp.points);
}

template <class M>
void set_zero(SelectedPoints<M>& sp)
{
  TIMER("set_zero(SelectedPoints)");
  set_zero(get_data(sp));
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
  qacc_for(idx, n_points, {
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

}  // namespace qlat
