#pragma once

#include <qlat/config.h>
#include <qlat/coordinate-d.h>
#include <qlat/field.h>

#include <map>

QLAT_START_NAMESPACE

inline void set_g_rand_anti_hermitian_matrix_field(Field<ColorMatrix>& fc, const RngState& rs, const double sigma)
{
  TIMER("set_g_rand_anti_hermitian_matrix_field");
  const Geometry& geo = fc.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate& xl = geo.coordinate_from_index(index);
    const Coordinate& xg = geo.coordinate_g_from_l(xl);
    const long gindex = geo.g_index_from_g_coordinate(xg);
    RngState rsi(rs, gindex);
    Vector<ColorMatrix> v = fc.get_elems(xl);
    for (int m = 0; m < v.size(); ++m) {
      v[m] = make_g_rand_anti_hermitian_matrix(rsi, sigma);
    }
  }
}

inline void set_g_rand_color_matrix_field(Field<ColorMatrix>& fc, const RngState& rs, const double sigma)
{
  TIMER("set_g_rand_color_matrix_field");
  const Geometry& geo = fc.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate& xl = geo.coordinate_from_index(index);
    const Coordinate& xg = geo.coordinate_g_from_l(xl);
    const long gindex = geo.g_index_from_g_coordinate(xg);
    RngState rsi(rs, gindex);
    Vector<ColorMatrix> v = fc.get_elems(xl);
    for (int m = 0; m < v.size(); ++m) {
      v[m] = make_color_matrix_exp(make_g_rand_anti_hermitian_matrix(rsi, sigma));
    }
  }
}

inline ColorMatrix gf_wilson_line_no_comm(const GaugeField& gf, const Coordinate& xl, const std::vector<int>& path)
{
  ColorMatrix ret;
  set_unit(ret);
  Coordinate xl1 = xl;
  for (int i = 0; i < path.size(); ++i) {
    const int dir = path[i];
    qassert(-DIM <= dir && dir < DIM);
    if (0 <= dir) {
      ret *= gf.get_elem(xl1,dir);
      xl1[dir] += 1;
    } else {
      xl1[-dir-1] -= 1;
      ret *= matrix_adjoint(gf.get_elem(xl1,-dir-1));
    }
  }
  return ret;
}

inline ColorMatrix gf_staple_no_comm_v1(const GaugeField& gf, const Coordinate& xl, const int mu)
{
  ColorMatrix ret;
  set_zero(ret);
  const Coordinate xl_mu = coordinate_shifts(xl,mu);
  for (int m = 0; m < DIM; ++m) {
    if (mu != m) {
      ret += gf.get_elem(xl, m) *
        gf.get_elem(coordinate_shifts(xl,m), mu) *
        matrix_adjoint(gf.get_elem(xl_mu, m));
      ret += matrix_adjoint(gf.get_elem(coordinate_shifts(xl,-m-1), m)) *
        gf.get_elem(coordinate_shifts(xl,-m-1), mu) *
        gf.get_elem(coordinate_shifts(xl_mu,-m-1), m);
    }
  }
  return ret;
}

inline ColorMatrix gf_staple_no_comm_v2(const GaugeField& gf, const Coordinate& xl, const int mu)
{
  ColorMatrix ret;
  set_zero(ret);
  std::vector<int> path(3);
  path[1] = mu;
  for (int m = 0; m < DIM; ++m) {
    if (mu != m) {
      path[0] = m;
      path[2] = -m-1;
      ret += gf_wilson_line_no_comm(gf, xl, path);
    }
  }
  for (int m = 0; m < DIM; ++m) {
    if (mu != m) {
      path[0] = -m-1;
      path[2] = m;
      ret += gf_wilson_line_no_comm(gf, xl, path);
    }
  }
  return ret;
}

inline ColorMatrix gf_staple_no_comm(const GaugeField& gf, const Coordinate& xl, const int mu)
{
  return gf_staple_no_comm_v1(gf, xl, mu);
  // return gf_staple_no_comm_v2(gf, xl, mu);
}

inline void set_multiply_wilson_line_field_partial_comm(FieldM<ColorMatrix,1>& wlf, FieldM<ColorMatrix,1>& wlf1, const GaugeField& gf1, const std::vector<int>& path)
  // gf1 need to be refresh_expanded.
  // wlf1 need to have correct size
  // wlf1 will be modified
  // wlf will be initialized
{
  TIMER("set_multiply_wilson_line_field_partial_comm");
  const Geometry geo = geo_reform(gf1.geo);
  const Geometry geo1 = geo_reform(gf1.geo, 1, 1);
  qassert(&wlf != &wlf1);
  wlf.init(geo);
  for (size_t i = 0; i < path.size(); ++i) {
    const int dir = path[i];
    qassert(-DIM <= dir && dir < DIM);
    refresh_expanded(wlf1);
#pragma omp parallel for
    for (long index = 0; index < geo.local_volume(); ++index) {
      Coordinate xl = geo.coordinate_from_index(index);
      ColorMatrix& l1 = wlf.get_elem(xl);
      if (0 <= dir) {
        xl[dir] -= 1;
        const ColorMatrix& link = gf1.get_elem(xl, dir);
        const ColorMatrix& l0 = wlf1.get_elem(xl);
        l1 = l0 * link;
      } else {
        const ColorMatrix& link = gf1.get_elem(xl, -dir-1);
        xl[-dir-1] += 1;
        const ColorMatrix& l0 = wlf1.get_elem(xl);
        l1 = l0 * matrix_adjoint(link);
      }
    }
    if (i != path.size() - 1) {
      wlf1 = wlf;
    }
  }
}

enum WilsonLinePathType
{
  WilsonLinePathType_simp,
  WilsonLinePathType_comp_seq,
  WilsonLinePathType_comp_par,
  WilsonLinePathType_comp_par_multi_target,
};

struct WilsonLinePath
{
  WilsonLinePathType type;
  std::vector<int> simp;
  std::vector<WilsonLinePath> comp;
  Coordinate target;
  //
  void init()
  {
    type = WilsonLinePathType_simp;
  }
  //
  WilsonLinePath()
  {
    init();
  }
};

struct WilsonLinePathStop
{
  Coordinate x;
  std::vector<std::vector<int> > paths;
  int num_origins;
  //
  WilsonLinePathStop()
  {
    num_origins = 0;
  }
};

struct WilsonLinePathSegment
{
  Coordinate target;
  std::map<Coordinate,WilsonLinePathStop> stops;
};

inline double coordinate_distance_from_wilson_line(const Coordinate& x, const Coordinate& target_wilson_line)
{
  const CoordinateD target(target_wilson_line);
  const double target_len = coordinate_len(target);
  const CoordinateD vdir = target / target_len;
  const CoordinateD xd(x);
  const CoordinateD dis = xd - dot_product(xd, vdir) * vdir;
  return coordinate_len(dis);
}

inline std::vector<int> find_next_dirs(const Coordinate& loc, const Coordinate& target_wilson_line)
{
  std::vector<int> dirs;
  for (int i = 0; i < DIM; ++i) {
    const int x = target_wilson_line[i];
    if (0 < x) {
      dirs.push_back(i);
    } else if (0 > x) {
      dirs.push_back(-i-1);
    }
  }
  qassert(0 < dirs.size());
  if (1 == dirs.size()) {
    return dirs;
  } else {
    const double eps_dis = 1.0e-8;
    std::vector<int> next_dirs;
    const double dis = coordinate_distance_from_wilson_line(loc, target_wilson_line);
    if (dis < eps_dis) {
      next_dirs = dirs;
    } else {
      double min_dis = 1e6; // a very large number
      for (int i = 0; i < dirs.size(); ++i) {
        const double next_dis = coordinate_distance_from_wilson_line(coordinate_shifts(loc, dirs[i]), target_wilson_line);
        if (next_dis < min_dis) {
          min_dis = next_dis;
        }
      }
      for (int i = 0; i < dirs.size(); ++i) {
        const double next_dis = coordinate_distance_from_wilson_line(coordinate_shifts(loc, dirs[i]), target_wilson_line);
        if (next_dis <= min_dis + eps_dis) {
          next_dirs.push_back(dirs[i]);
        }
      }
    }
    return next_dirs;
  }
}

inline WilsonLinePathSegment make_wilson_line_path(const Coordinate& target)
{
  WilsonLinePathSegment ret;
  ret.target = target;
  std::vector<Coordinate> cs;
  cs.push_back(Coordinate());
  for (int i = 0; i < cs.size(); ++i) {
    const Coordinate& c = cs[i];
    std::vector<int> dirs = find_next_dirs(c, target);
    WilsonLinePathStop& ps = ret.stops[c];
    ps.x = c;
    ps.paths.resize(dirs.size());
    for (int k = 0; k < dirs.size(); ++k) {
      int dir = dirs[k];
      ps.paths[k].push_back(dir);
      Coordinate nc = coordinate_shifts(c, dir);
      std::vector<int> ndirs = find_next_dirs(nc, target);
      while (ndirs.size() == 1 && ret.stops.find(nc) == ret.stops.end()) {
        dir = ndirs[0];
        ps.paths[k].push_back(dir);
        nc = coordinate_shifts(nc, dir);
        ndirs = find_next_dirs(nc, target);
      }
      qassert(nc == coordinate_shifts(c, ps.paths[k]));
      cs.push_back(nc);
    }
  }
  return ret;
}

inline void acc_wilson_line_path_segment(WilsonLinePathSegment& path)
{
  std::map<Coordinate,WilsonLinePathStop>& stops = path.stops;
  std::vector<Coordinate> cs;
  cs.push_back(Coordinate());
  for (int i = 0; i < cs.size(); ++i) {
    const Coordinate& c = cs[i];
    const WilsonLinePathStop& ps = path.stops[c];
    qassert(c == ps.x);
    for (int k = 0; k < ps.paths.size(); ++k) {
      const Coordinate nc = coordinate_shifts(c, ps.paths[k]);
      path.stops[nc].num_origins = 0;
    }
  }
  cs.clear();
  cs.push_back(Coordinate());
  for (int i = 0; i < cs.size(); ++i) {
    const Coordinate& c = cs[i];
    const WilsonLinePathStop& ps = path.stops[c];
    qassert(c == ps.x);
    for (int k = 0; k < ps.paths.size(); ++k) {
      const Coordinate nc = coordinate_shifts(c, ps.paths[k]);
      path.stops[nc].num_origins += 1;
    }
  }
}

inline void set_multiply_wilson_line_field_partial_comm(FieldM<ColorMatrix,1>& wlf, FieldM<ColorMatrix,1>& wlf1, const GaugeField& gf1, const WilsonLinePathSegment& path)
{
  TIMER("set_multiply_wilson_line_field_partial_comm");
  const Geometry geo = geo_reform(gf1.geo);
  WilsonLinePathSegment pacc = path;
  acc_wilson_line_path_segment(pacc);
  std::vector<Coordinate> cs;
  std::vector<FieldM<ColorMatrix,1> > fs;
  std::map<Coordinate,int> dict;
  cs.push_back(Coordinate());
  fs.push_back(FieldM<ColorMatrix,1>());
  fs.back().init(geo);
  fs.back() = wlf1;
  dict[cs[0]] = 0;
  FieldM<ColorMatrix,1> wlf_t;
  wlf_t.init(geo);
  wlf.init(geo);
  while (true) {
    for (int i = 0; i < cs.size(); ++i) {
      const Coordinate& c = cs[i];
      const WilsonLinePathStop& ps = pacc.stops[c];
      qassert(c == ps.x);
      if (is_initialized(fs[i]) && ps.num_origins == 0) {
        if (c == pacc.target) {
          qassert(ps.paths.size() == 0);
          wlf = fs[i];
          return;
        }
        wlf1 = fs[i];
        for (int k = 0; k < ps.paths.size(); ++k) {
          qassert(cs.size() == fs.size());
          const Coordinate nc = coordinate_shifts(c, ps.paths[k]);
          set_multiply_wilson_line_field_partial_comm(wlf, wlf1, gf1, ps.paths[k]);
          Handle<FieldM<ColorMatrix,1> > hf;
          if (dict.find(nc) == dict.end()) {
            cs.push_back(nc);
            fs.push_back(FieldM<ColorMatrix,1>());
            fs.back().init(geo);
            set_zero(fs.back());
            hf.init(fs.back());
          } else {
            hf.init(fs[dict[nc]]);
          }
          hf() += wlf_t;
          pacc.stops[nc].num_origins -= 1;
        }
        fs[i].init();
      }
    }
  }
}

inline void set_multiply_wilson_line_field_partial_comm(FieldM<ColorMatrix,1>& wlf, FieldM<ColorMatrix,1>& wlf1, const GaugeField& gf1, const WilsonLinePath& path)
{
  TIMER("set_multiply_wilson_line_field_partial_comm");
  if (WilsonLinePathType_simp == path.type) {
    qassert(0 == path.comp.size());
    set_multiply_wilson_line_field_partial_comm(wlf, wlf1, gf1, path.simp);
  } else if (WilsonLinePathType_comp_seq == path.type) {
    qassert(0 == path.simp.size());
    for (int i = 0; i < path.comp.size(); ++i) {
      set_multiply_wilson_line_field_partial_comm(wlf, wlf1, gf1, path.comp[i]);
      wlf1 = wlf;
    }
  } else if (WilsonLinePathType_comp_par == path.type) {
    qassert(0 == path.simp.size());
    const Geometry geo = geo_reform(wlf1.geo);
    FieldM<ColorMatrix,1> wlf_init, wlf_t;
    wlf_init.init(geo);
    wlf_t.init(geo);
    wlf_init = wlf1;
    set_zero(wlf);
    for (int i = 0; i < path.comp.size(); ++i) {
      wlf1 = wlf_init;
      set_multiply_wilson_line_field_partial_comm(wlf_t, wlf1, gf1, path.comp[i]);
      wlf += wlf_t;
    }
  } else {
    // cannot be WilsonLinePathType_comp_par_multi_target
    qassert(false);
  }
}

inline ColorMatrix gf_avg_wilson_line(const GaugeField& gf, const WilsonLinePath& path)
{
  TIMER("gf_avg_wilson_line");
  const Geometry geo = geo_reform(gf.geo);
  const Coordinate expansion_left(1, 1, 1, 1);
  const Coordinate expansion_right(0, 0, 0, 0);
  GaugeField gf1;
  gf1.init(geo_resize(geo, expansion_left, expansion_right));
  gf1 = gf;
  refresh_expanded(gf1);
  FieldM<ColorMatrix,1> wlf, wlf1;
  wlf1.init(geo_resize(geo, 1));
  wlf.init(geo);
  set_unit(wlf1);
  set_multiply_wilson_line_field_partial_comm(wlf, wlf1, gf1, path);
  return field_glb_sum_double(wlf)[0] / geo.total_volume();
}

inline ColorMatrix gf_avg_wilson_line(const GaugeField& gf, const std::vector<int>& spath)
{
  TIMER("gf_avg_wilson_line");
  WilsonLinePath path;
  path.type = WilsonLinePathType_simp;
  path.simp = spath;
  return gf_avg_wilson_line(gf, path);
}

inline std::vector<int> make_wilson_line_path(const int l1, const int l2, const int dir1, const int dir2)
{
  std::vector<int> path;
  for (int i = 0; i < l1; ++i) {
    path.push_back(dir1);
  }
  for (int i = 0; i < l2; ++i) {
    path.push_back(dir2);
  }
  for (int i = 0; i < l1; ++i) {
    path.push_back(-dir1-1);
  }
  for (int i = 0; i < l2; ++i) {
    path.push_back(-dir2-1);
  }
  return path;
}

inline ColorMatrix gf_avg_wilson_loop(const GaugeField& gf, const int l, const int t)
{
  TIMER("gf_avg_wilson_loop");
  ColorMatrix m =
    gf_avg_wilson_line(gf, make_wilson_line_path(l, t, 0, 3)) +
    gf_avg_wilson_line(gf, make_wilson_line_path(l, t, 1, 3)) +
    gf_avg_wilson_line(gf, make_wilson_line_path(l, t, 2, 3));
  return m / 3.0;
}

void gf_show_info(const GaugeField& gf, const int level = 0)
{
  TIMER_VERBOSE("gf_show_info");
  displayln_info(shows("plaq : ") + ::show(gf_avg_plaq(gf)));
  displayln_info(shows("trace: ") + ::show(gf_avg_link_trace(gf)));
  if (0 < level) {
    displayln_info(shows("plaq 1x1 : ") + ::show(matrix_trace(gf_avg_wilson_loop(gf, 1, 1)).real() / 3.0));
    displayln_info(shows("plaq 1x2 : ") + ::show(matrix_trace(gf_avg_wilson_loop(gf, 1, 2)).real() / 3.0));
    displayln_info(shows("plaq 2x1 : ") + ::show(matrix_trace(gf_avg_wilson_loop(gf, 2, 1)).real() / 3.0));
    displayln_info(shows("plaq 2x2 : ") + ::show(matrix_trace(gf_avg_wilson_loop(gf, 2, 2)).real() / 3.0));
    displayln_info(shows("plaq 1x3 : ") + ::show(matrix_trace(gf_avg_wilson_loop(gf, 1, 3)).real() / 3.0));
    displayln_info(shows("plaq 3x1 : ") + ::show(matrix_trace(gf_avg_wilson_loop(gf, 3, 1)).real() / 3.0));
    displayln_info(shows("plaq 2x3 : ") + ::show(matrix_trace(gf_avg_wilson_loop(gf, 2, 3)).real() / 3.0));
    displayln_info(shows("plaq 3x2 : ") + ::show(matrix_trace(gf_avg_wilson_loop(gf, 3, 2)).real() / 3.0));
    displayln_info(shows("plaq 3x3 : ") + ::show(matrix_trace(gf_avg_wilson_loop(gf, 3, 3)).real() / 3.0));
  }
}

inline WilsonLinePath add_to_path(const WilsonLinePath& path, const int dir)
{
  if (WilsonLinePathType_simp == path.type) {
    WilsonLinePath p = path;
    p.simp.push_back(dir);
    p.target = coordinate_shifts(p.target, dir);
    return p;
  } else if (WilsonLinePathType_comp_par == path.type) {
    WilsonLinePath ps;
    ps.simp.push_back(dir);
    ps.target = coordinate_shifts(ps.target, dir);
    WilsonLinePath p;
    p.type = WilsonLinePathType_comp_seq;
    p.comp.push_back(path);
    p.comp.push_back(ps);
    p.target = coordinate_shifts(path.target, dir);
    return p;
  } else if (WilsonLinePathType_comp_par_multi_target == path.type) {
    WilsonLinePath p = path;
    qassert(2 <= p.comp.size());
    for (int i = 0; i < p.comp.size(); ++i) {
      p.comp[i] = add_to_path(p.comp[i], dir);
    }
    bool same = true;
    for (int i = 1; i < p.comp.size(); ++i) {
      same = same && p.comp[i-1].target == p.comp[i].target;
    }
    if (same) {
      p.target = p.comp[0].target;
      p.type = WilsonLinePathType_comp_par;
    }
    return p;
  } else if (WilsonLinePathType_comp_seq == path.type) {
    WilsonLinePath p = path;
    qassert(0 != p.comp.size());
    if (p.comp.back().type == WilsonLinePathType_simp) {
      p.comp.back() = add_to_path(p.comp.back(), dir);
      p.target = coordinate_shifts(p.target, dir);
      return p;
    } else if (p.comp.back().type == WilsonLinePathType_comp_par) {
      WilsonLinePath ps;
      ps.simp.push_back(dir);
      ps.target = coordinate_shifts(ps.target, dir);
      p.comp.push_back(ps);
      p.target = coordinate_shifts(p.target, dir);
      return p;
    } else if (p.comp.back().type == WilsonLinePathType_comp_par_multi_target) {
      p.comp.back() = add_to_path(p.comp.back(), dir);
      if (WilsonLinePathType_comp_par == p.comp.back().type) {
        p.target = p.target + p.comp.back().target;
      }
      return p;
    } else {
      // comp_seq in comp_seq
      qassert(false);
    }
  }
}

inline WilsonLinePath make_init_wilson_line_path(const int dir)
{
  WilsonLinePath path;
  path.simp.push_back(dir);
  path.target = coordinate_shifts(Coordinate(), dir);
  return path;
}

inline WilsonLinePath make_init_wilson_line_path(const std::vector<int> dirs)
{
  qassert(0 < dirs.size());
  if (dirs.size() == 1) {
    return make_init_wilson_line_path(dirs[0]);
  } else {
    WilsonLinePath path;
    path.type = WilsonLinePathType_comp_par_multi_target;
    for (int i = 0; i < dirs.size(); ++i) {
      path.comp.push_back(make_init_wilson_line_path(dirs[i]));
    }
    return path;
  }
}

void extend_wilson_line_path_from_start_to_target(WilsonLinePath& path, const Coordinate& start_loc, const Coordinate& target_wilson_line)
{
  if (path.type == WilsonLinePathType_simp) {
    qassert(false);
  } else if (path.type == WilsonLinePathType_comp_seq) {
    if (path.comp.size() == 0) {
    }
  } else if (path.type == WilsonLinePathType_comp_par) {
  } else if (path.type == WilsonLinePathType_comp_par_multi_target) {
  }
}

inline WilsonLinePath make_wilson_line_path_todo(const Coordinate& target_wilson_line)
{
  int steps = 0;
  for (int i = 0; i < DIM; ++i) {
    const int x = target_wilson_line[i];
    steps += std::abs(x);
  }
  const double eps_dis = 1.0e-8;
  WilsonLinePath path;
  for (int i = 0; i < steps; ++i) {
    WilsonLinePath* ppath = &path;
    bool is_multi_target = false;
    while (ppath != NULL) {
      if (WilsonLinePathType_comp_par_multi_target == ppath->type) {
        is_multi_target = true;
        break;
      } else if (WilsonLinePathType_comp_par == ppath->type) {
        is_multi_target = false;
        break;
      } else if (WilsonLinePathType_comp_seq == ppath->type) {
        ppath = &path.comp.back();
      } else if (WilsonLinePathType_simp == ppath->type) {
        is_multi_target = false;
        break;
      }
    }
    if (!is_multi_target) { // FIXME
      const Coordinate& loc = path.target;
      std::vector<int> next_dirs = find_next_dirs(path.target, target_wilson_line);
      // update path;
      qassert(1 <= next_dirs.size());
      if (1 == next_dirs.size()) {
        path = add_to_path(path, next_dirs[0]);
      } else {
        WilsonLinePath pp;
        pp.type = WilsonLinePathType_comp_par_multi_target;
        pp.comp.resize(next_dirs.size());
        for (int i = 0; i < pp.comp.size(); ++i) {
          pp.comp[i] = add_to_path(WilsonLinePath(), next_dirs[i]);
        }
        if (WilsonLinePathType_comp_seq == path.type) {
          path.comp.push_back(pp);
        } else if (WilsonLinePathType_simp == path.type || WilsonLinePathType_comp_par == path.type) {
          WilsonLinePath old_path = path;
          path.type = WilsonLinePathType_comp_seq;
          path.simp.clear();
          path.comp.push_back(old_path);
          path.comp.push_back(pp);
        } else {
          //
        }
      }
    } else {
      //
    }
  }
  return path;
}

QLAT_END_NAMESPACE
