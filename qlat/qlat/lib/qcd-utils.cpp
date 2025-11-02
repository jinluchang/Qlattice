#include <qlat/qcd-acc.h>
#include <qlat/qcd-utils.h>

namespace qlat
{  //

// BUG ! path with std::vector using in qacc
void gf_wilson_line_no_comm(Field<ColorMatrix>& wilson_line_field,
                            const Int wilson_line_field_m,
                            const GaugeField& gf_ext,
                            const std::vector<Int>& path)
// wilson_line_field needs to be initialized before hand
// 0 <= wilson_line_field_m < wilson_line_field.multiplicity
{
  TIMER("gf_wilson_line_no_comm")
  const Geometry& geo = wilson_line_field.geo();
  const Int multiplicity = wilson_line_field.multiplicity;
  Qassert(check_matching_geo(geo, gf_ext.geo()));
  Qassert(0 <= wilson_line_field_m and wilson_line_field_m < multiplicity);
  qthread_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = wilson_line_field.get_elems(xl);
    v[wilson_line_field_m] = gf_wilson_line_no_comm(gf_ext, xl, path);
  });
}

void gf_wilson_line_no_comm(Field<ColorMatrix>& wilson_line_field,
                            const Int wilson_line_field_m,
                            const GaugeField& gf_ext,
                            const std::vector<Int>& path,
                            const std::vector<Int>& path_n)
// wilson_line_field needs to be initialized before hand
// 0 <= wilson_line_field_m < wilson_line_field.multiplicity
{
  TIMER("gf_wilson_line_no_comm")
  const Geometry& geo = wilson_line_field.geo();
  const Int multiplicity = wilson_line_field.multiplicity;
  Qassert(check_matching_geo(geo, gf_ext.geo()));
  Qassert(0 <= wilson_line_field_m and wilson_line_field_m < multiplicity);
  qthread_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = wilson_line_field.get_elems(xl);
    v[wilson_line_field_m] = gf_wilson_line_no_comm(gf_ext, xl, path, path_n);
  });
}

void set_g_rand_anti_hermitian_matrix_field(Field<ColorMatrix>& fc,
                                            const RngState& rs,
                                            const double sigma)
//  Creates a field of anti-hermitian 3x3 complex matrices with each complex
//  element drawn at random from a Gaussian distribution with zero mean.
//  Hence the matrices are distributed according to
//
//  exp[- Tr(mat^2)/(2 sigma**2)]
{
  TIMER("set_g_rand_anti_hermitian_matrix_field");
  const Geometry& geo = fc.geo();
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Long gindex = geo.g_index_from_g_coordinate(xg);
    RngState rsi(rs, gindex);
    Vector<ColorMatrix> v = fc.get_elems(xl);
    for (Int m = 0; m < (int)v.size(); ++m) {
      v[m] = make_g_rand_anti_hermitian_matrix(rsi, sigma);
    }
  }
}

void set_g_rand_color_matrix_field(Field<ColorMatrix>& fc, const RngState& rs,
                                   const double sigma, const Int n_step)
{
  TIMER("set_g_rand_color_matrix_field");
  const Geometry& geo = fc.geo();
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Long gindex = geo.g_index_from_g_coordinate(xg);
    RngState rsi(rs, gindex);
    Vector<ColorMatrix> v = fc.get_elems(xl);
    for (Int m = 0; m < (int)v.size(); ++m) {
      v[m] =
          make_color_matrix_exp(make_g_rand_anti_hermitian_matrix(rsi, sigma));
      for (Int k = 1; k < n_step; ++k) {
        v[m] *= make_color_matrix_exp(
            make_g_rand_anti_hermitian_matrix(rsi, sigma));
      }
    }
  }
}

void set_local_current_from_props(FieldM<WilsonMatrix, 4>& cf,
                                  const Propagator4d& prop1,
                                  const Propagator4d& prop2)
// ->- prop1 ->- gamma_mu ->- gamma5 prop2^+ gamma5 ->-
{
  TIMER_VERBOSE("set_local_current_from_props");
  const Geometry geo = geo_resize(prop1.geo());
  Qassert(geo == geo_resize(prop2.geo()));
  const array<SpinMatrix, 4>& gammas =
      SpinMatrixConstants::get_cps_gammas();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  cf.init(geo);
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const WilsonMatrix& m1 = prop1.get_elem(xl);
    const WilsonMatrix& m2 = prop2.get_elem(xl);
    Vector<WilsonMatrix> v = cf.get_elems(xl);
    const WilsonMatrix m2rev =
        gamma5 * (WilsonMatrix)matrix_adjoint(m2) * gamma5;
    for (Int m = 0; m < 4; ++m) {
      v[m] = m2rev * gammas[m] * m1;
    }
  }
}

double coordinate_distance_from_wilson_line(
    const Coordinate& x, const Coordinate& target_wilson_line)
{
  const CoordinateD target(target_wilson_line);
  const double target_len = coordinate_len(target);
  const CoordinateD vdir = target / target_len;
  const CoordinateD xd(x);
  const CoordinateD dis = xd - dot_product(xd, vdir) * vdir;
  return coordinate_len(dis);
}

std::vector<int> find_next_dirs(const Coordinate& loc,
                                const Coordinate& target_wilson_line)
{
  TIMER("find_next_dirs");
  if (target_wilson_line == loc) {
    return std::vector<int>();
  }
  std::vector<int> dirs;
  for (Int i = 0; i < DIMN; ++i) {
    const Int x = target_wilson_line[i];
    if (0 < x) {
      dirs.push_back(i);
    } else if (0 > x) {
      dirs.push_back(-i - 1);
    }
  }
  Qassert(0 < dirs.size());
  if (1 == dirs.size()) {
    return dirs;
  } else {
    const double eps_dis = 1.0e-8;
    std::vector<int> next_dirs;
    const double dis =
        coordinate_distance_from_wilson_line(loc, target_wilson_line);
    // ADJUST ME ; whether to choose more route to average
    if (false && dis < eps_dis) {
      next_dirs = dirs;
    } else {
      double min_dis = 1e6;  // a very large number
      for (Int i = 0; i < (int)dirs.size(); ++i) {
        const double next_dis = coordinate_distance_from_wilson_line(
            coordinate_shifts(loc, dirs[i]), target_wilson_line);
        if (next_dis < min_dis) {
          min_dis = next_dis;
        }
      }
      for (Int i = 0; i < (int)dirs.size(); ++i) {
        const double next_dis = coordinate_distance_from_wilson_line(
            coordinate_shifts(loc, dirs[i]), target_wilson_line);
        if (next_dis <= min_dis + eps_dis) {
          next_dirs.push_back(dirs[i]);
        }
      }
    }
    return next_dirs;
  }
}

void acc_wilson_line_path_segment(WilsonLinePathSegment& path)
{
  TIMER("acc_wilson_line_path_segment");
  std::map<Coordinate, WilsonLinePathStop>& stops = path.stops;
  std::set<Coordinate> cset;
  std::vector<Coordinate> cs;
  cs.push_back(Coordinate());
  for (Int i = 0; i < (int)cs.size(); ++i) {
    const Coordinate c = cs[i];
    cset.insert(c);
    const WilsonLinePathStop& ps = stops[c];
    Qassert(c == ps.x);
    for (Int k = 0; k < (int)ps.paths.size(); ++k) {
      const Coordinate nc = coordinate_shifts(c, ps.paths[k]);
      cs.push_back(nc);
      stops[nc].num_origins = 0;
    }
  }
  clear(cs);
  for (std::set<Coordinate>::iterator it = cset.begin(); it != cset.end();
       ++it) {
    cs.push_back(*it);
  }
  for (Int i = 0; i < (int)cs.size(); ++i) {
    const Coordinate& c = cs[i];
    const WilsonLinePathStop& ps = stops[c];
    Qassert(c == ps.x);
    for (Int k = 0; k < (int)ps.paths.size(); ++k) {
      const Coordinate nc = coordinate_shifts(c, ps.paths[k]);
      stops[nc].num_origins += 1;
    }
  }
}

WilsonLinePathSegment make_wilson_line_path_segment(const Coordinate& target)
{
  TIMER("make_wilson_line_path_segment");
  WilsonLinePathSegment ret;
  ret.target = target;
  std::vector<Coordinate> cs;
  cs.push_back(Coordinate());
  for (Int i = 0; i < (int)cs.size(); ++i) {
    const Coordinate c = cs[i];
    std::vector<int> dirs = find_next_dirs(c, target);
    WilsonLinePathStop& ps = ret.stops[c];
    ps.x = c;
    ps.paths.resize(dirs.size());
    for (Int k = 0; k < (int)dirs.size(); ++k) {
      Qassert(ps.paths[k].size() == 0);
      Int dir = dirs[k];
      ps.paths[k].push_back(dir);
      Coordinate nc = coordinate_shifts(c, dir);
      Qassert(nc == coordinate_shifts(c, ps.paths[k]));
      std::vector<int> ndirs = find_next_dirs(nc, target);
      while (ndirs.size() == 1 && ret.stops.find(nc) == ret.stops.end()) {
        dir = ndirs[0];
        ps.paths[k].push_back(dir);
        nc = coordinate_shifts(nc, dir);
        Qassert(nc == coordinate_shifts(c, ps.paths[k]));
        ndirs = find_next_dirs(nc, target);
      }
      Qassert(nc == coordinate_shifts(c, ps.paths[k]));
      if (ret.stops.find(nc) == ret.stops.end()) {
        ret.stops[nc].x = nc;
        cs.push_back(nc);
      }
    }
  }
  acc_wilson_line_path_segment(ret);
  return ret;
}

void set_multiply_simple_wilson_line_field_partial_comm(
    FieldM<ColorMatrix, 1>& wlf, FieldM<ColorMatrix, 1>& wlf1,
    const GaugeField& gf1, const std::vector<int>& path)
// gf1 need to be refresh_expanded_1.
// wlf1 need to have correct size: wlf1.init(geo_resize(geo, 1));
// wlf1 will be modified
// wlf will be initialized
{
  TIMER("set_multiply_simple_wilson_line_field_partial_comm");
  const Geometry geo = geo_resize(gf1.geo());
  Qassert(&wlf != &wlf1);
  wlf.init(geo);
  for (size_t i = 0; i < path.size(); ++i) {
    const Int dir = path[i];
    Qassert(-DIMN <= dir && dir < DIMN);
    refresh_expanded_1(wlf1);
#pragma omp parallel for
    for (Long index = 0; index < geo.local_volume(); ++index) {
      Coordinate xl = geo.coordinate_from_index(index);
      ColorMatrix& l1 = wlf.get_elem(xl);
      if (0 <= dir) {
        xl[dir] -= 1;
        const ColorMatrix& link = gf1.get_elem(xl, dir);
        const ColorMatrix& l0 = wlf1.get_elem(xl);
        l1 = l0 * link;
      } else {
        const ColorMatrix& link = gf1.get_elem(xl, -dir - 1);
        xl[-dir - 1] += 1;
        const ColorMatrix& l0 = wlf1.get_elem(xl);
        l1 = l0 * matrix_adjoint(link);
      }
    }
    if (i != path.size() - 1) {
      wlf1 = wlf;
    }
  }
}

void set_multiply_wilson_line_field_partial_comm(
    FieldM<ColorMatrix, 1>& wlf, FieldM<ColorMatrix, 1>& wlf1,
    const GaugeField& gf1, const WilsonLinePathSegment& path)
// gf1 need to be refresh_expanded_1.
// wlf1 need to have correct size: wlf1.init(geo_resize(geo, 1));
// wlf1 will be modified
// wlf will be initialized
{
  TIMER("set_multiply_wilson_line_field_partial_comm");
  WilsonLinePathSegment pacc = path;
  if (pacc.stops[Coordinate()].paths.size() == 1 && pacc.stops.size() == 2) {
    set_multiply_simple_wilson_line_field_partial_comm(
        wlf, wlf1, gf1, pacc.stops[Coordinate()].paths[0]);
    return;
  }
  const Geometry geo = geo_resize(gf1.geo());
  std::vector<Coordinate> cs;
  std::vector<FieldM<ColorMatrix, 1> > fs(pacc.stops.size());
  std::map<Coordinate, int> dict;
  cs.push_back(Coordinate());
  fs[0].init(geo);
  fs[0] = wlf1;
  dict[cs[0]] = 0;
  wlf.init(geo);
  while (true) {
    for (Int i = 0; i < (int)cs.size(); ++i) {
      const Coordinate c = cs[i];
      const WilsonLinePathStop& ps = pacc.stops[c];
      Qassert(c == ps.x);
      if (is_initialized(fs[i]) && ps.num_origins == 0) {
        if (c == pacc.target) {
          Qassert(ps.paths.size() == 0);
          wlf = fs[i];
          return;
        }
        for (Int k = 0; k < (int)ps.paths.size(); ++k) {
          wlf1 = fs[i];
          const Coordinate nc = coordinate_shifts(c, ps.paths[k]);
          pacc.stops[nc].num_origins -= 1;
          Handle<FieldM<ColorMatrix, 1> > hf;
          if (dict.find(nc) == dict.end()) {
            cs.push_back(nc);
            dict[nc] = cs.size() - 1;
            hf.init(fs[dict[nc]]);
            set_multiply_simple_wilson_line_field_partial_comm(hf(), wlf1, gf1,
                                                               ps.paths[k]);
          } else {
            set_multiply_simple_wilson_line_field_partial_comm(wlf, wlf1, gf1,
                                                               ps.paths[k]);
            hf.init(fs[dict[nc]]);
            hf() += wlf;
          }
          Qassert(cs[dict[nc]] == nc);
        }
        fs[i].init();
      }
    }
  }
}

void set_left_expanded_gauge_field(GaugeField& gf1, const GaugeField& gf)
{
  TIMER_VERBOSE("set_left_expanded_gauge_field");
  const Coordinate expansion_left(1, 1, 1, 1);
  const Coordinate expansion_right(0, 0, 0, 0);
  const Geometry geo1 = geo_resize(gf.geo(), expansion_left, expansion_right);
  gf1.init(geo1);
  Qassert(gf1.geo() == geo1);
  gf1 = gf;
  refresh_expanded_1(gf1);
}

ColorMatrix gf_avg_wilson_line(const GaugeField& gf, const WilsonLinePath& path)
{
  TIMER("gf_avg_wilson_line");
  const Geometry geo = geo_resize(gf.geo());
  const Coordinate expansion_left(1, 1, 1, 1);
  const Coordinate expansion_right(0, 0, 0, 0);
  GaugeField gf1;
  set_left_expanded_gauge_field(gf1, gf);
  FieldM<ColorMatrix, 1> wlf, wlf1;
  wlf1.init(geo_resize(geo, 1));
  wlf.init(geo);
  set_unit(wlf1);
  for (Int i = 0; i < (int)path.ps.size(); ++i) {
    set_multiply_wilson_line_field_partial_comm(wlf, wlf1, gf1, path.ps[i]);
    if (i != (int)path.ps.size() - 1) {
      wlf1 = wlf;
    }
  }
  return field_glb_sum(wlf)[0] / (ComplexD)geo.total_volume();
}

WilsonLinePath make_wilson_loop_path(const Coordinate& target_l, const Int t)
{
  const Coordinate target_t(0, 0, 0, t);
  WilsonLinePath path;
  path.ps.push_back(make_wilson_line_path_segment(target_l));
  path.ps.push_back(make_wilson_line_path_segment(target_t));
  path.ps.push_back(make_wilson_line_path_segment(-target_l));
  path.ps.push_back(make_wilson_line_path_segment(-target_t));
  return path;
}

ColorMatrix gf_avg_wilson_loop(const GaugeField& gf, const Int l, const Int t)
{
  TIMER_VERBOSE("gf_avg_wilson_loop");
  ColorMatrix m =
      gf_avg_wilson_line(gf, make_wilson_loop_path(Coordinate(l, 0, 0, 0), t)) +
      gf_avg_wilson_line(gf, make_wilson_loop_path(Coordinate(0, l, 0, 0), t)) +
      gf_avg_wilson_line(gf, make_wilson_loop_path(Coordinate(0, 0, l, 0), t));
  return m / ComplexD(3.0);
}

std::vector<Coordinate> spatial_permute_direction(const Coordinate& l)
{
  const Int x = l[0];
  const Int y = l[1];
  const Int z = l[2];
  std::set<Coordinate> cset;
  cset.insert(Coordinate(x, y, z, 0));
  cset.insert(Coordinate(y, z, x, 0));
  cset.insert(Coordinate(z, x, y, 0));
  cset.insert(Coordinate(z, y, x, 0));
  cset.insert(Coordinate(y, x, z, 0));
  cset.insert(Coordinate(x, z, y, 0));
  cset.insert(Coordinate(-x, y, z, 0));
  cset.insert(Coordinate(-y, z, x, 0));
  cset.insert(Coordinate(-z, x, y, 0));
  cset.insert(Coordinate(-z, y, x, 0));
  cset.insert(Coordinate(-y, x, z, 0));
  cset.insert(Coordinate(-x, z, y, 0));
  cset.insert(Coordinate(x, -y, z, 0));
  cset.insert(Coordinate(y, -z, x, 0));
  cset.insert(Coordinate(z, -x, y, 0));
  cset.insert(Coordinate(z, -y, x, 0));
  cset.insert(Coordinate(y, -x, z, 0));
  cset.insert(Coordinate(x, -z, y, 0));
  cset.insert(Coordinate(x, y, -z, 0));
  cset.insert(Coordinate(y, z, -x, 0));
  cset.insert(Coordinate(z, x, -y, 0));
  cset.insert(Coordinate(z, y, -x, 0));
  cset.insert(Coordinate(y, x, -z, 0));
  cset.insert(Coordinate(x, z, -y, 0));
  std::vector<Coordinate> cs;
  for (std::set<Coordinate>::iterator it = cset.begin(); it != cset.end();
       ++it) {
    cs.push_back(-*it);
  }
  for (Int i = 0; i < (int)cs.size(); ++i) {
    cset.insert(cs[i]);
  }
  clear(cs);
  for (std::set<Coordinate>::iterator it = cset.begin(); it != cset.end();
       ++it) {
    cs.push_back(*it);
  }
  return cs;
}

ColorMatrix gf_avg_wilson_loop(const GaugeField& gf, const Coordinate& l,
                               const Int t)
{
  TIMER_VERBOSE("gf_avg_wilson_loop(Coordinate&)");
  ColorMatrix m;
  set_zero(m);
  std::vector<Coordinate> cs = spatial_permute_direction(l);
  for (Int i = 0; i < (int)cs.size(); ++i) {
    m += gf_avg_wilson_line(gf, make_wilson_loop_path(cs[i], t));
  }
  return m;
}

void gf_show_info(const GaugeField& gf, const Int level)
{
  TIMER_VERBOSE("gf_show_info");
  displayln_info(shows("plaq : ") + show(gf_avg_plaq(gf)));
  displayln_info(shows("spatial plaq : ") + show(gf_avg_spatial_plaq(gf)));
  displayln_info(shows("trace: ") + show(gf_avg_link_trace(gf)));
  if (0 < level) {
    displayln_info(
        shows("plaq 1x1 : ") +
        show(matrix_trace(gf_avg_wilson_loop(gf, 1, 1)).real() / 3.0));
    displayln_info(
        shows("plaq 1x2 : ") +
        show(matrix_trace(gf_avg_wilson_loop(gf, 1, 2)).real() / 3.0));
    displayln_info(
        shows("plaq 2x1 : ") +
        show(matrix_trace(gf_avg_wilson_loop(gf, 2, 1)).real() / 3.0));
    displayln_info(
        shows("plaq 2x2 : ") +
        show(matrix_trace(gf_avg_wilson_loop(gf, 2, 2)).real() / 3.0));
    displayln_info(
        shows("plaq 1x3 : ") +
        show(matrix_trace(gf_avg_wilson_loop(gf, 1, 3)).real() / 3.0));
    displayln_info(
        shows("plaq 3x1 : ") +
        show(matrix_trace(gf_avg_wilson_loop(gf, 3, 1)).real() / 3.0));
    displayln_info(
        shows("plaq 2x3 : ") +
        show(matrix_trace(gf_avg_wilson_loop(gf, 2, 3)).real() / 3.0));
    displayln_info(
        shows("plaq 3x2 : ") +
        show(matrix_trace(gf_avg_wilson_loop(gf, 3, 2)).real() / 3.0));
    displayln_info(
        shows("plaq 3x3 : ") +
        show(matrix_trace(gf_avg_wilson_loop(gf, 3, 3)).real() / 3.0));
    displayln_info(
        shows("plaq (1,1,0,0)x2 : ") +
        show(matrix_trace(gf_avg_wilson_loop(gf, Coordinate(1, 1, 0, 0), 2))
                 .real() /
             3.0));
    displayln_info(
        shows("plaq (2,1,0,0)x1 : ") +
        show(matrix_trace(gf_avg_wilson_loop(gf, Coordinate(1, 1, 0, 0), 1))
                 .real() /
             3.0));
  }
}

}  // namespace qlat
