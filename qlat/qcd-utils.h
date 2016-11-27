#pragma once

#include <qlat/config.h>
#include <qlat/field.h>

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

inline ColorMatrix gf_avg_wilson_line(const GaugeField& gf, const std::vector<int>& path)
{
  TIMER("gf_avg_wilson_line");
  const Geometry geo = geo_reform(gf.geo);
  const Geometry geo1 = geo_reform(gf.geo, 1, 1);
  const Coordinate expansion_left(1, 1, 1, 1);
  const Coordinate expansion_right(0, 0, 0, 0);
  GaugeField gf1;
  gf1.init(geo_resize(geo, expansion_left, expansion_right));
  gf1 = gf;
  refresh_expanded(gf1);
  FieldM<ColorMatrix,1> wlf0, wlf1;
  wlf0.init(geo1);
  wlf1.init(geo);
  set_unit(wlf0);
  for (size_t i = 0; i < path.size(); ++i) {
    const int dir = path[i];
    qassert(-DIM <= dir && dir < DIM);
    refresh_expanded(wlf0);
#pragma omp parallel for
    for (long index = 0; index < geo.local_volume(); ++index) {
      Coordinate xl = geo.coordinate_from_index(index);
      ColorMatrix& l1 = wlf1.get_elem(xl);
      if (0 <= dir) {
        xl[dir] -= 1;
        const ColorMatrix& link = gf1.get_elem(xl, dir);
        const ColorMatrix& l0 = wlf0.get_elem(xl);
        l1 = l0 * link;
      } else {
        const ColorMatrix& link = gf1.get_elem(xl, -dir-1);
        xl[-dir-1] += 1;
        const ColorMatrix& l0 = wlf0.get_elem(xl);
        l1 = l0 * matrix_adjoint(link);
      }
    }
    wlf0 = wlf1;
  }
  return field_glb_sum_double(wlf0)[0] / geo.total_volume();
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
  displayln_info(ssprintf("%s: l=%3d; t=%3d.", fname, l, t));
  ColorMatrix m =
    gf_avg_wilson_line(gf, make_wilson_line_path(l, t, 0, 3)) +
    gf_avg_wilson_line(gf, make_wilson_line_path(l, t, 1, 3)) +
    gf_avg_wilson_line(gf, make_wilson_line_path(l, t, 2, 3));
  return m / 3.0;
}

QLAT_END_NAMESPACE
