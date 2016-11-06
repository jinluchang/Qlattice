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

QLAT_END_NAMESPACE
