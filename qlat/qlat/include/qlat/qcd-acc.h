#pragma once

#include <qlat-utils/matrix.h>
#include <qlat/qcd-utils.h>
#include <qlat/qcd.h>

namespace qlat
{  //

qacc ColorMatrix gf_get_link(const GaugeField& gf, const Coordinate& xl,
                             const Int mu)
// mu can be negative
{
  if (0 <= mu) {
    return gf.get_elem(xl, mu);
  } else {
    const Coordinate xl1 = coordinate_shifts(xl, mu);
    return matrix_adjoint(gf.get_elem(xl1, -mu - 1));
  }
}

template <class Vec>
qacc ColorMatrix gf_wilson_line_no_comm(const GaugeField& gf,
                                        const Coordinate& xl, const Vec& path)
{
  ColorMatrix ret;
  set_unit(ret);
  Coordinate xl1 = xl;
  for (Int i = 0; i < (int)path.size(); ++i) {
    const Int dir = path[i];
    qassert(-DIMN <= dir && dir < DIMN);
    if (0 <= dir) {
      ret *= gf.get_elem(xl1, dir);
      xl1[dir] += 1;
    } else {
      xl1[-dir - 1] -= 1;
      ret *= matrix_adjoint(gf.get_elem(xl1, -dir - 1));
    }
  }
  return ret;
}

template <class Vec>
qacc ColorMatrix gf_wilson_line_no_comm(const GaugeField& gf,
                                        const Coordinate& xl, const Vec& path,
                                        const Vec& path_n)
{
  qassert((Long)path.size() == (Long)path_n.size());
  ColorMatrix ret;
  set_unit(ret);
  Coordinate xl1 = xl;
  for (Int i = 0; i < (int)path.size(); ++i) {
    const Int dir = path[i];
    qassert(-DIMN <= dir && dir < DIMN);
    for (Int j = 0; j < (int)path_n[i]; ++j) {
      if (0 <= dir) {
        ret *= gf.get_elem(xl1, dir);
        xl1[dir] += 1;
      } else {
        xl1[-dir - 1] -= 1;
        ret *= matrix_adjoint(gf.get_elem(xl1, -dir - 1));
      }
    }
  }
  return ret;
}

qacc ColorMatrix gf_staple_no_comm_v1(const GaugeField& gf,
                                      const Coordinate& xl, const Int mu)
{
  ColorMatrix ret;
  set_zero(ret);
  const Coordinate xl_mu = coordinate_shifts(xl, mu);
  for (Int m = 0; m < DIMN; ++m) {
    if (mu != m) {
      ret += gf.get_elem(xl, m) * gf.get_elem(coordinate_shifts(xl, m), mu) *
             matrix_adjoint(gf.get_elem(xl_mu, m));
      ret += matrix_adjoint(gf.get_elem(coordinate_shifts(xl, -m - 1), m)) *
             gf.get_elem(coordinate_shifts(xl, -m - 1), mu) *
             gf.get_elem(coordinate_shifts(xl_mu, -m - 1), m);
    }
  }
  return ret;
}

qacc ColorMatrix gf_staple_no_comm_v2(const GaugeField& gf,
                                      const Coordinate& xl, const Int mu)
{
  ColorMatrix ret;
  set_zero(ret);
  array<int, 3> path;
  path[1] = mu;
  for (Int m = 0; m < DIMN; ++m) {
    if (mu != m) {
      path[0] = m;
      path[2] = -m - 1;
      ret += gf_wilson_line_no_comm(gf, xl, path);
    }
  }
  for (Int m = 0; m < DIMN; ++m) {
    if (mu != m) {
      path[0] = -m - 1;
      path[2] = m;
      ret += gf_wilson_line_no_comm(gf, xl, path);
    }
  }
  return ret;
}

qacc ColorMatrix gf_staple_no_comm(const GaugeField& gf, const Coordinate& xl,
                                   const Int mu)
{
  return gf_staple_no_comm_v1(gf, xl, mu);
  // return gf_staple_no_comm_v2(gf, xl, mu);
}

qacc ColorMatrix gf_spatial_staple_no_comm(const GaugeField& gf,
                                           const Coordinate& xl, const Int mu)
{
  ColorMatrix ret;
  set_zero(ret);
  const Coordinate xl_mu = coordinate_shifts(xl, mu);
  for (Int m = 0; m < 3; ++m) {
    if (mu != m) {
      ret += gf.get_elem(xl, m) * gf.get_elem(coordinate_shifts(xl, m), mu) *
             matrix_adjoint(gf.get_elem(xl_mu, m));
      ret += matrix_adjoint(gf.get_elem(coordinate_shifts(xl, -m - 1), m)) *
             gf.get_elem(coordinate_shifts(xl, -m - 1), mu) *
             gf.get_elem(coordinate_shifts(xl_mu, -m - 1), m);
    }
  }
  return ret;
}

qacc ColorMatrix gf_clover_leaf_no_comm(const GaugeField& gf1,
                                        const Coordinate& xl, const Int mu,
                                        const Int nu)
{
  ColorMatrix m;
  set_zero(m);
  m += gf_wilson_line_no_comm(gf1, xl,
                              make_array<int>(mu, nu, -mu - 1, -nu - 1));
  m += gf_wilson_line_no_comm(gf1, xl,
                              make_array<int>(-mu - 1, -nu - 1, mu, nu));
  m += gf_wilson_line_no_comm(gf1, xl,
                              make_array<int>(nu, -mu - 1, -nu - 1, mu));
  m += gf_wilson_line_no_comm(gf1, xl,
                              make_array<int>(-nu - 1, mu, nu, -mu - 1));
  return (ComplexD)0.25 * m;
}

qacc ColorMatrix gf_plaq_staple_no_comm(const GaugeField& gf,
                                        const Coordinate& xl, const Int mu)
// transpose the same way as gf.get_elem(xl, mu)
{
  ColorMatrix acc;
  set_zero(acc);
  for (Int nu = -4; nu < 4; ++nu) {
    if (nu == mu or -nu - 1 == mu) {
      continue;
    }
    acc += gf_wilson_line_no_comm(gf, xl, make_array<int>(nu, mu, -nu - 1));
  }
  return acc;
}

}  // namespace qlat
