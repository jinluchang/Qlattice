#pragma once

#include <qlat-utils/matrix.h>
#include <qlat/qcd-utils.h>
#include <qlat/qcd.h>

namespace qlat
{  //

qacc ColorMatrix gf_clover_leaf_no_comm(const GaugeField& gf1,
                                        const Coordinate& xl, const int mu,
                                        const int nu)
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
                                        const Coordinate& xl, const int mu)
// transpose the same way as gf.get_elem(xl, mu)
{
  ColorMatrix acc;
  set_zero(acc);
  for (int nu = -4; nu < 4; ++nu) {
    if (nu == mu or -nu - 1 == mu) {
      continue;
    }
    acc += gf_wilson_line_no_comm(gf, xl, make_array<int>(nu, mu, -nu - 1));
  }
  return acc;
}

}  // namespace qlat
