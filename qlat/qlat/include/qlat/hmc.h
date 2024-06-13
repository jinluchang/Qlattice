#pragma once

#include <qlat/gauge-action.h>
#include <qlat/scalar-action.h>
#include <qlat/qcd-smear.h>
#include <qlat/qcd-topology.h>
#include <qlat/qcd-utils.h>
#include <qlat/qcd.h>

#include <cmath>
#include <sstream>
#include <string>

namespace qlat
{  //

struct API GaugeMomentum : FieldM<ColorMatrix, 4> {
};

bool metropolis_accept(double& accept_prob, const double delta_h,
                       const int traj, const RngState& rs_);

void set_rand_gauge_momentum(GaugeMomentum& gm, const double sigma,
                             const RngState& rs);

double gm_hamilton_node(const GaugeMomentum& gm);

double gf_sum_re_tr_plaq_node_no_comm(const GaugeField& gf);

double gf_sum_re_tr_rect_node_no_comm(const GaugeField& gf);

double gf_hamilton_node_no_comm(const GaugeField& gf, const GaugeAction& ga);

double gf_hamilton_node(const GaugeField& gf, const GaugeAction& ga);

double gf_hamilton(const GaugeField& gf, const GaugeAction& ga);

void gf_evolve(GaugeField& gf, const GaugeMomentum& gm, const double step_size);

void gf_evolve_dual(GaugeField& gf, const GaugeMomentum& gm_dual,
                    const double step_size);

void set_gm_force_no_comm(GaugeMomentum& gm_force, const GaugeField& gf,
                          const GaugeAction& ga);

void set_gm_force(GaugeMomentum& gm_force, const GaugeField& gf,
                  const GaugeAction& ga);

void set_gm_force_dual(GaugeMomentum& gm_force_dual, const GaugeField& gf,
                       const GaugeMomentum& gm_force);

// -------------------

qacc double gf_re_tr_plaq_no_comm(const GaugeField& gf, const Coordinate& xl,
                                  const int mu, const int nu)
{
  const ColorMatrix m =
      gf_wilson_line_no_comm(gf, xl, make_array<int>(mu, nu, -mu - 1, -nu - 1));
  return matrix_trace(m).real();
}

qacc double gf_re_tr_rect_no_comm(const GaugeField& gf, const Coordinate& xl,
                                  const int mu, const int nu)
{
  const ColorMatrix m = gf_wilson_line_no_comm(
      gf, xl, make_array<int>(mu, mu, nu, -mu - 1, -mu - 1, -nu - 1));
  return matrix_trace(m).real();
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

qacc ColorMatrix gf_rect_staple_no_comm(const GaugeField& gf,
                                        const Coordinate& xl, const int mu)
// transpose the same way as gf.get_elem(xl, mu)
{
  ColorMatrix acc;
  set_zero(acc);
  for (int nu = -4; nu < 4; ++nu) {
    if (nu == mu or -nu - 1 == mu) {
      continue;
    }
    acc += gf_wilson_line_no_comm(
        gf, xl, make_array<int>(nu, nu, mu, -nu - 1, -nu - 1));
    acc += gf_wilson_line_no_comm(
        gf, xl, make_array<int>(nu, mu, mu, -nu - 1, -mu - 1));
    acc += gf_wilson_line_no_comm(
        gf, xl, make_array<int>(-mu - 1, nu, mu, mu, -nu - 1));
  }
  return acc;
}

qacc ColorMatrix gf_all_staple_no_comm(const GaugeField& gf,
                                       const GaugeAction& ga,
                                       const Coordinate& xl, const int mu)
// transpose the same way as gf.get_elem(xl, mu)
{
  ColorMatrix acc;
  set_zero(acc);
  const double c1 = ga.c1;
  acc += (ComplexD)(1.0 - 8.0 * c1) * gf_plaq_staple_no_comm(gf, xl, mu);
  if (c1 != 0.0) {
    acc += (ComplexD)c1 * gf_rect_staple_no_comm(gf, xl, mu);
  }
  return acc;
}

qacc ColorMatrix gf_force_site_no_comm(const GaugeField& gf,
                                       const GaugeAction& ga,
                                       const Coordinate& xl, const int mu)
{
  const double beta = ga.beta;
  const ColorMatrix ad_staple =
      matrix_adjoint(gf_all_staple_no_comm(gf, ga, xl, mu));
  const ColorMatrix force =
      (ComplexD)(-beta / 3.0) * (gf.get_elem(xl, mu) * ad_staple);
  return make_tr_less_anti_herm_matrix(force);
}

}  // namespace qlat
