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

struct GaugeMomentum : FieldM<ColorMatrix, 4> {
};

inline bool metropolis_accept(double& accept_prob, const double delta_h,
                              const int traj, const RngState& rs_)
// only compute at get_id_node() == 0
// broad_cast the result to all nodes
{
  TIMER_VERBOSE("metropolis_accept");
  double flag_d = 0.0;
  accept_prob = 0;
  if (get_id_node() == 0) {
    if (delta_h <= 0.0) {
      flag_d = 1.0;
      accept_prob = 1.0;
    } else {
      RngState rs = rs_;
      const double rand_num = u_rand_gen(rs, 0.0, 1.0);
      accept_prob = std::exp(-delta_h);
      if (rand_num <= accept_prob) {
        flag_d = 1.0;
      }
    }
  }
  bcast(get_data_one_elem(accept_prob));
  bcast(get_data_one_elem(flag_d));
  const bool flag = flag_d > 0.5;
  displayln_info(fname + ssprintf(": accept flag = %d with prob accept = "
                                  "%.1f%% deltaH = %.16f traj = %d",
                                  flag, accept_prob * 100, delta_h, traj));
  return flag;
}

inline void set_rand_gauge_momentum(GaugeMomentum& gm, const double sigma,
                                    const RngState& rs)
//  Creates a field of antihermitian 3x3 complex matrices with each complex
//  element drawn at random from a gaussian distribution with zero mean.
//  Hence the matrices are distributed according to
//
//  exp[- Tr(mat^2)/(2 sigma**2)]
{
  TIMER_VERBOSE("set_rand_gauge_momentum");
  set_g_rand_anti_hermitian_matrix_field(gm, rs, sigma);
}

inline double gm_hamilton_node(const GaugeMomentum& gm)
{
  TIMER("gm_hamilton_node");
  const Geometry geo = geo_reform(gm.geo());
  FieldM<double, 1> fd;
  fd.init(geo);
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<ColorMatrix> gm_v = gm.get_elems_const(xl);
    double s = 0.0;
    qassert(gm_v.size() == 4);
    for (int mu = 0; mu < 4; ++mu) {
      s += neg_half_tr_square(gm_v[mu]);
    }
    fd.get_elem(index) = s;
  });
  double sum = 0.0;
  for (long index = 0; index < geo.local_volume(); ++index) {
    sum += fd.get_elem(index);
  }
  return sum;
}

inline double gf_re_tr_plaq_no_comm(const GaugeField& gf, const Coordinate& xl,
                                    const int mu, const int nu)
{
  const ColorMatrix m =
      gf_wilson_line_no_comm(gf, xl, make_array<int>(mu, nu, -mu - 1, -nu - 1));
  return matrix_trace(m).real();
}

inline double gf_re_tr_rect_no_comm(const GaugeField& gf, const Coordinate& xl,
                                    const int mu, const int nu)
{
  const ColorMatrix m = gf_wilson_line_no_comm(
      gf, xl, make_array<int>(mu, mu, nu, -mu - 1, -mu - 1, -nu - 1));
  return matrix_trace(m).real();
}

inline double gf_sum_re_tr_plaq_node_no_comm(const GaugeField& gf)
// subtract the free field value
{
  TIMER("gf_sum_re_tr_plaq_node_no_comm");
  const Geometry geo = geo_reform(gf.geo());
  FieldM<double, 1> fd;
  fd.init(geo);
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    double s = 0.0;
    for (int mu = 0; mu < 3; ++mu) {
      for (int nu = mu + 1; nu < 4; ++nu) {
        s += gf_re_tr_plaq_no_comm(gf, xl, mu, nu) - 3.0;
      }
    }
    fd.get_elem(index) = s;
  });
  double sum = 0.0;
  for (long index = 0; index < geo.local_volume(); ++index) {
    sum += fd.get_elem(index);
  }
  return sum;
}

inline double gf_sum_re_tr_rect_node_no_comm(const GaugeField& gf)
// subtract the free field value
{
  TIMER("gf_sum_re_tr_rect_node_no_comm");
  const Geometry geo = geo_reform(gf.geo());
  FieldM<double, 1> fd;
  fd.init(geo);
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    double s = 0.0;
    for (int mu = 0; mu < 3; ++mu) {
      for (int nu = mu + 1; nu < 4; ++nu) {
        s += gf_re_tr_rect_no_comm(gf, xl, mu, nu) - 3.0;
        s += gf_re_tr_rect_no_comm(gf, xl, nu, mu) - 3.0;
      }
    }
    fd.get_elem(index) = s;
  });
  double sum = 0.0;
  for (long index = 0; index < geo.local_volume(); ++index) {
    sum += fd.get_elem(index);
  }
  return sum;
}

inline double gf_hamilton_node_no_comm(const GaugeField& gf,
                                       const GaugeAction& ga)
{
  TIMER("gf_hamilton_node_no_comm");
  const double beta = ga.beta;
  const double c1 = ga.c1;
  const double sum_plaq = gf_sum_re_tr_plaq_node_no_comm(gf);
  const double sum_rect = c1 == 0.0 ? 0.0 : gf_sum_re_tr_rect_node_no_comm(gf);
  return -beta / 3.0 * ((1.0 - 8.0 * c1) * sum_plaq + c1 * sum_rect);
}

inline double gf_hamilton_node(const GaugeField& gf, const GaugeAction& ga)
{
  TIMER("gf_hamilton_node");
  const Coordinate expand_left(0, 0, 0, 0);
  Coordinate expand_right(2, 2, 2, 2);
  if (ga.c1 == 0.0) {
    expand_right = Coordinate(1, 1, 1, 1);
  }
  const Geometry geo_ext = geo_resize(gf.geo(), expand_left, expand_right);
  GaugeField gf_ext;
  gf_ext.init(geo_ext);
  gf_ext = gf;
  const std::string tag_comm = ga.c1 == 0.0 ? "plaq" : "plaq+rect";
  const CommPlan& plan =
      get_comm_plan(set_marks_field_gf_hamilton, tag_comm, gf_ext.geo());
  refresh_expanded(gf_ext, plan);
  return gf_hamilton_node_no_comm(gf_ext, ga);
}

inline double gf_hamilton(const GaugeField& gf, const GaugeAction& ga)
// beta * ( (1-8*c1) n_plaq (1 - avg_plaq) + c1 n_rect (1 - avg_rect) )
// n_plaq = total_volume() * 6
// n_rect = total_volume() * 12
{
  TIMER("gf_hamilton");
  double energy = gf_hamilton_node(gf, ga);
  glb_sum(energy);
  return energy;
}

inline void gf_evolve(GaugeField& gf, const GaugeMomentum& gm,
                      const double step_size)
//  U(t+dt) = exp(i dt H) U(t)
{
  TIMER("gf_evolve");
  const Geometry& geo = gf.geo();
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> gf_v = gf.get_elems(xl);
    const Vector<ColorMatrix> gm_v = gm.get_elems_const(xl);
    qassert(gf_v.size() == 4);
    qassert(gm_v.size() == 4);
    for (int mu = 0; mu < 4; ++mu) {
      gf_v[mu] = matrix_evolve(gf_v[mu], gm_v[mu], step_size);
    }
  });
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
  acc += (Complex)(1.0 - 8.0 * c1) * gf_plaq_staple_no_comm(gf, xl, mu);
  if (c1 != 0.0) {
    acc += (Complex)c1 * gf_rect_staple_no_comm(gf, xl, mu);
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
      (Complex)(-beta / 3.0) * (gf.get_elem(xl, mu) * ad_staple);
  return make_tr_less_anti_herm_matrix(force);
}

inline void set_gm_force_no_comm(GaugeMomentum& gm_force, const GaugeField& gf,
                                 const GaugeAction& ga)
// gf need comm
{
  TIMER("set_gm_force_no_comm");
  const Geometry geo = geo_resize(gf.geo());
  gm_force.init(geo);
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> gm_force_v = gm_force.get_elems(xl);
    qassert(gm_force_v.size() == 4);
    for (int mu = 0; mu < 4; ++mu) {
      gm_force_v[mu] = gf_force_site_no_comm(gf, ga, xl, mu);
    }
  });
}

inline void set_gm_force(GaugeMomentum& gm_force, const GaugeField& gf,
                         const GaugeAction& ga)
{
  TIMER("set_gm_force");
  Coordinate expand_left(2, 2, 2, 2);
  Coordinate expand_right(2, 2, 2, 2);
  if (ga.c1 == 0.0) {
    expand_left = Coordinate(1, 1, 1, 1);
    expand_right = Coordinate(1, 1, 1, 1);
  }
  const Geometry geo_ext = geo_resize(gf.geo(), expand_left, expand_right);
  GaugeField gf_ext;
  gf_ext.init(geo_ext);
  gf_ext = gf;
  const std::string tag_comm = ga.c1 == 0.0 ? "plaq" : "plaq+rect";
  const CommPlan& plan =
      get_comm_plan(set_marks_field_gm_force, tag_comm, gf_ext.geo());
  refresh_expanded(gf_ext, plan);
  set_gm_force_no_comm(gm_force, gf_ext, ga);
}

}  // namespace qlat
