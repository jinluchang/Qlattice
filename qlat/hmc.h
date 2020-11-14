#pragma once

#include <qlat/qcd.h>

#include <cmath>
#include <sstream>
#include <string>

namespace qlat
{  //

struct GaugeMomentum : FieldM<ColorMatrix, 4> {
  virtual const std::string& cname()
  {
    static const std::string s = "GaugeMomentum";
    return s;
  }
};

inline bool metropolis_accept(double& accept_prob, const double delta_h,
                              const RngState& rs_)
// only compute at get_id_node() == 0
// broad_cast the result to all nodes
{
  TIMER_VERBOSE("metropolis_accept");
  double flag = 0.0;
  accept_prob = 0;
  if (get_id_node() == 0) {
    if (delta_h <= 0.0) {
      flag = 1.0;
      accept_prob = 1.0;
    } else {
      RngState rs = rs_;
      const double rand_num = u_rand_gen(rs, 0.0, 1.0);
      accept_prob = std::exp(-delta_h);
      if (rand_num <= accept_prob) {
        flag = 1.0;
      }
    }
  }
  bcast(get_data_one_elem(accept_prob));
  bcast(get_data_one_elem(flag));
  return flag > 0.5;
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
  const Geometry geo = geo_reform(gm.geo);
  FieldM<double, 1> fd;
  fd.init(geo);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<ColorMatrix> gm_v = gm.get_elems_const(xl);
    double s = 0.0;
    qassert(gm_v.size() == 4);
    for (int mu = 0; mu < 4; ++mu) {
      s += neg_half_tr_square(gm_v[mu]);
    }
    fd.get_elem(index) = s;
  }
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
{
  TIMER("gf_sum_re_tr_plaq_node_no_comm");
  const Geometry geo = geo_reform(gf.geo);
  FieldM<double, 1> fd;
  fd.init(geo);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    double s = 0.0;
    for (int mu = 0; mu < 3; ++mu) {
      for (int nu = mu + 1; nu < 4; ++nu) {
        s += gf_re_tr_plaq_no_comm(gf, xl, mu, nu);
      }
    }
    fd.get_elem(index) = s;
  }
  double sum = 0.0;
  for (long index = 0; index < geo.local_volume(); ++index) {
    sum += fd.get_elem(index);
  }
  return sum;
}

inline double gf_sum_re_tr_rect_node_no_comm(const GaugeField& gf)
{
  TIMER("gf_sum_re_tr_rect_node_no_comm");
  const Geometry geo = geo_reform(gf.geo);
  FieldM<double, 1> fd;
  fd.init(geo);
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    double s = 0.0;
    for (int mu = 0; mu < 3; ++mu) {
      for (int nu = mu + 1; nu < 4; ++nu) {
        s += gf_re_tr_rect_no_comm(gf, xl, mu, nu);
        s += gf_re_tr_rect_no_comm(gf, xl, nu, mu);
      }
    }
    fd.get_elem(index) = s;
  }
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
  const double sum_plaq = gf_sum_re_tr_plaq_node_no_comm(gf);
  const double sum_rect = gf_sum_re_tr_rect_node_no_comm(gf);
  const double beta = ga.beta;
  const double c1 = ga.c1;
  return -beta / 3.0 * ((1.0 - 8.0 * c1) * sum_plaq + c1 * sum_rect);
}

inline void gf_evolve(GaugeField& gf, const GaugeMomentum& gm, const double step_size)
//  U(t+dt) = exp(i dt H) U(t)
{
  TIMER("gf_evolve");
  const Geometry& geo = gf.geo;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> gf_v = gf.get_elems(xl);
    const Vector<ColorMatrix> gm_v = gm.get_elems_const(xl);
    qassert(gf_v.size() == 4);
    qassert(gm_v.size() == 4);
    for (int mu = 0; mu < 4; ++mu) {
      gf_v[mu] = matrix_evolve(gf_v[mu], gm_v[mu], step_size);
    }
  }
}

}  // namespace qlat
