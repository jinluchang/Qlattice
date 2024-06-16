#include <qlat/hmc.h>

namespace qlat
{  //

bool metropolis_accept(double& accept_prob, const double delta_h,
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

void set_rand_gauge_momentum(GaugeMomentum& gm, const double sigma,
                             const RngState& rs)
//  Creates a field of antihermitian 3x3 complex matrices with each complex
//  element drawn at random from a gaussian distribution with zero mean.
//  Hence the matrices are distributed according to
//
// default sigma=1.0
//
//  exp[- Tr(mat^2)/(2 sigma**2)]
{
  TIMER_VERBOSE("set_rand_gauge_momentum");
  set_g_rand_anti_hermitian_matrix_field(gm, rs, sigma);
}

double gm_hamilton_node(const GaugeMomentum& gm)
{
  TIMER("gm_hamilton_node");
  const Geometry geo = geo_resize(gm.geo());
  FieldM<double, 1> fd;
  fd.init(geo);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = fd.geo();
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
  for (Long index = 0; index < geo.local_volume(); ++index) {
    sum += fd.get_elem(index);
  }
  return sum;
}

double gf_sum_re_tr_plaq_node_no_comm(const GaugeField& gf)
// subtract the free field value
{
  TIMER("gf_sum_re_tr_plaq_node_no_comm");
  const Geometry geo = geo_resize(gf.geo());
  FieldM<double, 1> fd;
  fd.init(geo);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = gf.geo();
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
  for (Long index = 0; index < geo.local_volume(); ++index) {
    sum += fd.get_elem(index);
  }
  return sum;
}

double gf_sum_re_tr_rect_node_no_comm(const GaugeField& gf)
// subtract the free field value
{
  TIMER("gf_sum_re_tr_rect_node_no_comm");
  const Geometry geo = geo_resize(gf.geo());
  FieldM<double, 1> fd;
  fd.init(geo);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = gf.geo();
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
  for (Long index = 0; index < geo.local_volume(); ++index) {
    sum += fd.get_elem(index);
  }
  return sum;
}

double gf_hamilton_node_no_comm(const GaugeField& gf, const GaugeAction& ga)
{
  TIMER("gf_hamilton_node_no_comm");
  const double beta = ga.beta;
  const double c1 = ga.c1;
  const double sum_plaq = gf_sum_re_tr_plaq_node_no_comm(gf);
  const double sum_rect = c1 == 0.0 ? 0.0 : gf_sum_re_tr_rect_node_no_comm(gf);
  return -beta / 3.0 * ((1.0 - 8.0 * c1) * sum_plaq + c1 * sum_rect);
}

double gf_hamilton_node(const GaugeField& gf, const GaugeAction& ga)
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
  const CommPlan& plan = get_comm_plan(set_marks_field_gf_hamilton, tag_comm,
                                       gf_ext.geo(), gf_ext.multiplicity);
  refresh_expanded(gf_ext, plan);
  return gf_hamilton_node_no_comm(gf_ext, ga);
}

double gf_hamilton(const GaugeField& gf, const GaugeAction& ga)
// beta * ( (1-8*c1) n_plaq (1 - avg_plaq) + c1 n_rect (1 - avg_rect) )
// n_plaq = total_volume() * 6
// n_rect = total_volume() * 12
{
  TIMER("gf_hamilton");
  double energy = gf_hamilton_node(gf, ga);
  glb_sum(energy);
  return energy;
}

void gf_evolve(GaugeField& gf, const GaugeMomentum& gm, const double step_size)
//  U(t+dt) = exp(pi * dt) U(t)
{
  TIMER("gf_evolve");
  qassert(gf.multiplicity == 4);
  qassert(gm.multiplicity == 4);
  qacc_for(index, gf.geo().local_volume(), {
    const Geometry& geo = gf.geo();
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

void gf_evolve_dual(GaugeField& gf, const GaugeMomentum& gm_dual,
                    const double step_size)
//  U(t+dt) = U(t) exp(-pi_dual * dt)
{
  TIMER("gf_evolve_dual");
  qassert(gf.multiplicity == 4);
  qassert(gm_dual.multiplicity == 4);
  qacc_for(index, gf.geo().local_volume(), {
    const Geometry& geo = gf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> gf_v = gf.get_elems(xl);
    const Vector<ColorMatrix> gm_v = gm_dual.get_elems_const(xl);
    qassert(gf_v.size() == 4);
    qassert(gm_v.size() == 4);
    for (int mu = 0; mu < 4; ++mu) {
      gf_v[mu] = matrix_evolve_dual(gf_v[mu], gm_v[mu], step_size);
    }
  });
}

void set_gm_force_no_comm(GaugeMomentum& gm_force, const GaugeField& gf,
                          const GaugeAction& ga)
// gf need comm
{
  TIMER("set_gm_force_no_comm");
  const Geometry geo = geo_resize(gf.geo());
  gm_force.init(geo);
  qassert(gf.multiplicity == 4);
  qassert(gm_force.multiplicity == 4);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = gm_force.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> gm_force_v = gm_force.get_elems(xl);
    qassert(gm_force_v.size() == 4);
    for (int mu = 0; mu < 4; ++mu) {
      gm_force_v[mu] = gf_force_site_no_comm(gf, ga, xl, mu);
    }
  });
}

void set_gm_force(GaugeMomentum& gm_force, const GaugeField& gf,
                  const GaugeAction& ga)
{
  TIMER("set_gm_force");
  qassert(gf.multiplicity == 4);
  Coordinate expand_left(2, 2, 2, 2);
  Coordinate expand_right(2, 2, 2, 2);
  if (ga.c1 == 0.0) {
    expand_left = Coordinate(1, 1, 1, 1);
    expand_right = Coordinate(1, 1, 1, 1);
  }
  const Geometry geo_ext = geo_resize(gf.geo(), expand_left, expand_right);
  GaugeField gf_ext;
  gf_ext.init(geo_ext);
  qassert(gf_ext.multiplicity == 4);
  gf_ext = gf;
  const std::string tag_comm = ga.c1 == 0.0 ? "plaq" : "plaq+rect";
  const CommPlan& plan = get_comm_plan(set_marks_field_gm_force, tag_comm,
                                       gf_ext.geo(), gf_ext.multiplicity);
  refresh_expanded(gf_ext, plan);
  set_gm_force_no_comm(gm_force, gf_ext, ga);
  qassert(gm_force.multiplicity == 4);
}

void set_gm_force_dual(GaugeMomentum& gm_force_dual, const GaugeField& gf,
                       const GaugeMomentum& gm_force)
{
  TIMER("set_gm_force_dual");
  gm_force_dual.init(gf.geo());
  qassert(gm_force_dual.multiplicity == 4);
  qassert(gf.multiplicity == 4);
  qassert(gm_force.multiplicity == 4);
  qacc_for(index, gf.geo().local_volume(), {
    const Geometry& geo = gf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> gm_dual_v = gm_force_dual.get_elems(xl);
    const Vector<ColorMatrix> gf_v = gf.get_elems_const(xl);
    const Vector<ColorMatrix> gm_v = gm_force.get_elems_const(xl);
    qassert(gm_dual_v.size() == 4);
    qassert(gf_v.size() == 4);
    qassert(gm_v.size() == 4);
    for (int mu = 0; mu < 4; ++mu) {
      gm_dual_v[mu] = -matrix_adjoint(gf_v[mu]) * gm_v[mu] * gf_v[mu];
    }
  });
}

}  // namespace qlat
