#include <qlat/hmc.h>
#include <qlat/qcd-acc.h>

namespace qlat
{  //

static qacc RealD gf_re_tr_plaq_no_comm(const GaugeField& gf,
                                         const Coordinate& xl, const Int mu,
                                         const Int nu)
{
  const ColorMatrix m =
      gf_wilson_line_no_comm(gf, xl, make_array<int>(mu, nu, -mu - 1, -nu - 1));
  return matrix_trace(m).real();
}

static qacc RealD gf_re_tr_rect_no_comm(const GaugeField& gf,
                                         const Coordinate& xl, const Int mu,
                                         const Int nu)
{
  const ColorMatrix m = gf_wilson_line_no_comm(
      gf, xl, make_array<int>(mu, mu, nu, -mu - 1, -mu - 1, -nu - 1));
  return matrix_trace(m).real();
}

static qacc ColorMatrix gf_rect_staple_no_comm(const GaugeField& gf,
                                               const Coordinate& xl,
                                               const Int mu)
// transpose the same way as gf.get_elem(xl, mu)
{
  ColorMatrix acc;
  set_zero(acc);
  for (Int nu = -4; nu < 4; ++nu) {
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

static qacc ColorMatrix gf_all_staple_no_comm(const GaugeField& gf,
                                              const GaugeAction& ga,
                                              const Coordinate& xl,
                                              const Int mu)
// transpose the same way as gf.get_elem(xl, mu)
{
  ColorMatrix acc;
  set_zero(acc);
  const RealD c1 = ga.c1;
  acc += (ComplexD)(1.0 - 8.0 * c1) * gf_plaq_staple_no_comm(gf, xl, mu);
  if (c1 != 0.0) {
    acc += (ComplexD)c1 * gf_rect_staple_no_comm(gf, xl, mu);
  }
  return acc;
}

static qacc ColorMatrix gf_force_site_no_comm(const GaugeField& gf,
                                              const GaugeAction& ga,
                                              const Coordinate& xl,
                                              const Int mu)
{
  const RealD beta = ga.beta;
  const ColorMatrix ad_staple =
      matrix_adjoint(gf_all_staple_no_comm(gf, ga, xl, mu));
  const ColorMatrix force =
      (ComplexD)(-beta / 3.0) * (gf.get_elem(xl, mu) * ad_staple);
  return make_tr_less_anti_herm_matrix(force);
}

bool metropolis_accept(RealD& accept_prob, const RealD delta_h,
                       const Int traj, const RngState& rs_)
// only compute at get_id_node() == 0
// broad_cast the result to all nodes
{
  TIMER_VERBOSE("metropolis_accept");
  RealD flag_d = 0.0;
  accept_prob = 0;
  if (get_id_node() == 0) {
    if (delta_h <= 0.0) {
      flag_d = 1.0;
      accept_prob = 1.0;
    } else {
      RngState rs = rs_;
      const RealD rand_num = u_rand_gen(rs, 0.0, 1.0);
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

void set_rand_gauge_momentum(GaugeMomentum& gm, const RealD sigma,
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

void set_rand_gauge_momentum(GaugeMomentum& gm, const Field<RealD>& mf,
                             const RngState& rs)
//  Creates a field of antihermitian 3x3 complex matrices with each complex
//  element drawn at random from a gaussian distribution with zero mean.
//  Hence the matrices are distributed according to
//
// default mf=1.0
//
//  exp[- Tr(mat^2)/(2 mf)]
{
  TIMER_VERBOSE("set_rand_gauge_momentum");
  set_g_rand_anti_hermitian_matrix_field(gm, rs, 1.0);
  qassert(gm.multiplicity == 4);
  qassert(mf.multiplicity == 4);
  const Geometry& geo = gm.geo();
  qassert(check_matching_geo(geo, mf.geo()));
  qacc_for(index, geo.local_volume(), {
    const Vector<RealD> vm = mf.get_elems_const(index);
    Vector<ColorMatrix> v = gm.get_elems(index);
    const RealD inf = std::numeric_limits<RealD>::infinity();
    for (Int m = 0; m < 4; ++m) {
      if (vm[m] != inf) {
        v[m] *= std::sqrt(vm[m]);
      } else {
        set_zero(v[m]);
      }
    }
  });
}

RealD gm_hamilton_node(const GaugeMomentum& gm)
{
  TIMER("gm_hamilton_node");
  const Geometry geo = geo_resize(gm.geo());
  FieldM<RealD, 1> fd;
  fd.init(geo);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = fd.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<ColorMatrix> gm_v = gm.get_elems_const(xl);
    RealD s = 0.0;
    qassert(gm_v.size() == 4);
    for (Int mu = 0; mu < 4; ++mu) {
      s += neg_half_tr_square(gm_v[mu]);
    }
    fd.get_elem(index) = s;
  });
  RealD sum = 0.0;
  for (Long index = 0; index < geo.local_volume(); ++index) {
    sum += fd.get_elem(index);
  }
  return sum;
}

RealD gm_hamilton_node(const GaugeMomentum& gm, const Field<RealD>& mf)
{
  TIMER("gm_hamilton_node");
  const Geometry geo = geo_resize(gm.geo());
  qassert(check_matching_geo(geo, mf.geo()));
  qassert(gm.multiplicity == 4);
  qassert(mf.multiplicity == 4);
  FieldM<RealD, 1> fd;
  fd.init(geo);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = fd.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<ColorMatrix> gm_v = gm.get_elems_const(xl);
    const Vector<RealD> mf_v = mf.get_elems_const(xl);
    RealD s = 0.0;
    qassert(gm_v.size() == 4);
    for (Int mu = 0; mu < 4; ++mu) {
      s += neg_half_tr_square(gm_v[mu]) / mf_v[mu];
    }
    fd.get_elem(index) = s;
  });
  RealD sum = 0.0;
  for (Long index = 0; index < geo.local_volume(); ++index) {
    sum += fd.get_elem(index);
  }
  return sum;
}

RealD gf_sum_re_tr_plaq_node_no_comm(const GaugeField& gf)
// subtract the free field value
{
  TIMER("gf_sum_re_tr_plaq_node_no_comm");
  const Geometry geo = geo_resize(gf.geo());
  FieldM<RealD, 1> fd;
  fd.init(geo);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = gf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    RealD s = 0.0;
    for (Int mu = 0; mu < 3; ++mu) {
      for (Int nu = mu + 1; nu < 4; ++nu) {
        s += gf_re_tr_plaq_no_comm(gf, xl, mu, nu) - 3.0;
      }
    }
    fd.get_elem(index) = s;
  });
  RealD sum = 0.0;
  for (Long index = 0; index < geo.local_volume(); ++index) {
    sum += fd.get_elem(index);
  }
  return sum;
}

RealD gf_sum_re_tr_rect_node_no_comm(const GaugeField& gf)
// subtract the free field value
{
  TIMER("gf_sum_re_tr_rect_node_no_comm");
  const Geometry geo = geo_resize(gf.geo());
  FieldM<RealD, 1> fd;
  fd.init(geo);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = gf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    RealD s = 0.0;
    for (Int mu = 0; mu < 3; ++mu) {
      for (Int nu = mu + 1; nu < 4; ++nu) {
        s += gf_re_tr_rect_no_comm(gf, xl, mu, nu) - 3.0;
        s += gf_re_tr_rect_no_comm(gf, xl, nu, mu) - 3.0;
      }
    }
    fd.get_elem(index) = s;
  });
  RealD sum = 0.0;
  for (Long index = 0; index < geo.local_volume(); ++index) {
    sum += fd.get_elem(index);
  }
  return sum;
}

RealD gf_hamilton_node_no_comm(const GaugeField& gf, const GaugeAction& ga)
// number of plaq: 6 * number of site
// number of rect: 12 * number of site
/*
  \ba
  S_\text{gauge}
  =&
  \frac{ \beta }{ 3 }
  \Big[
  (1-8 c_1) \sum_P \mathrm{Re}\mathrm{Tr} (1 - U_P)
  +
  c_1 \sum_R \mathrm{Re}\mathrm{Tr} (1 - U_R)
  \Big]
  \\
  \ea
*/
{
  TIMER("gf_hamilton_node_no_comm");
  const RealD beta = ga.beta;
  const RealD c1 = ga.c1;
  const RealD sum_plaq = gf_sum_re_tr_plaq_node_no_comm(gf);
  const RealD sum_rect = c1 == 0.0 ? 0.0 : gf_sum_re_tr_rect_node_no_comm(gf);
  return -beta / 3.0 * ((1.0 - 8.0 * c1) * sum_plaq + c1 * sum_rect);
}

RealD gf_hamilton_node(const GaugeField& gf, const GaugeAction& ga)
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
  QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
  const CommPlan& plan = get_comm_plan(set_marks_field_gf_hamilton, tag_comm,
                                       gf_ext.geo(), gf_ext.multiplicity);
  QLAT_DIAGNOSTIC_POP;
  refresh_expanded(gf_ext, plan);
  return gf_hamilton_node_no_comm(gf_ext, ga);
}

RealD gf_hamilton(const GaugeField& gf, const GaugeAction& ga)
// beta * ( (1-8*c1) n_plaq (1 - avg_plaq) + c1 n_rect (1 - avg_rect) )
// n_plaq = total_volume() * 6
// n_rect = total_volume() * 12
{
  TIMER("gf_hamilton");
  RealD energy = gf_hamilton_node(gf, ga);
  glb_sum(energy);
  return energy;
}

void gf_evolve(GaugeField& gf, const GaugeMomentum& gm, const RealD step_size)
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
    for (Int mu = 0; mu < 4; ++mu) {
      gf_v[mu] = matrix_evolve(gf_v[mu], gm_v[mu], step_size);
    }
  });
}

void gf_evolve_dual(GaugeField& gf, const GaugeMomentum& gm_dual,
                    const RealD step_size)
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
    for (Int mu = 0; mu < 4; ++mu) {
      gf_v[mu] = matrix_evolve_dual(gf_v[mu], gm_v[mu], step_size);
    }
  });
}

void gf_evolve(GaugeField& gf, const GaugeMomentum& gm, const Field<RealD>& mf,
               const RealD step_size)
//  U(t+dt) = exp(pi * dt / mass) U(t)
{
  TIMER("gf_evolve");
  qassert(gf.multiplicity == 4);
  qassert(gm.multiplicity == 4);
  qassert(mf.multiplicity == 4);
  qacc_for(index, gf.geo().local_volume(), {
    const Geometry& geo = gf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> gf_v = gf.get_elems(xl);
    const Vector<ColorMatrix> gm_v = gm.get_elems_const(xl);
    const Vector<RealD> mf_v = mf.get_elems_const(xl);
    qassert(gf_v.size() == 4);
    qassert(gm_v.size() == 4);
    qassert(mf_v.size() == 4);
    for (Int mu = 0; mu < 4; ++mu) {
      const RealD dt_m = step_size / mf_v[mu];
      gf_v[mu] = matrix_evolve(gf_v[mu], gm_v[mu], dt_m);
    }
  });
}

void gf_evolve_dual(GaugeField& gf, const GaugeMomentum& gm_dual,
                    const Field<RealD>& mf_dual, const RealD step_size)
//  U(t+dt) = U(t) exp(-pi_dual * dt / mass)
{
  TIMER("gf_evolve_dual");
  qassert(gf.multiplicity == 4);
  qassert(gm_dual.multiplicity == 4);
  qassert(mf_dual.multiplicity == 4);
  qacc_for(index, gf.geo().local_volume(), {
    const Geometry& geo = gf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> gf_v = gf.get_elems(xl);
    const Vector<ColorMatrix> gm_v = gm_dual.get_elems_const(xl);
    const Vector<RealD> mf_v = mf_dual.get_elems_const(xl);
    qassert(gf_v.size() == 4);
    qassert(gm_v.size() == 4);
    qassert(mf_v.size() == 4);
    for (Int mu = 0; mu < 4; ++mu) {
      const RealD dt_m = step_size / mf_v[mu];
      gf_v[mu] = matrix_evolve_dual(gf_v[mu], gm_v[mu], dt_m);
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
    for (Int mu = 0; mu < 4; ++mu) {
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
  QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
  const CommPlan& plan = get_comm_plan(set_marks_field_gm_force, tag_comm,
                                       gf_ext.geo(), gf_ext.multiplicity);
  QLAT_DIAGNOSTIC_POP;
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
    for (Int mu = 0; mu < 4; ++mu) {
      gm_dual_v[mu] = -matrix_adjoint(gf_v[mu]) * gm_v[mu] * gf_v[mu];
    }
  });
}

RealD project_gauge_transform(GaugeMomentum& gm, GaugeMomentum& gm_dual,
                              const Field<RealD>& mf,
                              const Field<RealD>& mf_dual)
// return qnorm of the change / geo.total_volume()
// (up to an overall coefficient)
{
  TIMER("project_gauge_transform");
  const Geometry& geo = gm.geo();
  qassert(check_matching_geo(geo, gm_dual.geo()));
  qassert(check_matching_geo(geo, mf.geo()));
  qassert(check_matching_geo(geo, mf_dual.geo()));
  qassert(gm.multiplicity == 4);
  qassert(gm_dual.multiplicity == 4);
  qassert(mf.multiplicity == 4);
  qassert(mf_dual.multiplicity == 4);
  const Coordinate expand_left = Coordinate(1, 1, 1, 1);
  const Geometry geo_ext = geo_resize(geo, expand_left, Coordinate());
  GaugeMomentum gm2_ext;
  gm2_ext.init(geo_ext, 4);
  gm2_ext = gm_dual;
  Field<RealD> mf2_ext;
  mf2_ext.init(geo_ext, 4);
  mf2_ext = mf_dual;
  refresh_expanded_1(gm2_ext);
  refresh_expanded_1(mf2_ext);
  const Coordinate expand_right = Coordinate(1, 1, 1, 1);
  const Geometry geo_extr = geo_resize(geo, Coordinate(), expand_right);
  Field<ColorMatrix> gm_ag_ext;
  gm_ag_ext.init(geo_extr, 1);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = gm.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<ColorMatrix> v1 = gm.get_elems_const(xl);
    const Vector<RealD> vm1 = mf.get_elems_const(xl);
    const RealD inf = std::numeric_limits<RealD>::infinity();
    ColorMatrix v_sum;
    set_zero(v_sum);
    for (Int m = 0; m < 4; ++m) {
      const Coordinate xl_m = coordinate_shifts(xl, -m - 1);
      const RealD m1 = vm1[m];
      const RealD m2m = mf2_ext.get_elem(xl_m, m);
      const ColorMatrix& c1 = v1[m];
      const ColorMatrix& c2m = gm2_ext.get_elem(xl_m, m);
      if (m1 != inf) {
        v_sum += c1 * (1.0 / std::sqrt(m1));
      }
      if (m2m != inf) {
        v_sum += c2m * (1.0 / std::sqrt(m2m));
      }
    }
    ColorMatrix& ag = gm_ag_ext.get_elem(xl);
    ag = (1.0 / 8.0) * v_sum;
  });
  refresh_expanded_1(gm_ag_ext);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = gm.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v1 = gm.get_elems(index);
    Vector<ColorMatrix> v2 = gm_dual.get_elems(index);
    const Vector<RealD> vm1 = mf.get_elems_const(xl);
    const Vector<RealD> vm2 = mf_dual.get_elems_const(xl);
    const ColorMatrix& ag1 = gm_ag_ext.get_elem(xl);
    const RealD inf = std::numeric_limits<RealD>::infinity();
    for (Int m = 0; m < 4; ++m) {
      const Coordinate xl_p = coordinate_shifts(xl, m);
      if (vm1[m] != inf) {
        v1[m] -= std::sqrt(vm1[m]) * ag1;
      }
      if (vm2[m] != inf) {
        v2[m] -= std::sqrt(vm2[m]) * gm_ag_ext.get_elem(xl_p);
      }
    }
  });
  return qnorm(gm_ag_ext) / geo.total_volume();
}

void dot_gauge_momentum(Field<RealD>& f, const GaugeMomentum& gm1,
                        const GaugeMomentum& gm2)
{
  TIMER("dot_gauge_momentum(f,gm1,gm2)");
  qassert(check_matching_geo(gm1.geo(), gm2.geo()));
  qassert(gm1.multiplicity == 4);
  qassert(gm2.multiplicity == 4);
  f.init(gm1.geo(), 4);
  qacc_for(index, f.geo().local_volume(), {
    // const Geometry& geo = f1.geo();
    const Vector<ColorMatrix> v1 = gm1.get_elems_const(index);
    const Vector<ColorMatrix> v2 = gm2.get_elems_const(index);
    Vector<RealD> v = f.get_elems(index);
    for (Int m = 0; m < v.size(); ++m) {
      v[m] = qnorm(v1[m], v2[m]);
    }
  });
}

}  // namespace qlat
