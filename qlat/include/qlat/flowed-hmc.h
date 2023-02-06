#pragma once

#include <qlat/field-expand.h>
#include <qlat/hmc.h>

namespace qlat
{  //

// mask -> mask
// flow_size -> flow_type

// Masking scheme: mask

// get_mask_node_site: node_site -> mask_node_site
//
// get_num_mask : flow_type -> num_mask
//
// coordinate_from_mask_coordinate : mxg, mask, flow_type -> xg
//
// mask_from_coordinate : xg, flow_type -> mask
//
// mask_coordinate_from_coordinate : xg, flow_type -> mxg

struct FlowStepInfo {
  int mask;
  int mu;
  double epsilon;
  int flow_size;
  //
  FlowStepInfo() {}
  FlowStepInfo(const int mask_, const int mu_, const double epsilon_,
               const int flow_size_ = 1)
  {
    mask = mask_;
    mu = mu_;
    epsilon = epsilon_;
    flow_size = flow_size_;
  }
};

struct FlowInfo {
  std::vector<FlowStepInfo> v;
};

inline std::string show(const FlowInfo& fi);

inline void gf_flow(GaugeField& gf, const GaugeField& gf0, const FlowInfo& fi);

inline void gf_flow_inv(GaugeField& gf, const GaugeField& gf1,
                        const FlowInfo& fi);

inline double gf_hamilton_flowed_node(const GaugeField& gf0,
                                      const GaugeAction& ga,
                                      const FlowInfo& fi);

inline void set_gm_force_flowed(GaugeMomentum& gm_force, const GaugeField& gf0,
                                const GaugeAction& ga, const FlowInfo& fi);

inline void set_gm_force_flowed_no_det(GaugeMomentum& gm_force,
                                       const GaugeMomentum& gm_force_pre,
                                       const GaugeField& gf0,
                                       const FlowInfo& fi);

inline double gf_flow_and_ln_det_node(GaugeField& gf, const GaugeField& gf0,
                                      const FlowInfo& fi);

inline void set_flowed_gauge_fields(std::vector<GaugeField>& gf_ext_vec,
                                    const GaugeField& gf0, const FlowInfo& fi);

inline void set_gm_force_propagated_det_from_flow(
    GaugeMomentum& gm_force, const GaugeMomentum& gm_force_pre,
    const std::vector<GaugeField>& gf_ext_vec, const FlowInfo& fi);

inline void set_gm_force_propagated_and_gm_force_det_from_flow_step(
    GaugeMomentum& gm_force_propagated, GaugeMomentum& gm_force_det,
    const GaugeMomentum& gm_force_pre, const GaugeField& gf0_ext,
    const FlowStepInfo& fsi);

inline void set_gm_force_propagated_no_det_from_flow(
    GaugeMomentum& gm_force, const GaugeMomentum& gm_force_pre,
    const std::vector<GaugeField>& gf_ext_vec, const FlowInfo& fi);

inline void set_gm_force_propagated_from_flow_step(
    GaugeMomentum& gm_force_propagated, const GaugeMomentum& gm_force_pre,
    const GaugeField& gf0_ext, const FlowStepInfo& fsi);

// -------------------------------------------------------------------------

inline std::string show(const FlowInfo& fi)
{
  std::ostringstream out;
  for (int i = 0; i < (int)fi.v.size(); ++i) {
    const FlowStepInfo& fsi = fi.v[i];
    out << ssprintf("fi.v[%d]: mask=%d, mu=%d, epsilon=%.4f, flow_size=%d.", i,
                    fsi.mask, fsi.mu, fsi.epsilon, fsi.flow_size);
    if (i != (int)fi.v.size() - 1) {
      out << std::endl;
    }
  }
  return out.str();
}

inline FlowInfo mk_flow_info_step(const RngState& rs, const double epsilon)
{
  FlowInfo fi;
  fi.v.push_back(FlowStepInfo(2, 3, epsilon));
  fi.v.push_back(FlowStepInfo(1, 3, epsilon));
  fi.v.push_back(FlowStepInfo(2, 2, epsilon));
  fi.v.push_back(FlowStepInfo(1, 2, epsilon));
  fi.v.push_back(FlowStepInfo(2, 1, epsilon));
  fi.v.push_back(FlowStepInfo(1, 1, epsilon));
  fi.v.push_back(FlowStepInfo(2, 0, epsilon));
  fi.v.push_back(FlowStepInfo(1, 0, epsilon));
  random_permute(fi.v, rs);
  return fi;
}

inline FlowInfo mk_flow_info_step(const RngState& rs, const double epsilon,
                                  const double epsilon2)
{
  std::vector<std::vector<FlowStepInfo> > fi1s(4), fi2s(4);
  if (epsilon != 0.0) {
    fi1s[3].push_back(FlowStepInfo(2, 3, epsilon));
    fi1s[3].push_back(FlowStepInfo(1, 3, epsilon));
    fi1s[2].push_back(FlowStepInfo(2, 2, epsilon));
    fi1s[2].push_back(FlowStepInfo(1, 2, epsilon));
    fi1s[1].push_back(FlowStepInfo(2, 1, epsilon));
    fi1s[1].push_back(FlowStepInfo(1, 1, epsilon));
    fi1s[0].push_back(FlowStepInfo(2, 0, epsilon));
    fi1s[0].push_back(FlowStepInfo(1, 0, epsilon));
  }
  if (epsilon2 != 0.0) {
    fi2s[3].push_back(FlowStepInfo(2, 3, epsilon2, 2));
    fi2s[3].push_back(FlowStepInfo(1, 3, epsilon2, 2));
    fi2s[2].push_back(FlowStepInfo(2, 2, epsilon2, 2));
    fi2s[2].push_back(FlowStepInfo(1, 2, epsilon2, 2));
    fi2s[1].push_back(FlowStepInfo(2, 1, epsilon2, 2));
    fi2s[1].push_back(FlowStepInfo(1, 1, epsilon2, 2));
    fi2s[0].push_back(FlowStepInfo(2, 0, epsilon2, 2));
    fi2s[0].push_back(FlowStepInfo(1, 0, epsilon2, 2));
  }
  for (int i = 0; i < 4; ++i) {
    random_permute(fi1s[i], rs.split(ssprintf("fi1s-%d", i)));
    random_permute(fi2s[i], rs.split(ssprintf("fi2s-%d", i)));
    vector_append(fi1s[i], fi2s[i]);
  }
  random_permute(fi1s, rs);
  FlowInfo fi;
  fi.v = vector_concat(fi1s);
  return fi;
}

qacc int is_same_link(const Coordinate& xl, const int mu,
                      const Coordinate& xl_ref, const int mu_ref)
// return 0: not the same link
// return 1: the same link, same direction
// return -1: the same link, opposite direction
{
  if (0 <= mu) {
    return mu == mu_ref and xl == xl_ref ? 1 : 0;
  } else {
    return -mu - 1 == mu_ref and coordinate_shifts(xl, mu) == xl_ref ? -1 : 0;
  }
}

// qacc int mask_block2_node_from_geo(const Geometry& geo)
// {
//   const Coordinate& coor_node = geo.geon.coor_node;
//   const Coordinate& node_site = geo.node_site;
//   const Coordinate xl = coor_node * node_site;
//   return 2 - (xl[0] / 2 + xl[1] / 2 + xl[2] / 2 + xl[3] / 2) % 2;
// }

qacc int mask_block2_from_coordinate(const Coordinate& xl,
                                     const int mask_block2_node)
{
  return 2 - (mask_block2_node + (xl[0] + 1024) / 2 + (xl[1] + 1024) / 2 +
              (xl[2] + 1024) / 2 + (xl[3] + 1024) / 2) %
                 2;
}

qacc int mask_from_coordinate(const Coordinate& xg, const int flow_size)
{
  if (1 == flow_size) {
    return eo_from_coordinate(xg);
  } else if (2 == flow_size) {
    return mask_block2_from_coordinate(xg, 0);
  } else {
    qassert(false);
    return 0;
  }
}

qacc int multiplicity_flow_hmc_plaq(const bool is_same_mask_as_flow)
{
  if (is_same_mask_as_flow) {
    return 1 + 3 * 2;
  } else {
    return 6 + 3 * 2;
  }
}

qacc int multiplicity_flow_hmc_srect(const bool is_same_mask_as_flow)
// srect stand for special rectangular plaq
// (the short edge of the rectangule is the flowed link)
{
  if (is_same_mask_as_flow) {
    return 1 + 3 * 4;
  } else {
    return 6 + 3 * 4;
  }
}

qacc int multiplicity_flow_hmc_max(const int flow_size)
{
  if (1 == flow_size) {
    return multiplicity_flow_hmc_plaq(false);
  } else if (2 == flow_size) {
    return multiplicity_flow_hmc_srect(false);
  } else {
    qassert(false);
    return 0;
  }
}

qacc int multiplicity_flow_hmc(const bool is_same_mask_as_flow,
                               const int flow_size)
{
  if (1 == flow_size) {
    return multiplicity_flow_hmc_plaq(is_same_mask_as_flow);
  } else if (2 == flow_size) {
    return multiplicity_flow_hmc_srect(is_same_mask_as_flow);
  } else {
    qassert(false);
    return 0;
  }
}

inline const vector<long>& get_flowed_hmc_indices_mask_flow_size(
    const Geometry& geo, const int mask, const int flow_size)
{
  static Cache<std::string, vector<long> > cache("flowed_hmc_indices_cache", 8,
                                                 2);
  const std::string key =
      ssprintf("%s-%d-%d", show(geo.node_site).c_str(), mask, flow_size);
  vector<long>& vec = cache[key];
  if (vec.size() == 0) {
    long count = 0;
    qfor(index, geo.local_volume(), {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      const int mask_xl = mask_from_coordinate(xg, flow_size);
      if (mask_xl == mask) {
        count += 1;
      }
    });
    vec.resize(count);
    count = 0;
    qfor(index, geo.local_volume(), {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Coordinate xg = geo.coordinate_g_from_l(xl);
      const int mask_xl = mask_from_coordinate(xg, flow_size);
      if (mask_xl == mask) {
        vec[count] = index;
        count += 1;
      }
    });
  }
  return cache[key];
}

qacc void set_xl_nu_from_mask_mu_yl_m_plaq(Coordinate& xl, int& nu,
                                           const int mask, const int mu,
                                           const Coordinate& yl, const int m,
                                           const Geometry& geo,
                                           const int flow_type)
// mask_from_coordinate(xg, 1) = mask
// mask, mu are the flow parameters
//
// xl, mu is the link to be flowed
// yl, nu is the link to calculate derivative
//
// given yl, m is the index of the pair: (xl, mu) (yl, nu)
//
// m: mu(1 xl), (mu+1)%4(2 xl), (mu+2)%4(2 xl), (mu+3)%4(2 xl)
// m: mu(6 xl), (mu+1)%4(2 xl), (mu+2)%4(2 xl), (mu+3)%4(2 xl)
{
  qassert(flow_type == 1);
  const Coordinate yg = geo.coordinate_g_from_l(yl);
  const int mask_yl = mask_from_coordinate(yg, flow_type);
  if (mask == mask_yl) {
    const int nu_size = 1;
    if (0 <= m and m < nu_size) {
      nu = mu;
      xl = yl;
    } else {
      qassert(nu_size <= m and m < nu_size + 3 * 2);
      const int k = (m - nu_size) / 2;
      const int l = (m - nu_size) % 2;
      nu = mod(mu + k + 1, 4);
      if (0 == l) {
        xl = yl;
      } else if (1 == l) {
        xl = coordinate_shifts(yl, nu, -mu - 1);
      } else {
        qassert(false);
      }
    }
  } else if (mask == 3 - mask_yl) {
    const int nu_size = 6;
    if (0 <= m and m < nu_size) {
      nu = mu;
      const int k = m / 2;
      const int l = m % 2;
      int dir = mod(nu + k + 1, 4);
      if (0 == l) {
        xl = coordinate_shifts(yl, dir);
      } else if (1 == l) {
        xl = coordinate_shifts(yl, -dir - 1);
      } else {
        qassert(false);
      }
    } else {
      qassert(nu_size <= m and m < nu_size + 3 * 2);
      const int k = (m - nu_size) / 2;
      const int l = (m - nu_size) % 2;
      nu = mod(mu + k + 1, 4);
      if (0 == l) {
        xl = coordinate_shifts(yl, -mu - 1);
      } else if (1 == l) {
        xl = coordinate_shifts(yl, nu);
      } else {
        qassert(false);
      }
    }
  } else {
    qassert(false);
  }
  const Coordinate xg = geo.coordinate_g_from_l(xl);
  qassert(mask == mask_from_coordinate(xg, flow_type));
}

qacc void set_xg_nu_from_mask_mu_yg_m_srect_nu_neq_mu(
    Coordinate& xg, int& nu, const int mask_yl, const int nu_size,
    const int mask, const int mu, const Coordinate& yg, const int m,
    const int flow_type)
{
  (void)mask_yl;
  qassert(flow_type == 2);
  qassert(nu_size <= m and m < nu_size + 3 * 4);
  const int k = (m - nu_size) / 4;
  const int l = (m - nu_size) % 4;
  nu = mod(mu + k + 1, 4);
  if (0 == l) {
    xg = coordinate_shifts(yg, -nu - 1);
    while (mask != mask_from_coordinate(xg, flow_type)) {
      xg = coordinate_shifts(xg, nu);
    }
  } else if (1 == l) {
    xg = coordinate_shifts(yg, nu, nu);
    while (mask != mask_from_coordinate(xg, flow_type)) {
      xg = coordinate_shifts(xg, -nu - 1);
    }
  } else if (2 == l) {
    xg = coordinate_shifts(yg, -mu - 1, -nu - 1);
    while (mask != mask_from_coordinate(xg, flow_type)) {
      xg = coordinate_shifts(xg, nu);
    }
  } else if (3 == l) {
    xg = coordinate_shifts(yg, -mu - 1, nu, nu);
    while (mask != mask_from_coordinate(xg, flow_type)) {
      xg = coordinate_shifts(xg, -nu - 1);
    }
  } else {
    qassert(false);
  }
}

qacc void set_xl_nu_from_mask_mu_yl_m_srect(Coordinate& xl, int& nu,
                                            const int mask, const int mu,
                                            const Coordinate& yl, const int m,
                                            const Geometry& geo,
                                            const int flow_type)
// mask_from_coordinate(xg, 2) = mask
// mask, mu are the flow parameters
//
// xl, mu is the link to be flowed
// yl, nu is the link to calculate derivative
//
// given yl, m is the index of the pair: (xl, mu) (yl, nu)
//
// m: mu(1 xl), (mu+1)%4(4 xl), (mu+2)%4(4 xl), (mu+3)%4(4 xl)
// m: mu(6 xl), (mu+1)%4(4 xl), (mu+2)%4(4 xl), (mu+3)%4(4 xl)
{
  qassert(flow_type == 2);
  const Coordinate yg = geo.coordinate_g_from_l(yl);
  const int mask_yl = mask_from_coordinate(yg, flow_type);
  Coordinate xg;
  if (mask == mask_yl) {
    const int nu_size = 1;
    if (0 <= m and m < nu_size) {
      nu = mu;
      xg = yg;
    } else {
      set_xg_nu_from_mask_mu_yg_m_srect_nu_neq_mu(xg, nu, mask_yl, nu_size,
                                                  mask, mu, yg, m, flow_type);
    }
  } else if (mask == 3 - mask_yl) {
    const int nu_size = 6;
    if (0 <= m and m < nu_size) {
      nu = mu;
      const int k = m / 2;
      const int l = m % 2;
      int dir = mod(nu + k + 1, 4);
      if (0 == l) {
        xg = coordinate_shifts(yg, dir, dir);
      } else if (1 == l) {
        xg = coordinate_shifts(yg, -dir - 1, -dir - 1);
      } else {
        qassert(false);
      }
    } else {
      set_xg_nu_from_mask_mu_yg_m_srect_nu_neq_mu(xg, nu, mask_yl, nu_size,
                                                  mask, mu, yg, m, flow_type);
    }
  } else {
    qassert(false);
  }
  qassert(mask == mask_from_coordinate(xg, flow_type));
  xl = geo.coordinate_l_from_g(xg);
}

qacc void set_xl_nu_from_mask_mu_yl_m(Coordinate& xl, int& nu, const int mask,
                                      const int mu, const Coordinate& yl,
                                      const int m, const Geometry& geo,
                                      const int flow_type)
{
  if (1 == flow_type) {
    set_xl_nu_from_mask_mu_yl_m_plaq(xl, nu, mask, mu, yl, m, geo, flow_type);
  } else if (2 == flow_type) {
    set_xl_nu_from_mask_mu_yl_m_srect(xl, nu, mask, mu, yl, m, geo, flow_type);
  } else {
    qassert(false);
  }
}

qacc ColorMatrix gf_srect_staple_no_comm(const GaugeField& gf,
                                         const Coordinate& xl, const int mu)
// transpose the same way as gf.get_elem(xl, mu)
//
// only compute subset of rectangular staple (link on the shorter edge)
{
  ColorMatrix acc;
  set_zero(acc);
  for (int nu = -4; nu < 4; ++nu) {
    if (nu == mu or -nu - 1 == mu) {
      continue;
    }
    acc += gf_wilson_line_no_comm(
        gf, xl, make_array<int>(nu, nu, mu, -nu - 1, -nu - 1));
  }
  return acc;
}

qacc ColorMatrix gf_flow_staple_no_comm(const GaugeField& gf,
                                        const Coordinate& xl, const int mu,
                                        const int flow_size)
// transpose the same way as gf.get_elem(xl, mu)
//
// only compute subset of rectangular staple (link on the shorter edge)
{
  if (1 == flow_size) {
    return gf_plaq_staple_no_comm(gf, xl, mu);
  } else if (2 == flow_size) {
    return gf_srect_staple_no_comm(gf, xl, mu);
  } else {
    qassert(false);
    return ColorMatrix();
  }
}

inline void set_flow_staple_mask_mu_no_comm(FieldM<ColorMatrix, 1>& cf,
                                            const GaugeField& gf_ext,
                                            const int mask, const int mu,
                                            const int flow_size)
// cf will be initialized
// gf_ext need proper communication
// mask: flow 1:odd / 2:even site
// mu: flow link direction
{
  TIMER("set_flow_staple_mask_mu_no_comm");
  qassert(is_initialized(gf_ext));
  qassert(mask == 1 or mask == 2);
  qassert(0 <= mu and mu < 4);
  qassert(flow_size == 1 or flow_size == 2);
  const Geometry geo = geo_reform(gf_ext.geo());
  const vector<long>& flowed_indices =
      get_flowed_hmc_indices_mask_flow_size(geo, mask, flow_size);
  cf.init(geo);
  qacc_for(idx, flowed_indices.size(), {
    const long index = flowed_indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    cf.get_elem(xl) = gf_flow_staple_no_comm(gf_ext, xl, mu, flow_size);
  });
}

inline void gf_flow_plaq_mask_mu_no_comm(GaugeField& gf,
                                         const GaugeField& gf0_ext,
                                         const int mask, const int mu,
                                         const double epsilon,
                                         const int flow_size)
// mask: flow 1:odd / 2:even site
// mu: flow link direction
// epsilon: is the flow step size
{
  TIMER("gf_flow_plaq_mask_mu_no_comm");
  qassert(is_initialized(gf0_ext));
  qassert(mask == 1 or mask == 2);
  qassert(0 <= mu and mu < 4);
  const Geometry geo = geo_reform(gf0_ext.geo());
  const vector<long>& flowed_indices =
      get_flowed_hmc_indices_mask_flow_size(geo, mask, flow_size);
  gf.init(geo);
  gf = gf0_ext;
  FieldM<ColorMatrix, 1> cf;
  set_flow_staple_mask_mu_no_comm(cf, gf0_ext, mask, mu, flow_size);
  qacc_for(idx, flowed_indices.size(), {
    const long index = flowed_indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const ColorMatrix& c_x_mu = cf.get_elem(xl);
    const ColorMatrix& u0_x_mu = gf0_ext.get_elem(xl, mu);
    const ColorMatrix z_u_x_mu =
        -make_tr_less_anti_herm_matrix(u0_x_mu * matrix_adjoint(c_x_mu));
    const ColorMatrix e_z_u_x_mu = (Complex)epsilon * z_u_x_mu;
    const ColorMatrix e_u_x_mu = make_matrix_exp(e_z_u_x_mu) * u0_x_mu;
    ColorMatrix& u_x_mu = gf.get_elem(xl, mu);
    u_x_mu = e_u_x_mu;
  });
}

inline void gf_flow_inv_plaq_mask_mu_no_comm(
    GaugeField& gf, const GaugeField& gf1_ext, const int mask, const int mu,
    const double epsilon, const int flow_size, const int n_iter = 50)
// mask: flow 1:odd / 2:even site
// mu: flow link direction
// epsilon: is the flow step size
{
  TIMER("gf_flow_inv_plaq_mask_mu_no_comm");
  qassert(is_initialized(gf1_ext));
  qassert(mask == 1 or mask == 2);
  qassert(0 <= mu and mu < 4);
  const Geometry geo = geo_reform(gf1_ext.geo());
  const vector<long>& flowed_indices =
      get_flowed_hmc_indices_mask_flow_size(geo, mask, flow_size);
  gf.init(geo);
  gf = gf1_ext;
  FieldM<ColorMatrix, 1> cf;
  set_flow_staple_mask_mu_no_comm(cf, gf1_ext, mask, mu, flow_size);
  qacc_for(idx, flowed_indices.size(), {
    const long index = flowed_indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const ColorMatrix c_x_mu_dagger = matrix_adjoint(cf.get_elem(xl));
    const ColorMatrix& u1_x_mu = gf1_ext.get_elem(xl, mu);
    ColorMatrix u0_x_mu = u1_x_mu;
    for (int n = 0; n < n_iter; ++n) {
      const ColorMatrix x_u_x_mu =
          -make_tr_less_anti_herm_matrix(u0_x_mu * c_x_mu_dagger);
      const ColorMatrix e_x_u_x_mu = (Complex)epsilon * x_u_x_mu;
      u0_x_mu = make_matrix_exp(-e_x_u_x_mu) * u1_x_mu;
    }
    ColorMatrix& u_x_mu = gf.get_elem(xl, mu);
    u_x_mu = u0_x_mu;
  });
}

inline void set_marks_flow_plaq_mask_mu(CommMarks& marks, const Geometry& geo,
                                        const std::string& tag)
{
  TIMER_VERBOSE("set_marks_flow_plaq_mask_mu");
  qassert(geo.multiplicity == 4);
  marks.init();
  marks.init(geo);
  set_zero(marks);
  const std::vector<std::string> words = split_line_with_spaces(tag);
  qassert(words.size() == 3);
  const int mask = read_long(words[0]);
  const int mu = read_long(words[1]);
  const int flow_size = read_long(words[2]);
  qassert(mask == 1 or mask == 2);
  qassert(0 <= mu and mu < 4);
  const vector<long>& flowed_indices =
      get_flowed_hmc_indices_mask_flow_size(geo, mask, flow_size);
  qacc_for(idx, flowed_indices.size(), {
    const long index = flowed_indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    for (int nu = -4; nu < 4; ++nu) {
      if (nu == mu or -nu - 1 == mu) {
        continue;
      }
      if (1 == flow_size) {
        set_marks_field_path(marks, xl, make_array<int>(nu, mu, -nu - 1));
      } else if (2 == flow_size) {
        set_marks_field_path(marks, xl,
                             make_array<int>(nu, nu, mu, -nu - 1, -nu - 1));
      } else {
        qassert(false);
      }
    }
  });
}

inline void refresh_expanded_gf_flow_plaq_mask_mu(GaugeField& gf_ext,
                                                  const int mask, const int mu,
                                                  const int flow_size)
{
  TIMER("refresh_expanded_gf_flow_plaq_mask_mu");
  const CommPlan& plan =
      get_comm_plan(set_marks_flow_plaq_mask_mu,
                    ssprintf("%d %d %d", mask, mu, flow_size), gf_ext.geo());
  refresh_expanded(gf_ext, plan);
}

inline void gf_flow(GaugeField& gf, const GaugeField& gf0, const FlowInfo& fi)
//
// gf0 is the gauge field which we perform HMC evolution
// gf is supposed to be the desired configuration
//
// The flow operation typically smooth the gauge field.
//
// Normally call this at the end of Flowed-HMC evolution step.
{
  TIMER("gf_flow");
  const Coordinate expand_left(2, 2, 2, 2);
  const Coordinate expand_right(2, 2, 2, 2);
  const Geometry geo = geo_reform(gf0.geo());
  const Geometry geo_ext = geo_reform(gf0.geo(), 1, expand_left, expand_right);
  GaugeField gf_ext;
  gf_ext.init(geo_ext);
  gf_ext = gf0;
  for (int i = 0; i < (int)fi.v.size(); ++i) {
    const FlowStepInfo& fsi = fi.v[i];
    refresh_expanded_gf_flow_plaq_mask_mu(gf_ext, fsi.mask, fsi.mu,
                                          fsi.flow_size);
    gf_flow_plaq_mask_mu_no_comm(gf_ext, gf_ext, fsi.mask, fsi.mu, fsi.epsilon,
                                 fsi.flow_size);
  }
  gf.init(geo);
  gf = gf_ext;
}

inline void gf_flow_inv(GaugeField& gf, const GaugeField& gf1,
                        const FlowInfo& fi)
//
// gf1 is supposed to be the desired configuration
// gf is the gauge field which we perform HMC evolution
//
// The flow operation typically smooth the gauge field
//
// Normally call this at the beginning of Flowed-HMC evolution step.
{
  TIMER("gf_flow_inv");
  const Coordinate expand_left(2, 2, 2, 2);
  const Coordinate expand_right(2, 2, 2, 2);
  const Geometry geo = geo_reform(gf1.geo());
  const Geometry geo_ext = geo_reform(gf1.geo(), 1, expand_left, expand_right);
  GaugeField gf_ext;
  gf_ext.init(geo_ext);
  gf_ext = gf1;
  for (int i = (int)fi.v.size() - 1; i >= 0; --i) {
    const FlowStepInfo& fsi = fi.v[i];
    refresh_expanded_gf_flow_plaq_mask_mu(gf_ext, fsi.mask, fsi.mu,
                                          fsi.flow_size);
    gf_flow_inv_plaq_mask_mu_no_comm(gf_ext, gf_ext, fsi.mask, fsi.mu,
                                     fsi.epsilon, fsi.flow_size);
  }
  gf.init(geo);
  gf = gf_ext;
}

qacc array<ColorMatrix, 2> d_uc_mat_plaq_site_no_comm(
    const FieldM<ColorMatrix, 1>& cf, const GaugeField& gf_ext,
    const Coordinate& xl, const int mu, const Coordinate& yl, const int nu)
{
  array<ColorMatrix, 2> uc_mats;
  ColorMatrix& uc_pre = uc_mats[0];
  ColorMatrix& uc_post = uc_mats[1];
  const ColorMatrix& u = gf_ext.get_elem(xl, mu);
  if (mu == nu and xl == yl) {
    const ColorMatrix c_dagger = matrix_adjoint(cf.get_elem(xl));
    set_unit(uc_pre);
    uc_post = u * c_dagger;
    return uc_mats;
  }
  const Coordinate& xl0 = xl;
  int dir1, dir2, dir3;
  Coordinate xl1, xl2;
  int n_step = 0;
  int is_pos_dir = 0;
  for (int k = -4; k < 4; ++k) {
    dir1 = k;
    dir2 = mu;
    dir3 = -k - 1;
    xl1 = coordinate_shifts(xl0, dir1);
    xl2 = coordinate_shifts(xl1, dir2);
    is_pos_dir = is_same_link(xl0, dir1, yl, nu);
    if (0 != is_pos_dir) {
      n_step = 1;
      break;
    }
    is_pos_dir = is_same_link(xl1, dir2, yl, nu);
    if (0 != is_pos_dir) {
      n_step = 2;
      break;
    }
    is_pos_dir = is_same_link(xl2, dir3, yl, nu);
    if (0 != is_pos_dir) {
      n_step = 3;
      break;
    }
  }
  ColorMatrix d_c_d_s_pre, d_c_d_s_u, d_c_d_s_post;
  if (1 == n_step) {
    set_unit(d_c_d_s_pre);
    d_c_d_s_u = gf_get_link(gf_ext, xl0, dir1);
    d_c_d_s_post =
        gf_get_link(gf_ext, xl1, dir2) * gf_get_link(gf_ext, xl2, dir3);
  } else if (2 == n_step) {
    d_c_d_s_pre = gf_get_link(gf_ext, xl0, dir1);
    d_c_d_s_u = gf_get_link(gf_ext, xl1, dir2);
    d_c_d_s_post = gf_get_link(gf_ext, xl2, dir3);
  } else if (3 == n_step) {
    d_c_d_s_pre =
        gf_get_link(gf_ext, xl0, dir1) * gf_get_link(gf_ext, xl1, dir2);
    d_c_d_s_u = gf_get_link(gf_ext, xl2, dir3);
    set_unit(d_c_d_s_post);
  } else {
    qassert(false);
  }
  if (1 == is_pos_dir) {
    d_c_d_s_post = d_c_d_s_u * d_c_d_s_post;
  } else if (-1 == is_pos_dir) {
    d_c_d_s_pre =
        -d_c_d_s_pre * d_c_d_s_u;  // - sign due to dagger of d_c_d_s_u for T^b
  } else {
    qassert(false);
  }
  uc_pre =
      -u * matrix_adjoint(d_c_d_s_post);  // - sign due to c^\dagger for T^b
  uc_post = matrix_adjoint(d_c_d_s_pre);
  return uc_mats;
}

qacc array<ColorMatrix, 2> d_uc_mat_srect_site_no_comm(
    const FieldM<ColorMatrix, 1>& cf, const GaugeField& gf_ext,
    const Coordinate& xl, const int mu, const Coordinate& yl, const int nu)
{
  array<ColorMatrix, 2> uc_mats;
  ColorMatrix& uc_pre = uc_mats[0];
  ColorMatrix& uc_post = uc_mats[1];
  const ColorMatrix& u = gf_ext.get_elem(xl, mu);
  if (mu == nu and xl == yl) {
    const ColorMatrix c_dagger = matrix_adjoint(cf.get_elem(xl));
    set_unit(uc_pre);
    uc_post = u * c_dagger;
    return uc_mats;
  }
  const Coordinate& xl0 = xl;
  int dir1, dir2, dir3;
  Coordinate xl1, xl2, xl3, xl4;
  int n_step = 0;
  int is_pos_dir = 0;
  for (int k = -4; k < 4; ++k) {
    dir1 = k;
    dir2 = mu;
    dir3 = -k - 1;
    xl1 = coordinate_shifts(xl0, dir1);
    xl2 = coordinate_shifts(xl1, dir1);
    xl3 = coordinate_shifts(xl2, dir2);
    xl4 = coordinate_shifts(xl3, dir3);
    is_pos_dir = is_same_link(xl0, dir1, yl, nu);
    if (0 != is_pos_dir) {
      n_step = 1;
      break;
    }
    is_pos_dir = is_same_link(xl1, dir1, yl, nu);
    if (0 != is_pos_dir) {
      n_step = 2;
      break;
    }
    is_pos_dir = is_same_link(xl2, dir2, yl, nu);
    if (0 != is_pos_dir) {
      n_step = 3;
      break;
    }
    is_pos_dir = is_same_link(xl3, dir3, yl, nu);
    if (0 != is_pos_dir) {
      n_step = 4;
      break;
    }
    is_pos_dir = is_same_link(xl4, dir3, yl, nu);
    if (0 != is_pos_dir) {
      n_step = 5;
      break;
    }
  }
  if (n_step == 0) {
#ifndef QLAT_USE_ACC
    displayln(ssprintf("xl=%s ; mu=%d ; yl=%s ; nu=%d", show(xl).c_str(), mu,
                       show(yl).c_str(), nu));
#endif
    qassert(false);
  }
  ColorMatrix d_c_d_s_pre, d_c_d_s_u, d_c_d_s_post;
  if (1 == n_step) {
    set_unit(d_c_d_s_pre);
    d_c_d_s_u = gf_get_link(gf_ext, xl0, dir1);
    d_c_d_s_post =
        gf_get_link(gf_ext, xl1, dir1) * gf_get_link(gf_ext, xl2, dir2) *
        gf_get_link(gf_ext, xl3, dir3) * gf_get_link(gf_ext, xl4, dir3);
  } else if (2 == n_step) {
    d_c_d_s_pre = gf_get_link(gf_ext, xl0, dir1);
    d_c_d_s_u = gf_get_link(gf_ext, xl1, dir1);
    d_c_d_s_post = gf_get_link(gf_ext, xl2, dir2) *
                   gf_get_link(gf_ext, xl3, dir3) *
                   gf_get_link(gf_ext, xl4, dir3);
  } else if (3 == n_step) {
    d_c_d_s_pre =
        gf_get_link(gf_ext, xl0, dir1) * gf_get_link(gf_ext, xl1, dir1);
    d_c_d_s_u = gf_get_link(gf_ext, xl2, dir2);
    d_c_d_s_post =
        gf_get_link(gf_ext, xl3, dir3) * gf_get_link(gf_ext, xl4, dir3);
  } else if (4 == n_step) {
    d_c_d_s_pre = gf_get_link(gf_ext, xl0, dir1) *
                  gf_get_link(gf_ext, xl1, dir1) *
                  gf_get_link(gf_ext, xl2, dir2);
    d_c_d_s_u = gf_get_link(gf_ext, xl3, dir3);
    d_c_d_s_post = gf_get_link(gf_ext, xl4, dir3);
  } else if (5 == n_step) {
    d_c_d_s_pre =
        gf_get_link(gf_ext, xl0, dir1) * gf_get_link(gf_ext, xl1, dir1) *
        gf_get_link(gf_ext, xl2, dir2) * gf_get_link(gf_ext, xl3, dir3);
    d_c_d_s_u = gf_get_link(gf_ext, xl4, dir3);
    set_unit(d_c_d_s_post);
  } else {
    qassert(false);
  }
  if (1 == is_pos_dir) {
    d_c_d_s_post = d_c_d_s_u * d_c_d_s_post;
  } else if (-1 == is_pos_dir) {
    d_c_d_s_pre =
        -d_c_d_s_pre * d_c_d_s_u;  // - sign due to dagger of d_c_d_s_u for T^b
  } else {
    qassert(false);
  }
  uc_pre =
      -u * matrix_adjoint(d_c_d_s_post);  // - sign due to c^\dagger for T^b
  uc_post = matrix_adjoint(d_c_d_s_pre);
  return uc_mats;
}

qacc array<ColorMatrix, 2> d_uc_mat_site_no_comm(
    const FieldM<ColorMatrix, 1>& cf, const GaugeField& gf_ext,
    const Coordinate& xl, const int mu, const Coordinate& yl, const int nu,
    const int flow_size)
{
  if (1 == flow_size) {
    return d_uc_mat_plaq_site_no_comm(cf, gf_ext, xl, mu, yl, nu);
  } else if (2 == flow_size) {
    return d_uc_mat_srect_site_no_comm(cf, gf_ext, xl, mu, yl, nu);
  } else {
    qassert(false);
    return array<ColorMatrix, 2>();
  }
}

inline void set_d_uc_mat_plaq_mask_mu_no_comm(
    Field<array<ColorMatrix, 2> >& ducf, const FieldM<ColorMatrix, 1>& cf,
    const GaugeField& gf_ext, const int mask, const int mu, const int flow_size)
// ducf does NOT need to be initialized.
// It will be initialized with no expansion
//
// See set_xl_nu_from_mask_mu_yl_m, n_mat_plaq_site_no_comm
//
// gf_ext need proper communication
// mask: flow 1:odd / 2:even site
// mu: flow link direction
{
  TIMER("set_d_uc_mat_plaq_mask_mu_no_comm");
  qassert(is_initialized(cf));
  qassert(cf.geo().eo == 0);
  qassert(is_initialized(gf_ext));
  qassert(mask == 1 or mask == 2);
  qassert(0 <= mu and mu < 4);
  const Geometry geo = geo_reform(gf_ext.geo());
  ducf.init(geo, multiplicity_flow_hmc_max(flow_size));
  qacc_for(index, geo.local_volume(), {
    const Coordinate yl = geo.coordinate_from_index(index);
    const Coordinate yg = geo.coordinate_g_from_l(yl);
    const int mask_yl = mask_from_coordinate(yg, flow_size);
    Vector<array<ColorMatrix, 2> > ducfv = ducf.get_elems(yl);
    for (int m = 0; m < multiplicity_flow_hmc(mask_yl == mask, flow_size);
         ++m) {
      Coordinate xl;
      int nu;
      set_xl_nu_from_mask_mu_yl_m(xl, nu, mask, mu, yl, m, geo, flow_size);
      ducfv[m] = d_uc_mat_site_no_comm(cf, gf_ext, xl, mu, yl, nu, flow_size);
    }
  });
}

qacc AdjointColorMatrix n_mat_plaq_site_no_comm(
    const array<ColorMatrix, 2>& uc_mats, const ColorMatrixConstants& cmcs)
{
  const array<ColorMatrix, 8>& ts = cmcs.ts;
  const ColorMatrix& uc_pre = uc_mats[0];
  const ColorMatrix& uc_post = uc_mats[1];
  AdjointColorMatrix n_mat;
  for (int b = 0; b < 8; ++b) {
    const ColorMatrix d_c_d_s = uc_pre * ts[b] * uc_post;
    const ColorMatrix n_b = make_tr_less_anti_herm_matrix(d_c_d_s);
    const array<double, 8> basis_b =
        basis_projection_anti_hermitian_matrix(n_b);
    for (int a = 0; a < 8; ++a) {
      n_mat(a, b) = basis_b[a];
    }
  }
  return n_mat;
}

inline void set_n_mat_plaq_mask_mu_no_comm(
    Field<AdjointColorMatrix>& nf, const Field<array<ColorMatrix, 2> >& ducf,
    const int mask, const int flow_size)
// nf does NOT need to be initialized.
// It will be initialized with ducf geometry.
{
  TIMER("set_n_mat_plaq_mask_mu_no_comm");
  qassert(is_initialized(ducf));
  qassert(ducf.geo().multiplicity == multiplicity_flow_hmc_max(flow_size));
  const Geometry& geo = ducf.geo();
  nf.init(geo, multiplicity_flow_hmc_max(flow_size));
  const box<ColorMatrixConstants>& cmcs =
      ColorMatrixConstants::get_instance_box();
  qacc_for(index, geo.local_volume(), {
    const Coordinate yl = geo.coordinate_from_index(index);
    const Coordinate yg = geo.coordinate_g_from_l(yl);
    const int mask_yl = mask_from_coordinate(yg, flow_size);
    const Vector<array<ColorMatrix, 2> > ducfv = ducf.get_elems_const(yl);
    Vector<AdjointColorMatrix> nfv = nf.get_elems(yl);
    for (int m = 0; m < multiplicity_flow_hmc(mask_yl == mask, flow_size);
         ++m) {
      nfv[m] = n_mat_plaq_site_no_comm(ducfv[m], cmcs());
    }
  });
}

qacc AdjointColorMatrix n_mat_plaq_site_no_comm(
    const FieldM<ColorMatrix, 1>& cf, const GaugeField& gf,
    const Coordinate& xl, const int mu, const ColorMatrixConstants& cmcs)
{
  const array<ColorMatrix, 8>& ts = cmcs.ts;
  const ColorMatrix& u = gf.get_elem(xl, mu);
  const ColorMatrix c_dagger = matrix_adjoint(cf.get_elem(xl));
  AdjointColorMatrix n_mat;
  for (int b = 0; b < 8; ++b) {
    const ColorMatrix n_b = make_tr_less_anti_herm_matrix(ts[b] * u * c_dagger);
    const array<double, 8> basis_b =
        basis_projection_anti_hermitian_matrix(n_b);
    for (int a = 0; a < 8; ++a) {
      n_mat(a, b) = basis_b[a];
    }
  }
  return n_mat;
}

inline void set_n_mat_plaq_mask_mu_no_comm(FieldM<AdjointColorMatrix, 1>& nf,
                                           const FieldM<ColorMatrix, 1>& cf,
                                           const GaugeField& gf_ext,
                                           const int mask, const int mu,
                                           const int flow_size)
// nf does NOT need to be initialized.
// It will be initialized with no expansion
//
// nf: only xl, mu
//
// gf_ext need proper communication
// mask: flow 1:odd / 2:even site
// mu: flow link direction
{
  TIMER("set_n_mat_plaq_mask_mu_no_comm(nf)");
  qassert(is_initialized(cf));
  qassert(is_initialized(gf_ext));
  qassert(mask == 1 or mask == 2);
  qassert(0 <= mu and mu < 4);
  const box<ColorMatrixConstants>& cmcs =
      ColorMatrixConstants::get_instance_box();
  const Geometry geo = geo_reform(gf_ext.geo());
  const vector<long>& flowed_indices =
      get_flowed_hmc_indices_mask_flow_size(geo, mask, flow_size);
  nf.init(geo);
  qacc_for(idx, flowed_indices.size(), {
    const long index = flowed_indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    AdjointColorMatrix& n_mat = nf.get_elem(xl);
    n_mat = n_mat_plaq_site_no_comm(cf, gf_ext, xl, mu, cmcs());
  });
}

inline void set_ad_x_and_j_n_x_plaq_mask_mu_no_comm(
    FieldM<AdjointColorMatrix, 2>& f_ad_x_and_j_n_x,
    const FieldM<ColorMatrix, 1>& cf, const GaugeField& gf, const int mask,
    const int mu, const double epsilon, const int flow_size)
{
  TIMER("set_ad_x_and_j_n_x_plaq_mask_mu_no_comm");
  qassert(is_initialized(cf));
  qassert(is_initialized(gf));
  qassert(mask == 1 or mask == 2);
  qassert(0 <= mu and mu < 4);
  const box<ColorMatrixConstants>& cmcs =
      ColorMatrixConstants::get_instance_box();
  const Geometry geo = geo_reform(gf.geo());
  const vector<long>& flowed_indices =
      get_flowed_hmc_indices_mask_flow_size(geo, mask, flow_size);
  f_ad_x_and_j_n_x.init(geo);
  qacc_for(idx, flowed_indices.size(), {
    const long index = flowed_indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const ColorMatrix& u = gf.get_elem(xl, mu);
    const ColorMatrix c_dagger = matrix_adjoint(cf.get_elem(xl));
    const ColorMatrix x_mat =
        (Complex)(-epsilon) * make_tr_less_anti_herm_matrix(u * c_dagger);
    Vector<AdjointColorMatrix> ad_x_and_j_n_x = f_ad_x_and_j_n_x.get_elems(xl);
    AdjointColorMatrix& ad_x_mat = ad_x_and_j_n_x[0];
    AdjointColorMatrix& j_n_x_mat = ad_x_and_j_n_x[1];
    ad_x_mat = make_adjoint_representation(x_mat, cmcs());
    j_n_x_mat = make_diff_exp_map(-ad_x_mat);
  });
}

qacc AdjointColorMatrix m_mat_plaq_site_no_comm(
    const AdjointColorMatrix& n_mat,
    const Vector<AdjointColorMatrix> ad_x_and_j_n_x, const Coordinate& xl,
    const int mu, const Coordinate& yl, const int nu, const double epsilon)
{
  const AdjointColorMatrix& ad_x_mat = ad_x_and_j_n_x[0];
  const AdjointColorMatrix& j_n_x_mat = ad_x_and_j_n_x[1];
  if (mu == nu and xl == yl) {
    return make_matrix_exp(ad_x_mat) - epsilon * j_n_x_mat * n_mat;
  } else {
    return -epsilon * j_n_x_mat * n_mat;
  }
}

inline void set_m_mat_plaq_mask_mu_no_comm(
    Field<AdjointColorMatrix>& mf, const Field<AdjointColorMatrix>& nf,
    const FieldM<AdjointColorMatrix, 2>& f_ad_x_and_j_n_x_ext, const int mask,
    const int mu, const double epsilon, const int flow_size)
// mf and nf have the same structure
//
// See set_xl_nu_from_mask_mu_yl_m
//
// f_ad_x_and_j_n_x_ext need proper communication
// mask: flow 1:odd / 2:even site
// mu: flow link direction
{
  TIMER("set_m_mat_plaq_mask_mu_no_comm");
  qassert(is_initialized(nf));
  qassert(nf.geo().multiplicity == multiplicity_flow_hmc_max(flow_size));
  qassert(is_initialized(f_ad_x_and_j_n_x_ext));
  qassert(mask == 1 or mask == 2);
  qassert(0 <= mu and mu < 4);
  const Geometry& geo = nf.geo();
  mf.init(geo, multiplicity_flow_hmc_max(flow_size));
  qacc_for(index, geo.local_volume(), {
    const Coordinate yl = geo.coordinate_from_index(index);
    const Coordinate yg = geo.coordinate_g_from_l(yl);
    const int mask_yl = mask_from_coordinate(yg, flow_size);
    Vector<AdjointColorMatrix> mfv = mf.get_elems(yl);
    const Vector<AdjointColorMatrix> nfv = nf.get_elems_const(yl);
    qassert(mfv.size() == nfv.size());
    for (int m = 0; m < multiplicity_flow_hmc(mask == mask_yl, flow_size);
         ++m) {
      Coordinate xl;
      int nu;
      set_xl_nu_from_mask_mu_yl_m(xl, nu, mask, mu, yl, m, geo, flow_size);
      const Vector<AdjointColorMatrix> ad_x_and_j_n_x =
          f_ad_x_and_j_n_x_ext.get_elems_const(xl);
      mfv[m] = m_mat_plaq_site_no_comm(nfv[m], ad_x_and_j_n_x, xl, mu, yl, nu,
                                       epsilon);
    }
  });
}

qacc AdjointColorMatrix mp_mat_plaq_site_no_comm(
    const AdjointColorMatrix& n_mat, const FieldM<ColorMatrix, 1>& cf,
    const GaugeField& gf, const Coordinate& xl, const int mu,
    const double epsilon, const ColorMatrixConstants& cmcs)
{
  const ColorMatrix& u = gf.get_elem(xl, mu);
  const ColorMatrix c_dagger = matrix_adjoint(cf.get_elem(xl));
  const ColorMatrix x_mat =
      (Complex)(-epsilon) * make_tr_less_anti_herm_matrix(u * c_dagger);
  const AdjointColorMatrix j_x_mat = make_diff_exp_map(x_mat, cmcs);
  AdjointColorMatrix m_mat;
  set_unit(m_mat);
  m_mat -= epsilon * j_x_mat * n_mat;
  return m_mat;
}

inline void set_mp_mat_plaq_mask_mu_no_comm(
    FieldM<AdjointColorMatrix, 1>& mpf, const Field<AdjointColorMatrix>& nf,
    const FieldM<ColorMatrix, 1>& cf, const GaugeField& gf, const int mask,
    const int mu, const double epsilon, const int flow_size)
// only use the first elem in each site of nf
//
// mask: flow 1:odd / 2:even site
// mu: flow link direction
{
  TIMER("set_mp_mat_plaq_no_comm(mpf)");
  qassert(is_initialized(nf));
  qassert(is_initialized(cf));
  qassert(is_initialized(gf));
  qassert(mask == 1 or mask == 2);
  qassert(0 <= mu and mu < 4);
  const box<ColorMatrixConstants>& cmcs =
      ColorMatrixConstants::get_instance_box();
  const Geometry geo = geo_reform(gf.geo());
  const vector<long>& flowed_indices =
      get_flowed_hmc_indices_mask_flow_size(geo, mask, flow_size);
  mpf.init(geo);
  qacc_for(idx, flowed_indices.size(), {
    const long index = flowed_indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const AdjointColorMatrix& n_mat = nf.get_elem(xl, 0);
    AdjointColorMatrix& m_mat = mpf.get_elem(xl);
    m_mat = mp_mat_plaq_site_no_comm(n_mat, cf, gf, xl, mu, epsilon, cmcs());
  });
}

qacc void set_gm_force_from_flow_site_no_comm(
    GaugeMomentum& gm_force, const GaugeMomentum& gm_force_pre_ext,
    const Vector<AdjointColorMatrix> m_mat_vec, const int mask, const int mu,
    const int flow_size, const Coordinate& yl, const Geometry& geo)
{
  const Coordinate yg = geo.coordinate_g_from_l(yl);
  const int mask_yl = mask_from_coordinate(yg, flow_size);
  Vector<ColorMatrix> fv = gm_force.get_elems(yl);
  set_zero(fv);
  for (int m = 0; m < multiplicity_flow_hmc(mask == mask_yl, flow_size); ++m) {
    Coordinate xl;
    int nu;
    set_xl_nu_from_mask_mu_yl_m(xl, nu, mask, mu, yl, m, geo, flow_size);
    const AdjointColorMatrix& m_mat = m_mat_vec[m];
    const ColorMatrix& force_x_mu = gm_force_pre_ext.get_elem(xl, mu);
    const array<double, 8> basis_b =
        basis_projection_anti_hermitian_matrix(force_x_mu);
    array<double, 8> basis_a;
    set_zero(basis_a);
    for (int b = 0; b < 8; ++b) {
      for (int a = 0; a < 8; ++a) {
        basis_a[a] += basis_b[b] * m_mat(b, a);
      }
    }
    fv[nu] += make_anti_hermitian_matrix(basis_a);
  }
  const Vector<ColorMatrix> f_pre_v = gm_force_pre_ext.get_elems_const(yl);
  for (int nu = 0; nu < 4; ++nu) {
    if (not(mask_yl == mask and mu == nu)) {
      fv[nu] += f_pre_v[nu];
    }
  }
}

inline void set_gm_force_from_flow_no_comm(
    GaugeMomentum& gm_force, const GaugeMomentum& gm_force_pre_ext,
    const Field<AdjointColorMatrix>& mf, const int mask, const int mu,
    const int flow_size)
// mask: flow 1:odd / 2:even site
// mu: flow link direction
{
  TIMER("set_gm_force_from_flow_no_comm");
  qassert(is_initialized(gm_force_pre_ext));
  qassert(is_initialized(mf));
  qassert(mf.geo().multiplicity == multiplicity_flow_hmc_max(flow_size));
  qassert(mask == 1 or mask == 2);
  qassert(0 <= mu and mu < 4);
  const Geometry geo = geo_reform(gm_force_pre_ext.geo());
  gm_force.init(geo);
  qacc_for(index, geo.local_volume(), {
    const Coordinate yl = geo.coordinate_from_index(index);
    const Vector<AdjointColorMatrix> m_mat_vec = mf.get_elems_const(yl);
    set_gm_force_from_flow_site_no_comm(gm_force, gm_force_pre_ext, m_mat_vec,
                                        mask, mu, flow_size, yl, geo);
  });
}

inline void set_f_det_util_plaq_mask_mu(
    FieldM<array<double, 8>, 1>& f_e2_dj_x_n_mp_inv,
    FieldM<AdjointColorMatrix, 1>& f_n_e_mp_inv_j_x,
    const FieldM<AdjointColorMatrix, 1>& mpf,
    const Field<AdjointColorMatrix>& nf, const FieldM<ColorMatrix, 1>& cf,
    const GaugeField& gf, const int mask, const int mu, const double epsilon,
    const int flow_size)
// mask: flow 1:odd / 2:even site
// mu: flow link direction
{
  TIMER("set_f_det_util_plaq_mask_mu");
  qassert(is_initialized(mpf));
  qassert(is_initialized(nf));
  qassert(is_initialized(cf));
  qassert(is_initialized(gf));
  qassert(mask == 1 or mask == 2);
  qassert(0 <= mu and mu < 4);
  const box<ColorMatrixConstants>& cmcs =
      ColorMatrixConstants::get_instance_box();
  const Geometry geo = geo_reform(gf.geo());
  const vector<long>& flowed_indices =
      get_flowed_hmc_indices_mask_flow_size(geo, mask, flow_size);
  f_n_e_mp_inv_j_x.init(geo);
  f_e2_dj_x_n_mp_inv.init(geo);
  qacc_for(idx, flowed_indices.size(), {
    const long index = flowed_indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const AdjointColorMatrix& mp_mat = mpf.get_elem(xl);
    const AdjointColorMatrix mp_inv_mat = matrix_inverse(mp_mat);
    const AdjointColorMatrix& n_mat = nf.get_elem(xl, 0);
    const ColorMatrix& u = gf.get_elem(xl, mu);
    const ColorMatrix c_dagger = matrix_adjoint(cf.get_elem(xl));
    const ColorMatrix x_mat =
        (Complex)(-epsilon) * make_tr_less_anti_herm_matrix(u * c_dagger);
    const AdjointColorMatrix j_x_mat = make_diff_exp_map(x_mat, cmcs());
    const AdjointColorMatrix e2_n_mp_inv_mat =
        sqr(epsilon) * n_mat * mp_inv_mat;
    array<double, 8>& basis = f_e2_dj_x_n_mp_inv.get_elem(xl);
    for (int e = 0; e < 8; ++e) {
      basis[e] = matrix_trace(make_diff_exp_map_diff(x_mat, e, cmcs()),
                              e2_n_mp_inv_mat)
                     .real();
    }
    f_n_e_mp_inv_j_x.get_elem(xl) = (-epsilon) * mp_inv_mat * j_x_mat;
  });
}

inline void set_gm_force_from_flow_det_no_comm(
    GaugeMomentum& gm_force_det,
    const FieldM<array<double, 8>, 1>& f_e2_dj_x_n_mp_inv_ext,
    const FieldM<AdjointColorMatrix, 1>& f_n_e_mp_inv_j_x_ext,
    const Field<AdjointColorMatrix>& nf,
    const Field<array<ColorMatrix, 2> >& ducf, const int mask, const int mu,
    const int flow_size)
// See set_xl_nu_from_mask_mu_yl_m
// mask: flow 1:odd / 2:even site
// mu: flow link direction
{
  TIMER("set_gm_force_from_flow_det_no_comm");
  const box<ColorMatrixConstants>& cmcs =
      ColorMatrixConstants::get_instance_box();
  qassert(is_initialized(f_e2_dj_x_n_mp_inv_ext));
  qassert(is_initialized(f_n_e_mp_inv_j_x_ext));
  qassert(is_initialized(nf));
  qassert(nf.geo().multiplicity == multiplicity_flow_hmc_max(flow_size));
  qassert(is_initialized(ducf));
  qassert(ducf.geo().multiplicity == multiplicity_flow_hmc_max(flow_size));
  qassert(mask == 1 or mask == 2);
  qassert(0 <= mu and mu < 4);
  const Geometry geo = geo_reform(nf.geo());
  gm_force_det.init();
  gm_force_det.init(geo);
  qacc_for(index, geo.local_volume(), {
    const Coordinate yl = geo.coordinate_from_index(index);
    const Coordinate yg = geo.coordinate_g_from_l(yl);
    const int mask_yl = mask_from_coordinate(yg, flow_size);
    const array<ColorMatrix, 8>& ts = cmcs().ts;
    const Vector<AdjointColorMatrix> nfv = nf.get_elems_const(yl);
    const Vector<array<ColorMatrix, 2> > ducfv = ducf.get_elems_const(yl);
    Vector<ColorMatrix> gm_f_v = gm_force_det.get_elems(yl);
    set_zero(gm_f_v);
    for (int m = 0; m < multiplicity_flow_hmc(mask == mask_yl, flow_size);
         ++m) {
      Coordinate xl;
      int nu;
      set_xl_nu_from_mask_mu_yl_m(xl, nu, mask, mu, yl, m, geo, flow_size);
      const array<double, 8>& e2_dj_x_n_mp_inv =
          f_e2_dj_x_n_mp_inv_ext.get_elem(xl);
      const AdjointColorMatrix& n_e_mp_inv_j_x_mat =
          f_n_e_mp_inv_j_x_ext.get_elem(xl);
      const AdjointColorMatrix& n_mat = nfv[m];
      const array<ColorMatrix, 2>& uc_mats = ducfv[m];
      const ColorMatrix& uc_pre = uc_mats[0];
      const ColorMatrix& uc_post = uc_mats[1];
      array<double, 8> f_det_basis;
      set_zero(f_det_basis);
      for (int a = 0; a < 8; ++a) {
        for (int e = 0; e < 8; ++e) {
          f_det_basis[a] += n_mat(e, a) * e2_dj_x_n_mp_inv[e];
        }
        const ColorMatrix uc = uc_pre * ts[a] * uc_post;
        for (int c = 0; c < 8; ++c) {
          const ColorMatrix d_n = make_tr_less_anti_herm_matrix(ts[c] * uc);
          const array<double, 8> d_n_b =
              basis_projection_anti_hermitian_matrix(d_n);
          for (int b = 0; b < 8; ++b) {
            f_det_basis[a] += n_e_mp_inv_j_x_mat(c, b) * d_n_b[b];
          }
        }
      }
      const ColorMatrix f_det =
          (Complex)0.5 * make_anti_hermitian_matrix(f_det_basis);
      gm_f_v[nu] += f_det;
    }
  });
}

inline double mf_ln_det_sum(const Field<AdjointColorMatrix>& mpf,
                            const int mask, const int flow_size)
{
  TIMER("mf_ln_det_sum");
  const Geometry& geo = mpf.geo();
  const vector<long>& flowed_indices =
      get_flowed_hmc_indices_mask_flow_size(geo, mask, flow_size);
  FieldM<double, 1> f_ln_det;
  f_ln_det.init(geo);
  set_zero(f_ln_det);
  qacc_for(idx, flowed_indices.size(), {
    const long index = flowed_indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const AdjointColorMatrix& m_mat = mpf.get_elem(xl, 0);
    f_ln_det.get_elem(xl) = std::log(matrix_determinant(m_mat));
  });
  double sum = 0.0;
  for (long index = 0; index < geo.local_volume(); ++index) {
    sum += f_ln_det.get_elem(index);
  }
  return sum;
}

inline double gf_flow_and_ln_det_node(GaugeField& gf, const GaugeField& gf0,
                                      const FlowInfo& fi)
//
// Return ln(det(d gf / d gf0)) of the flow (on this node ONLY).
// And set gf to be the flowed gauge field from gf0.
//
// Normally call this function to compute the hamilton for this node.
//
// Note the glb_sum is NOT perform for the return value ln_det_node.
//
// The ln_det_node returned is to be subtracted from the orignal gauge action.
//
{
  TIMER("gf_flow_ln_det_node");
  const Coordinate expand_left(2, 2, 2, 2);
  const Coordinate expand_right(2, 2, 2, 2);
  const Geometry geo = geo_reform(gf0.geo());
  const Geometry geo_ext = geo_reform(gf0.geo(), 1, expand_left, expand_right);
  GaugeField gf_ext;
  gf_ext.init(geo_ext);
  gf_ext = gf0;
  double ln_det_node = 0.0;
  for (int i = 0; i < (int)fi.v.size(); ++i) {
    const FlowStepInfo& fsi = fi.v[i];
    refresh_expanded_gf_flow_plaq_mask_mu(gf_ext, fsi.mask, fsi.mu,
                                          fsi.flow_size);
    FieldM<ColorMatrix, 1> cf;
    set_flow_staple_mask_mu_no_comm(cf, gf_ext, fsi.mask, fsi.mu,
                                    fsi.flow_size);
    FieldM<AdjointColorMatrix, 1> nf;
    set_n_mat_plaq_mask_mu_no_comm(nf, cf, gf_ext, fsi.mask, fsi.mu,
                                   fsi.flow_size);
    FieldM<AdjointColorMatrix, 1> mpf;
    set_mp_mat_plaq_mask_mu_no_comm(mpf, nf, cf, gf_ext, fsi.mask, fsi.mu,
                                    fsi.epsilon, fsi.flow_size);
    ln_det_node += mf_ln_det_sum(mpf, fsi.mask, fsi.flow_size);
    gf_flow_plaq_mask_mu_no_comm(gf_ext, gf_ext, fsi.mask, fsi.mu, fsi.epsilon,
                                 fsi.flow_size);
  }
  gf.init(geo);
  gf = gf_ext;
  return ln_det_node;
}

inline double gf_hamilton_flowed_node(const GaugeField& gf0,
                                      const GaugeAction& ga, const FlowInfo& fi)
{
  TIMER("gf_hamilton_flowed_node");
  GaugeField gf;
  const double ln_det_node = gf_flow_and_ln_det_node(gf, gf0, fi);
  return gf_hamilton_node(gf, ga) - ln_det_node;
}

inline void gf_hamilton_flowed(double& flowed_action, double& ln_det,
                               const GaugeField& gf0, const GaugeAction& ga,
                               const FlowInfo& fi)
// effective hamilton for U_0 = flowed_action - ln_det
{
  TIMER("gf_hamilton_flowed");
  GaugeField gf;
  ln_det = gf_flow_and_ln_det_node(gf, gf0, fi);
  flowed_action = gf_hamilton_node(gf, ga);
  glb_sum(ln_det);
  glb_sum(flowed_action);
}

inline void set_gm_force_propagated_from_flow_step(
    GaugeMomentum& gm_force_propagated, const GaugeMomentum& gm_force_pre,
    const GaugeField& gf0_ext, const FlowStepInfo& fsi)
// The gf0_ext must be extended and refreshed and not flowed
// (for this flow step)
{
  TIMER("set_gm_force_propagated_from_flow_step");
  qassert(is_initialized(gm_force_pre));
  qassert(is_initialized(gf0_ext));
  const Coordinate expand_left(fsi.flow_size, fsi.flow_size, fsi.flow_size,
                               fsi.flow_size);
  const Coordinate expand_right(fsi.flow_size, fsi.flow_size, fsi.flow_size,
                                fsi.flow_size);
  const Geometry geo = geo_reform(gf0_ext.geo(), 1);
  const Geometry geo_ext = geo_reform(geo, 1, expand_left, expand_right);
  GaugeMomentum gm_force_pre_ext;
  gm_force_pre_ext.init(geo_ext);
  gm_force_pre_ext = gm_force_pre;
  refresh_expanded(gm_force_pre_ext);
  // clear the output forces
  gm_force_propagated.init();
  // To propagate force
  FieldM<ColorMatrix, 1> cf;
  set_flow_staple_mask_mu_no_comm(cf, gf0_ext, fsi.mask, fsi.mu, fsi.flow_size);
  Field<array<ColorMatrix, 2> > ducf;
  set_d_uc_mat_plaq_mask_mu_no_comm(ducf, cf, gf0_ext, fsi.mask, fsi.mu,
                                    fsi.flow_size);
  Field<AdjointColorMatrix> nf;
  set_n_mat_plaq_mask_mu_no_comm(nf, ducf, fsi.mask, fsi.flow_size);
  FieldM<AdjointColorMatrix, 2> f_ad_x_and_j_n_x_ext;
  f_ad_x_and_j_n_x_ext.init(geo_ext);
  set_ad_x_and_j_n_x_plaq_mask_mu_no_comm(f_ad_x_and_j_n_x_ext, cf, gf0_ext,
                                          fsi.mask, fsi.mu, fsi.epsilon,
                                          fsi.flow_size);
  refresh_expanded(f_ad_x_and_j_n_x_ext);
  Field<AdjointColorMatrix> mf;
  set_m_mat_plaq_mask_mu_no_comm(mf, nf, f_ad_x_and_j_n_x_ext, fsi.mask, fsi.mu,
                                 fsi.epsilon, fsi.flow_size);
  set_gm_force_from_flow_no_comm(gm_force_propagated, gm_force_pre_ext, mf,
                                 fsi.mask, fsi.mu, fsi.flow_size);
}

inline void set_gm_force_propagated_and_gm_force_det_from_flow_step(
    GaugeMomentum& gm_force_propagated, GaugeMomentum& gm_force_det,
    const GaugeMomentum& gm_force_pre, const GaugeField& gf0_ext,
    const FlowStepInfo& fsi)
//
// The gm_force_propagated and gm_force_det need to be summed together.
//
// The gf0_ext must be extended and refreshed and not flowed
// (for this flow step)
{
  TIMER("set_gm_force_propagated_and_gm_force_det_from_flow_step");
  qassert(is_initialized(gm_force_pre));
  qassert(is_initialized(gf0_ext));
  const Coordinate expand_left(fsi.flow_size, fsi.flow_size, fsi.flow_size,
                               fsi.flow_size);
  const Coordinate expand_right(fsi.flow_size, fsi.flow_size, fsi.flow_size,
                                fsi.flow_size);
  const Geometry geo = geo_reform(gf0_ext.geo(), 1);
  const Geometry geo_ext = geo_reform(geo, 1, expand_left, expand_right);
  GaugeMomentum gm_force_pre_ext;
  gm_force_pre_ext.init(geo_ext);
  gm_force_pre_ext = gm_force_pre;
  refresh_expanded(gm_force_pre_ext);
  // clear the output forces
  gm_force_det.init();
  gm_force_propagated.init();
  // To propagate force
  FieldM<ColorMatrix, 1> cf;
  set_flow_staple_mask_mu_no_comm(cf, gf0_ext, fsi.mask, fsi.mu, fsi.flow_size);
  Field<array<ColorMatrix, 2> > ducf;
  set_d_uc_mat_plaq_mask_mu_no_comm(ducf, cf, gf0_ext, fsi.mask, fsi.mu,
                                    fsi.flow_size);
  Field<AdjointColorMatrix> nf;
  set_n_mat_plaq_mask_mu_no_comm(nf, ducf, fsi.mask, fsi.flow_size);
  FieldM<AdjointColorMatrix, 2> f_ad_x_and_j_n_x_ext;
  f_ad_x_and_j_n_x_ext.init(geo_ext);
  set_ad_x_and_j_n_x_plaq_mask_mu_no_comm(f_ad_x_and_j_n_x_ext, cf, gf0_ext,
                                          fsi.mask, fsi.mu, fsi.epsilon,
                                          fsi.flow_size);
  refresh_expanded(f_ad_x_and_j_n_x_ext);
  Field<AdjointColorMatrix> mf;
  set_m_mat_plaq_mask_mu_no_comm(mf, nf, f_ad_x_and_j_n_x_ext, fsi.mask, fsi.mu,
                                 fsi.epsilon, fsi.flow_size);
  set_gm_force_from_flow_no_comm(gm_force_propagated, gm_force_pre_ext, mf,
                                 fsi.mask, fsi.mu, fsi.flow_size);
  // To include force from determinants below
  FieldM<AdjointColorMatrix, 1> mpf;
  set_mp_mat_plaq_mask_mu_no_comm(mpf, nf, cf, gf0_ext, fsi.mask, fsi.mu,
                                  fsi.epsilon, fsi.flow_size);
  FieldM<array<double, 8>, 1> f_e2_dj_x_n_mp_inv_ext;
  f_e2_dj_x_n_mp_inv_ext.init(geo_ext);
  FieldM<AdjointColorMatrix, 1> f_n_e_mp_inv_j_x_ext;
  f_n_e_mp_inv_j_x_ext.init(geo_ext);
  set_f_det_util_plaq_mask_mu(f_e2_dj_x_n_mp_inv_ext, f_n_e_mp_inv_j_x_ext, mpf,
                              nf, cf, gf0_ext, fsi.mask, fsi.mu, fsi.epsilon,
                              fsi.flow_size);
  refresh_expanded(f_n_e_mp_inv_j_x_ext);
  refresh_expanded(f_e2_dj_x_n_mp_inv_ext);
  set_gm_force_from_flow_det_no_comm(gm_force_det, f_e2_dj_x_n_mp_inv_ext,
                                     f_n_e_mp_inv_j_x_ext, nf, ducf, fsi.mask,
                                     fsi.mu, fsi.flow_size);
}

inline void set_flowed_gauge_fields(std::vector<GaugeField>& gf_ext_vec,
                                    const GaugeField& gf0, const FlowInfo& fi)
// Fill gf_ext_vec in the flow order (start with gf0)
// Note gf_ext_vec.size() == fi.v.size() + 1
//
// The gf_ext_vec.back() is the flowed (physical) gauge field.
//
// All gauge fields are refreshed except gf_ext_vec.back().
{
  TIMER("set_flowed_gauge_fields");
  const int n_steps = fi.v.size();
  clear(gf_ext_vec);
  gf_ext_vec.resize(n_steps + 1);
  const Coordinate expand_left(2, 2, 2, 2);
  const Coordinate expand_right(2, 2, 2, 2);
  const Geometry geo_ext = geo_reform(gf0.geo(), 1, expand_left, expand_right);
  gf_ext_vec[0].init(geo_ext);
  gf_ext_vec[0] = gf0;
  for (int i = 0; i < (int)fi.v.size(); ++i) {
    const FlowStepInfo& fsi = fi.v[i];
    refresh_expanded(gf_ext_vec[i]);
    gf_ext_vec[i + 1].init(geo_ext);
    gf_flow_plaq_mask_mu_no_comm(gf_ext_vec[i + 1], gf_ext_vec[i], fsi.mask,
                                 fsi.mu, fsi.epsilon, fsi.flow_size);
  }
}

inline void set_gm_force_propagated_det_from_flow(
    GaugeMomentum& gm_force, const GaugeMomentum& gm_force_pre,
    const std::vector<GaugeField>& gf_ext_vec, const FlowInfo& fi)
//
// Propagate the force act on the flowed (physical) gauge field to the unflowed
// gauge field.
//
// Force from the determinant is included.
//
// Call set_flowed_gauge_fields(gf_ext_vec, gf0, fi) to obtain gf_ext_vec.
{
  TIMER("set_gm_force_propagated_det_from_flow");
  const Geometry geo = geo_reform(gm_force_pre.geo());
  gm_force.init(geo);
  gm_force = gm_force_pre;
  for (int i = fi.v.size() - 1; i >= 0; --i) {
    GaugeMomentum gm_force_det;
    set_gm_force_propagated_and_gm_force_det_from_flow_step(
        gm_force, gm_force_det, gm_force, gf_ext_vec[i], fi.v[i]);
    gm_force += gm_force_det;
  }
}

inline void set_gm_force_propagated_no_det_from_flow(
    GaugeMomentum& gm_force, const GaugeMomentum& gm_force_pre,
    const std::vector<GaugeField>& gf_ext_vec, const FlowInfo& fi)
//
// Propagate the force act on the flowed (physical) gauge field to the unflowed
// gauge field.
//
// Force from the determinant is included.
//
// Call set_flowed_gauge_fields(gf_ext_vec, gf0, fi) to obtain gf_ext_vec.
{
  TIMER("set_gm_force_propagated_no_det_from_flow");
  const Geometry geo = geo_reform(gm_force_pre.geo());
  gm_force.init(geo);
  gm_force = gm_force_pre;
  for (int i = fi.v.size() - 1; i >= 0; --i) {
    set_gm_force_propagated_from_flow_step(gm_force, gm_force, gf_ext_vec[i],
                                           fi.v[i]);
  }
}

inline void set_gm_force_flowed(GaugeMomentum& gm_force, const GaugeField& gf0,
                                const GaugeAction& ga, const FlowInfo& fi)
{
  TIMER("set_gm_force_flowed");
  std::vector<GaugeField> gf_ext_vec;
  set_flowed_gauge_fields(gf_ext_vec, gf0, fi);
  gm_force.init();
  set_gm_force(gm_force, gf_ext_vec.back(), ga);
  set_gm_force_propagated_det_from_flow(gm_force, gm_force, gf_ext_vec, fi);
}

inline void set_gm_force_flowed_no_det(GaugeMomentum& gm_force,
                                       const GaugeMomentum& gm_force_pre,
                                       const GaugeField& gf0,
                                       const FlowInfo& fi)
{
  TIMER("set_gm_force_flowed_no_det");
  std::vector<GaugeField> gf_ext_vec;
  set_flowed_gauge_fields(gf_ext_vec, gf0, fi);
  set_gm_force_propagated_no_det_from_flow(gm_force, gm_force_pre, gf_ext_vec,
                                           fi);
}

}  // namespace qlat
