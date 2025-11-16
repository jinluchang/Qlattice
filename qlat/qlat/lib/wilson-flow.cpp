#include <qlat-utils/types.h>
#include <qlat/core.h>
#include <qlat/field-shuffle.h>
#include <qlat/qcd-acc.h>
#include <qlat/wilson-flow.h>

namespace qlat
{  //

void set_wilson_flow_z(GaugeMomentum& z, const GaugeField& gf, const RealD c1)
{
  TIMER("set_wilson_flow_z");
  const GaugeAction ga(3.0, c1);
  set_gm_force(z, gf, ga);
}

void gf_wilson_flow_step_euler(GaugeField& gf, const RealD epsilon,
                               const RealD c1)
{
  TIMER("gf_wilson_flow_step_euler");
  GaugeField& w = gf;
  GaugeMomentum z;
  set_wilson_flow_z(z, w, c1);
  gf_evolve(w, z, epsilon);
}

void gf_wilson_flow_step(GaugeField& gf, const RealD epsilon, const RealD c1)
// Runge-Kutta scheme
// http://arxiv.org/abs/1006.4518v3
{
  TIMER("gf_wilson_flow_step");
  GaugeField& w = gf;
  GaugeMomentum z, zp;
  set_wilson_flow_z(z, w, c1);
  z *= 1.0 / 4.0;
  gf_evolve(w, z, epsilon);
  qswap(z, zp);
  zp *= 17.0 / 9.0;
  set_wilson_flow_z(z, w, c1);
  z *= 8.0 / 9.0;
  z -= zp;
  gf_evolve(w, z, epsilon);
  qswap(z, zp);
  set_wilson_flow_z(z, w, c1);
  z *= 3.0 / 4.0;
  z -= zp;
  gf_evolve(w, z, epsilon);
}

// --------------------

static qacc RealD gf_energy_density_dir_site_no_comm(const GaugeField& gf,
                                                      const Coordinate& xl,
                                                      const Int mu,
                                                      const Int nu)
{
  const ColorMatrix g_mu_nu =
      make_tr_less_anti_herm_matrix(gf_clover_leaf_no_comm(gf, xl, mu, nu));
  const RealD s = -matrix_trace(g_mu_nu, g_mu_nu).real();
  return s;
}

static void gf_energy_density_dir_field_no_comm(Field<RealD>& fd,
                                                const GaugeField& gf)
{
  TIMER("gf_energy_density_dir_field_no_comm");
  const Geometry geo = geo_resize(gf.geo());
  qassert(geo.is_only_local);
  fd.init(geo, 6);
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<RealD> v = fd.get_elems(index);
    v[0] = gf_energy_density_dir_site_no_comm(gf, xl, 0, 1);
    v[1] = gf_energy_density_dir_site_no_comm(gf, xl, 0, 2);
    v[2] = gf_energy_density_dir_site_no_comm(gf, xl, 0, 3);
    v[3] = gf_energy_density_dir_site_no_comm(gf, xl, 1, 2);
    v[4] = gf_energy_density_dir_site_no_comm(gf, xl, 1, 3);
    v[5] = gf_energy_density_dir_site_no_comm(gf, xl, 2, 3);
  });
}

void gf_energy_density_dir_field(Field<RealD>& fd, const GaugeField& gf)
// Similar to `gf_plaq_feild`
// fd.init(geo, 6);
// https://arxiv.org/pdf/1006.4518.pdf Eq. (2.1) (Fig. 1) (approximate Eq. (3.1))
// https://arxiv.org/pdf/1203.4469.pdf
{
  TIMER("gf_energy_density_field");
  GaugeField gf1;
  gf1.init(geo_resize(gf.geo(), 1));
  gf1 = gf;
  refresh_expanded(gf1);
  gf_energy_density_dir_field_no_comm(fd, gf1);
}

// --------------------

static void gf_energy_density_field_no_comm(Field<RealD>& fd,
                                            const GaugeField& gf)
{
  TIMER("gf_energy_density_field_no_comm");
  const Geometry geo = geo_resize(gf.geo());
  fd.init(geo, 1);
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    RealD s = 0.0;
    for (Int mu = 0; mu < 3; ++mu) {
      for (Int nu = mu + 1; nu < 4; ++nu) {
        const ColorMatrix g_mu_nu = make_tr_less_anti_herm_matrix(
            gf_clover_leaf_no_comm(gf, xl, mu, nu));
        s += -matrix_trace(g_mu_nu, g_mu_nu).real();
      }
    }
    fd.get_elem(index) = s;
  });
}

void gf_energy_density_field(Field<RealD>& fd, const GaugeField& gf)
// https://arxiv.org/pdf/1006.4518.pdf Eq. (2.1) (Fig. 1) (approximate Eq. (3.1))
// https://arxiv.org/pdf/1203.4469.pdf
{
  TIMER("gf_energy_density_field");
  GaugeField gf1;
  gf1.init(geo_resize(gf.geo(), 1));
  gf1 = gf;
  refresh_expanded(gf1);
  gf_energy_density_field_no_comm(fd, gf1);
}

RealD gf_energy_density(const GaugeField& gf)
// https://arxiv.org/pdf/1006.4518.pdf Eq. (2.1) (Fig. 1) (approximate Eq. (3.1))
// https://arxiv.org/pdf/1203.4469.pdf
{
  TIMER("gf_energy_density");
  const Geometry& geo = gf.geo();
  FieldM<RealD, 1> fd;
  gf_energy_density_field(fd, gf);
  return field_glb_sum(fd)[0] / (RealD)geo.total_volume();
}

// --------------------

std::vector<RealD> gf_wilson_flow(GaugeField& gf,
                                   const RealD existing_flow_time,
                                   const RealD flow_time, const Int steps,
                                   const RealD c1)
{
  TIMER("gf_wilson_flow");
  std::vector<RealD> energy_density_list(steps, 0.0);
  const RealD epsilon = flow_time / (RealD)steps;
  for (Int i = 0; i < steps; ++i) {
    gf_wilson_flow_step(gf, epsilon, c1);
    const RealD t = (i + 1) * epsilon + existing_flow_time;
    const RealD energy_density = gf_energy_density(gf);
    energy_density_list[i] = energy_density;
    displayln_info(fname +
                   ssprintf(": t = %24.17E ; E = %24.17E ; t^2 E = %24.17E.", t,
                            energy_density, sqr(t) * energy_density));
  }
  return energy_density_list;
}

static qacc ColorMatrix gf_plaq_flow_staple_no_comm(
    const GaugeField& gf, const Field<RealD>& plaq_factor, const Coordinate& xl,
    const Int mu)
// transpose the same way as gf.get_elem(xl, mu)
{
  ColorMatrix acc;
  set_zero(acc);
  array<Int, 16> table;
  table.fill(-1);
  table[0 * 4 + 1] = 0;
  table[0 * 4 + 2] = 1;
  table[0 * 4 + 3] = 2;
  table[1 * 4 + 2] = 3;
  table[1 * 4 + 3] = 4;
  table[2 * 4 + 3] = 5;
  table[1 * 4 + 0] = 0;
  table[2 * 4 + 0] = 1;
  table[3 * 4 + 0] = 2;
  table[2 * 4 + 1] = 3;
  table[3 * 4 + 1] = 4;
  table[3 * 4 + 2] = 5;
  const Vector<RealD> vpf = plaq_factor.get_elems_const(xl);
  for (Int nu = 0; nu < 4; ++nu) {
    if (nu == mu) {
      continue;
    }
    const Vector<RealD> vpf_nu =
        plaq_factor.get_elems_const(coordinate_shifts(xl, -nu - 1));
    const Int m = table[mu * 4 + nu];
    qassert(0 <= m and m < 6);
    acc += vpf[m] *
           gf_wilson_line_no_comm(gf, xl, make_array<int>(nu, mu, -nu - 1));
    acc += vpf_nu[m] *
           gf_wilson_line_no_comm(gf, xl, make_array<int>(-nu - 1, mu, nu));
  }
  return acc;
}

static qacc ColorMatrix
gf_plaq_flow_site_no_comm(const GaugeField& gf, const Field<RealD>& plaq_factor,
                          const Coordinate& xl, const Int mu)
{
  const ColorMatrix ad_staple =
      matrix_adjoint(gf_plaq_flow_staple_no_comm(gf, plaq_factor, xl, mu));
  const ColorMatrix force = -gf.get_elem(xl, mu) * ad_staple;
  return make_tr_less_anti_herm_matrix(force);
}

static void set_plaq_flow_z_no_comm(GaugeMomentum& z,
                                    const GaugeField& gf,
                                    const Field<RealD>& plaq_factor)
// gf need comm
{
  TIMER("set_plaq_flow_z_no_comm");
  const Geometry geo = geo_resize(gf.get_geo());
  z.init(geo);
  qassert(gf.multiplicity == 4);
  qassert(z.multiplicity == 4);
  qassert(plaq_factor.multiplicity == 6);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = z.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> gm_force_v = z.get_elems(xl);
    qassert(gm_force_v.size() == 4);
    for (Int mu = 0; mu < 4; ++mu) {
      gm_force_v[mu] = gf_plaq_flow_site_no_comm(gf, plaq_factor, xl, mu);
    }
  });
}

void set_plaq_flow_z(GaugeMomentum& z, const GaugeField& gf,
                     const Field<RealD>& plaq_factor)
// Compute force with plaq dependent beta factor (relative to standard
// `set_wilson_flow_z`).
// `plaq_factor.multiplicity == 6`.
// Check `gf_plaq_field` for the order of plaq.
{
  TIMER("set_plaq_flow_z");
  qassert(gf.multiplicity == 4);
  qassert(plaq_factor.multiplicity == 6);
  const Coordinate expand_left = Coordinate(1, 1, 1, 1);
  const Coordinate expand_right = Coordinate(1, 1, 1, 1);
  const Geometry geo_ext = geo_resize(gf.geo(), expand_left, expand_right);
  GaugeField gf_ext;
  gf_ext.init(geo_ext);
  qassert(gf_ext.multiplicity == 4);
  gf_ext = gf;
  const std::string tag_comm = "plaq";
  QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
  const CommPlan& plan = get_comm_plan(set_marks_field_gm_force, tag_comm,
                                       gf_ext.geo(), gf_ext.multiplicity);
  QLAT_DIAGNOSTIC_POP;
  refresh_expanded(gf_ext, plan);
  const Geometry geo_pf_ext = geo_resize(gf.geo(), expand_left, Coordinate());
  Field<RealD> plaq_factor_ext;
  plaq_factor_ext.init(geo_pf_ext, 6);
  plaq_factor_ext = plaq_factor;
  refresh_expanded(plaq_factor_ext);
  set_plaq_flow_z_no_comm(z, gf_ext, plaq_factor_ext);
  qassert(z.multiplicity == 4);
}

static void set_marks_field_block_stout_smear(CommMarks& marks,
                                              const Geometry& geo,
                                              const Int multiplicity,
                                              const std::string& tag)
{
  TIMER_VERBOSE("set_marks_field_block_stout_smear");
  Qassert(multiplicity == 4);
  marks.init();
  marks.init(geo, multiplicity);
  set_zero(marks);
  const Coordinate block_site = read_coordinate(tag);
  const bool is_block_site = block_site != Coordinate();
  // const Coordinate size_node = geo.geon.size_node;
  // const Coordinate coor_node = geo.geon.coor_node;
  // const Coordinate node_site = geo.node_site;
  const Coordinate total_site = geo.total_site();
  Coordinate size_block;
  if (is_block_site) {
    size_block = total_site / block_site;
    Qassert(total_site == block_site * size_block);
  }
  qthread_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    Coordinate xl_block;
    if (is_block_site) {
      xl_block = xg % block_site;
    }
    for (Int mu = 0; mu < 4; ++mu) {
      if (is_block_site and (xl_block[mu] == block_site[mu] - 1)) {
        continue;
      }
      for (Int nu = -4; nu < 4; ++nu) {
        if (nu == mu or -nu - 1 == mu) {
          continue;
        }
        if (is_block_site and (nu >= 0 and xl_block[nu] == block_site[nu] - 1)) {
          continue;
        }
        if (is_block_site and (nu < 0 and xl_block[-nu - 1] == 0)) {
          continue;
        }
        set_marks_field_path(marks, xl, make_array<int>(nu, mu, -nu - 1));
      }
    }
  });
}

void gf_block_stout_smear(GaugeField& gf, const GaugeField& gf0,
                          const Coordinate& block_site, const RealD step_size)
// gf can be the same object as gf0
{
  TIMER("gf_block_stout_smear");
  Qassert(gf0.multiplicity == 4);
  const Geometry geo = gf0.geo();
  Qassert(geo.is_only_local);
  const bool is_block_site = block_site != Coordinate();
  // const Coordinate size_node = geo.geon.size_node;
  const Coordinate coor_node = geo.geon.coor_node;
  const Coordinate node_site = geo.node_site;
  const Coordinate total_site = geo.total_site();
  Coordinate size_block;
  Coordinate expand_left = Coordinate(1, 1, 1, 1);
  Coordinate expand_right = Coordinate(1, 1, 1, 1);
  if (is_block_site) {
    size_block = total_site / block_site;
    Qassert(total_site == block_site * size_block);
    for (Int i = 0; i < 4; ++i) {
      if ((coor_node[i] * node_site[i]) % block_site[i] == 0) {
        expand_left[i] = 0;
      }
      if (((coor_node[i] + 1) * node_site[i]) % block_site[i] == 0) {
        expand_right[i] = 0;
      }
    }
  }
  const Geometry geo_ext = geo_resize(geo, expand_left, expand_right);
  GaugeField gf_ext;
  gf_ext.init(geo_ext);
  Qassert(gf_ext.multiplicity == 4);
  gf_ext = gf0;
  const std::string tag_comm = show(block_site);
  QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
  const CommPlan& plan =
      get_comm_plan(set_marks_field_block_stout_smear, tag_comm, gf_ext.geo(),
                    gf_ext.multiplicity);
  QLAT_DIAGNOSTIC_POP;
  refresh_expanded(gf_ext, plan);
  GaugeMomentum gm;
  gm.init(geo);
  set_zero(gm);
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    Coordinate xl_block;
    if (is_block_site) {
      xl_block = xg % block_site;
    }
    for (Int mu = 0; mu < 4; ++mu) {
      if (is_block_site and (xl_block[mu] == block_site[mu] - 1)) {
        continue;
      }
      ColorMatrix acc;
      set_zero(acc);
      for (Int nu = -4; nu < 4; ++nu) {
        if (nu == mu or -nu - 1 == mu) {
          continue;
        }
        if (is_block_site and (nu >= 0 and xl_block[nu] == block_site[nu] - 1)) {
          continue;
        }
        if (is_block_site and (nu < 0 and xl_block[-nu - 1] == 0)) {
          continue;
        }
        acc += gf_wilson_line_no_comm(gf_ext, xl,
                                      make_array<int>(nu, mu, -nu - 1));
      }
      const ColorMatrix force = -gf_ext.get_elem(xl, mu) * matrix_adjoint(acc);
      gm.get_elem(xl, mu) = make_tr_less_anti_herm_matrix(force);
    }
  });
  gf = gf0;
  gf_evolve(gf, gm, step_size);
}

void gf_local_stout_smear(GaugeField& gf, const GaugeField& gf0,
                          const RealD step_size)
// gf can be the same object as gf0
// perform stout smear on each local object separately
// intended to be used after `shuffle_field(gf_vec, gf, new_size_node)`
{
  TIMER("gf_local_stout_smear");
  Qassert(gf0.multiplicity == 4);
  const Geometry geo = gf0.geo();
  Qassert(geo.is_only_local);
  const Coordinate node_site = geo.node_site;
  GaugeMomentum gm;
  gm.init(geo);
  set_zero(gm);
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    for (Int mu = 0; mu < 4; ++mu) {
      if (xl[mu] == node_site[mu] - 1) {
        continue;
      }
      ColorMatrix acc;
      set_zero(acc);
      for (Int nu = -4; nu < 4; ++nu) {
        if (nu == mu or -nu - 1 == mu) {
          continue;
        }
        if (nu >= 0 and xl[nu] == node_site[nu] - 1) {
          continue;
        }
        if (nu < 0 and xl[-nu - 1] == 0) {
          continue;
        }
        acc +=
            gf_wilson_line_no_comm(gf0, xl, make_array<int>(nu, mu, -nu - 1));
      }
      const ColorMatrix force = -gf0.get_elem(xl, mu) * matrix_adjoint(acc);
      gm.get_elem(xl, mu) = make_tr_less_anti_herm_matrix(force);
    }
  });
  gf = gf0;
  gf_evolve(gf, gm, step_size);
}

void set_local_tree_gauge_f_dir(Field<Int>& f_dir, const Geometry& geo,
                                const RngState& rs)
{
  TIMER("set_local_tree_gauge_f_dir");
  Qassert(geo.is_only_local);
  const Coordinate node_site = geo.node_site;
  const Coordinate node_center = node_site / 2;
  const Coordinate total_site = geo.total_site();
  f_dir.init(geo, 1);
  set_zero(f_dir);
  qthread_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    if (xl == node_center) {
      f_dir.get_elem(index) = 0;
      continue;
    }
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Long gindex = index_from_coordinate(xg, total_site);
    RngState rs_local = rs.newtype(gindex);
    while (true) {
      const Int dir = rand_gen(rs_local) % 4;
      if (xl[dir] > node_center[dir]) {
        f_dir.get_elem(index) = -dir - 1;
      } else if (xl[dir] < node_center[dir]) {
        f_dir.get_elem(index) = dir;
      } else {
        continue;
      }
      break;
    }
  });
}

void gt_local_tree_gauge(GaugeTransform& gt_inv, const GaugeField& gf,
                         const Field<Int>& f_dir)
// `f_dir` can be set by:
// set_local_tree_gauge_f_dir(f_dir, geo, rs);
// Note that to apply the local tree gauge, we need to invert the returned `gt` first.
// E.g.,
// gt_invert(gt, gt_inv);
// gf_apply_gauge_transformation(gf, gf, gt);
{
  TIMER("gt_local_tree_gauge");
  Qassert(gf.multiplicity == 4);
  Qassert(f_dir.multiplicity == 1);
  const Geometry geo = gf.geo();
  gt_inv.init(geo);
  set_unit(gt_inv);
  Qassert(geo.is_only_local);
  Field<Int> f_marks, f_marks_new;
  f_marks.init(geo, 1);
  f_marks_new.init(geo, 1);
  set_zero(f_marks);
  set_zero(f_marks_new);
  const Coordinate node_site = geo.node_site;
  const Coordinate node_center = node_site / 2;
  const Int dir = f_dir.get_elem(0);
  Int num_step = 0;
  for (Int mu = 0; mu < 4; ++mu) {
    num_step += node_site[mu] / 2;
  }
  for (Int i = 0; i < num_step + 1; ++i) {
    qacc_for(index, geo.local_volume(), {
      const Coordinate xl = geo.coordinate_from_index(index);
      const Int mark = f_marks.get_elem(index);
      if (mark != 0) {
        continue;
      }
      if (xl == node_center) {
        set_unit(gt_inv.get_elem(index));
        f_marks_new.get_elem(index) = 1;
        continue;
      }
      const Int dir = f_dir.get_elem(index);
      const Coordinate xl1 = coordinate_shifts(xl, dir);
      const Int mark1 = f_marks.get_elem(xl1);
      if (mark1 == 0) {
        continue;
      }
      gt_inv.get_elem(index) = gf_get_link(gf, xl, dir) * gt_inv.get_elem(xl1);
      f_marks_new.get_elem(index) = 1;
    });
    qacc_for(index, geo.local_volume(), {
      f_marks.get_elem(index) += f_marks_new.get_elem(index);
    });
    if (i == num_step) {
      qacc_for(index, geo.local_volume(),
               { qassert(f_marks.get_elem(index) >= 1); });
    }
  }
}

}  // namespace qlat
