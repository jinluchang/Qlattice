#pragma once

#include "hlbl-sl-table.h"
#include "muon-line.h"

namespace qlat
{  //

inline void set_m_z_field_tag(SelectedField<RealD>& smf_d,
                              const FieldSelection& fsel,
                              const Coordinate& xg_x, const Coordinate& xg_y,
                              const double a, const int tag)
// interface
// tag = 0 sub
// tag = 1 nosub
{
  TIMER_VERBOSE("set_m_z_field_tag(smf_d,fsel,xg_x,xg_y,a,tag)");
  SelectedField<ManyMagneticMoments>& smf =
      qcast<ManyMagneticMoments, RealD>(smf_d);
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate total_site = geo.total_site();
  smf.init(fsel, 1);
  qthread_for(idx, fsel.n_elems, {
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);  // z location
    ManyMagneticMoments& mmm = smf.get_elem(idx);
    mmm = get_muon_line_m_extra_lat(xg_x, xg_y, xg, total_site, a, tag);
  });
  qcast<RealD, ManyMagneticMoments>(smf);
}

// ------------------------------------------------------------------------

qacc ManyMagneticMoments simple_pion_projection(const CoordinateD& x,
                                                const CoordinateD& y,
                                                const CoordinateD& z)
{
  const CoordinateD mid_yz = (y + z) / 2;
  const CoordinateD x_mid_yz = x - mid_yz;
  const CoordinateD y_z = y - z;
  array<array<double, 4>, 3> eps_eps_xyz;
  set_zero(eps_eps_xyz);
  for (int i = 0; i < 3; ++i) {
    for (int rho = 0; rho < 4; ++rho) {
      for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < 3; ++k) {
          for (int m = 0; m < 4; ++m) {
            eps_eps_xyz[i][rho] += epsilon_tensor_acc(i, j, k) *
                                   epsilon_tensor_acc(rho, k, j, m) *
                                   x_mid_yz[m];
          }
        }
      }
    }
  }
  array<array<double, 4>, 4> eps_yz_xyz;
  set_zero(eps_yz_xyz);
  for (int sigma = 0; sigma < 4; ++sigma) {
    for (int lambda = 0; lambda < 4; ++lambda) {
      for (int r = 0; r < 4; ++r) {
        for (int s = 0; s < 4; ++s) {
          eps_yz_xyz[sigma][lambda] +=
              epsilon_tensor_acc(sigma, lambda, r, s) * y_z[r] * x_mid_yz[s];
        }
      }
    }
  }
  ManyMagneticMoments ret;
  set_zero(ret);
  for (int i = 0; i < 3; ++i) {
    for (int rho = 0; rho < 4; ++rho) {
      for (int sigma = 0; sigma < 4; ++sigma) {
        for (int lambda = 0; lambda < 4; ++lambda) {
          get_m_comp(ret, i, rho, sigma, lambda) =
              0.5 * eps_eps_xyz[i][rho] * eps_yz_xyz[sigma][lambda];
        }
      }
    }
  }
  return ret;
}

qacc ManyMagneticMoments pion_projection(const CoordinateD& x,
                                         const CoordinateD& y,
                                         const CoordinateD& z)
// y z are close
// x and y z are two ends of the pion
{
  const ManyMagneticMoments mmm = simple_pion_projection(x, y, z);
  double sum = 0.0;
  for (int i = 0; i < 3; ++i) {
    for (int rho = 0; rho < 4; ++rho) {
      for (int sigma = 0; sigma < 4; ++sigma) {
        for (int lambda = 0; lambda < 4; ++lambda) {
          sum += sqr(get_m_comp(mmm, i, rho, sigma, lambda));
        }
      }
    }
  }
  if (sum == 0.0) {
    return mmm;
  } else {
    return (1.0 / sqrt(sum)) * mmm;
  }
}

qacc ManyMagneticMoments pion_projection(const Coordinate& x,
                                         const Coordinate& y,
                                         const Coordinate& z,
                                         const Coordinate& total_site,
                                         const bool is_permute)
// if is_permute == false,
// y z should be located on one end of the pion,
// x is on the other end.
{
  ManyMagneticMoments mmm;
  set_zero(mmm);
  const Coordinate xy = smod(x - y, total_site);
  const Coordinate xz = smod(x - z, total_site);
  const Coordinate zz = smod(z - z, total_site);
  const Coordinate yz = smod(y - z, total_site);
  if (xy != xz - yz) {
    return mmm;
  }
  const long xy_len = sqr(xy);
  const long xz_len = sqr(xz);
  const long yz_len = sqr(yz);
  if (is_permute) {
    if (yz_len <= xz_len and yz_len <= xy_len) {
      return pion_projection(xz, yz, zz);
    } else if (xy_len <= xz_len and xy_len <= yz_len) {
      return permute_rho_sigma_nu(pion_projection(zz, xz, yz), 1, 2, 0);
    } else if (xz_len <= xy_len and xz_len <= yz_len) {
      return permute_rho_sigma_nu(pion_projection(yz, zz, xz), 2, 0, 1);
    }
  } else {
    return pion_projection(xz, yz, zz);
  }
  qassert(false);
  return mmm;
}

// ------------------------------------------------------------------------

inline void set_local_current_from_props(
    SelectedField<WilsonMatrix>& scf, const SelectedField<WilsonMatrix>& sprop1,
    const SelectedField<WilsonMatrix>& sprop2, const FieldSelection& fsel)
// ->- sprop1 ->- gamma_mu ->- gamma5 sprop2^+ gamma5 ->-
{
  TIMER_VERBOSE("set_local_current_from_props");
  const Geometry& geo = fsel.f_rank.geo();
  qassert(geo == sprop1.geo());
  qassert(geo == sprop2.geo());
  const array<SpinMatrix, 4>& gammas = SpinMatrixConstants::get_cps_gammas();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  scf.init(fsel, 4);
  set_zero(scf);
  qacc_for(idx, fsel.n_elems, {
    const WilsonMatrix& m1 = sprop1.get_elem(idx);
    const WilsonMatrix& m2 = sprop2.get_elem(idx);
    Vector<WilsonMatrix> v = scf.get_elems(idx);
    const WilsonMatrix m2rev =
        gamma5 * (WilsonMatrix)matrix_adjoint(m2) * gamma5;
    for (int m = 0; m < 4; ++m) {
      v[m] = m2rev * gammas[m] * m1;
    }
  });
}

template <class M>
struct CurrentMoments {
  vector_acc<array<M, 3 * 3>> d;
  // moment = 0.5 * epsilon_tensor(i, j, k) * smod_sym(xg[j] - ref[j],
  // total_site[j]) * d[ xg[j] ][3*j + k]
  //
  void init() { clear(d); }
  void init(const int lsize)
  {
    init();
    d.resize(lsize);
    set_zero(d);
  }
};

template <class M>
void set_current_moments_from_current(CurrentMoments<M>& cm,
                                      const SelectedField<M>& current,
                                      const FieldSelection& fsel)
{
  TIMER("set_current_moments_from_current");
  const Geometry& geo = current.geo();
  qassert(geo.multiplicity == 4);
  const Coordinate total_site = geo.total_site();
  const int lsize =
      std::max(total_site[0], std::max(total_site[1], total_site[2]));
  cm.init(lsize);
  qfor(idx, fsel.n_elems, {
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Vector<M> v = current.get_elems_const(idx);
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        cm.d[xg[j]][3 * j + k] += v[k];
      }
    }
  });
  glb_sum_double_vec(get_data(cm.d));
}

template <class M>
qacc array<M, 3> simple_moment(const CurrentMoments<M>& cm,
                               const CoordinateD& ref,
                               const Coordinate& total_site)
{
  const int lsize =
      std::max(total_site[0], std::max(total_site[1], total_site[2]));
  array<M, 3> ret;
  set_zero(ret);
  for (int x = 0; x < lsize; ++x) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        if (i == j) {
          continue;
        }
        for (int k = 0; k < 3; ++k) {
          if (i == k or j == k) {
            continue;
          }
          ret[i] += (Complex)(0.5 * epsilon_tensor_acc(i, j, k) *
                              smod_sym(x - ref[j], (double)total_site[j])) *
                    cm.d[x][3 * j + k];
        }
      }
    }
  }
  return ret;
}

template <class M>
qacc array<M, 3> simple_moment_with_contact_subtract(
    const CurrentMoments<M>& cm, const CoordinateD& ref,
    const Coordinate& total_site, const SelectedField<M>& current,
    const FieldSelection& fsel, const Coordinate& xg, const Complex& coef)
// xg should refer to a local site
// xg should be a selected point
// coef = (prob - 1) where prob is the probability of a point being selected,
// e.g. prob = 1.0/16.0, coef = -15.0/16.0
{
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate xl = geo.coordinate_l_from_g(xg);
  qassert(geo.is_local(xl));
  qassert(fsel.f_rank.get_elem(xl) >= 0);
  const long idx = fsel.f_local_idx.get_elem(xl);
  const Vector<M> cv = current.get_elems_const(idx);
  array<M, 3> ret = simple_moment(cm, ref, total_site);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (i == j) {
        continue;
      }
      for (int k = 0; k < 3; ++k) {
        if (i == k or j == k) {
          continue;
        }
        ret[i] += (coef * 0.5 * (Complex)epsilon_tensor_acc(i, j, k) *
                   (Complex)smod_sym(xg[j] - ref[j], (double)total_site[j])) *
                  cv[k];
      }
    }
  }
  return ret;
}

enum ChooseReferenceLabel {
  choose_reference_label_ref_far,
  choose_reference_label_ref_close,
  choose_reference_label_ref_center,
};

inline ChooseReferenceLabel choose_reference_label(const std::string& label)
{
  if (does_string_have_tag(label, "ref-far")) {
    return choose_reference_label_ref_far;
  } else if (does_string_have_tag(label, "ref-close")) {
    return choose_reference_label_ref_close;
  } else if (does_string_have_tag(label, "ref-center")) {
    return choose_reference_label_ref_center;
  } else {
    qassert(false);
  }
  qassert(false);
}

qacc CoordinateD choose_reference(const Coordinate& xg_x,
                                  const Coordinate& xg_y,
                                  const Coordinate& xg_z,
                                  const Coordinate& total_site,
                                  const ChooseReferenceLabel& label)
{
  const long dis2_xy = sqr(smod(xg_x - xg_y, total_site));
  const long dis2_xz = sqr(smod(xg_x - xg_z, total_site));
  const long dis2_yz = sqr(smod(xg_y - xg_z, total_site));
  if (choose_reference_label_ref_far == label) {
    if (dis2_xy < dis2_xz and dis2_xy < dis2_yz) {
      return CoordinateD(xg_z);
    } else if (dis2_xz < dis2_xy and dis2_xz < dis2_yz) {
      return CoordinateD(xg_y);
    } else if (dis2_yz < dis2_xy and dis2_yz < dis2_xz) {
      return CoordinateD(xg_x);
    } else {
      return mod(
          CoordinateD(xg_z) + 1.0 / 3.0 *
                                  CoordinateD(smod(xg_x - xg_z, total_site) +
                                              smod(xg_y - xg_z, total_site)),
          total_site);
    }
  } else if (choose_reference_label_ref_close == label) {
    if (dis2_xy < dis2_xz and dis2_xy < dis2_yz) {
      return middle_mod(CoordinateD(xg_x), CoordinateD(xg_y),
                        CoordinateD(total_site));
    } else if (dis2_xz < dis2_xy and dis2_xz < dis2_yz) {
      return middle_mod(CoordinateD(xg_x), CoordinateD(xg_z),
                        CoordinateD(total_site));
    } else if (dis2_yz < dis2_xy and dis2_yz < dis2_xz) {
      return middle_mod(CoordinateD(xg_y), CoordinateD(xg_z),
                        CoordinateD(total_site));
    } else {
      return mod(
          CoordinateD(xg_z) + 1.0 / 3.0 *
                                  CoordinateD(smod(xg_x - xg_z, total_site) +
                                              smod(xg_y - xg_z, total_site)),
          total_site);
    }
  } else if (choose_reference_label_ref_center == label) {
    return mod(
        CoordinateD(xg_z) + 1.0 / 3.0 *
                                CoordinateD(smod(xg_x - xg_z, total_site) +
                                            smod(xg_y - xg_z, total_site)),
        total_site);
  } else {
    qassert(false);
  }
  qassert(false);
  return CoordinateD();
}

inline void contract_four_loop(
    SelectedField<Complex>& f_loop_i_rho_sigma_lambda, const Complex& coef,
    const Coordinate& xg_x, const Coordinate& xg_y,
    const SelectedField<WilsonMatrix>& c_xy,
    const SelectedField<WilsonMatrix>& c_yx,
    const CurrentMoments<WilsonMatrix>& cm_xy,
    const CurrentMoments<WilsonMatrix>& cm_yx, const FieldSelection& fsel,
    const std::string& label)
{
  TIMER_VERBOSE("contract_four_loop");
  const box<SpinMatrixConstantsT<>>& smc = get_spin_matrix_constants();
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate total_site = geo.total_site();
  qassert(geo == geo_remult(c_yx.geo()));
  qassert(geo == geo_remult(c_xy.geo()));
  f_loop_i_rho_sigma_lambda.init(fsel, 3 * 4 * 4 * 4);
  set_zero(f_loop_i_rho_sigma_lambda);
  // TODO:
  const double sel_reduction_ratio = 32.0;
  const Complex sub_coef = -(sel_reduction_ratio - 1.0) / sel_reduction_ratio;
  const Complex sel_coef = sel_reduction_ratio * sel_reduction_ratio;
  const Complex final_coef = coef * sel_coef;
  const ChooseReferenceLabel cr_label = choose_reference_label(label);
  qacc_for(idx, fsel.n_elems, {
    const array<SpinMatrix, 4>& gammas = smc().cps_gammas;
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const CoordinateD xgref =
        choose_reference(xg_x, xg_y, xg, total_site, cr_label);
    const array<WilsonMatrix, 3> sm_yx = simple_moment_with_contact_subtract(
        cm_yx, xgref, total_site, c_yx, fsel, xg, sub_coef);
    const array<WilsonMatrix, 3> sm_xy = simple_moment_with_contact_subtract(
        cm_xy, xgref, total_site, c_xy, fsel, xg, sub_coef);
    const Vector<WilsonMatrix> vc_xy = c_xy.get_elems_const(idx);
    const Vector<WilsonMatrix> vc_yx = c_yx.get_elems_const(idx);
    Vector<Complex> v_loop = f_loop_i_rho_sigma_lambda.get_elems(idx);
    for (int i = 0; i < 3; ++i) {
      for (int rho = 0; rho < 4; ++rho) {
        const WilsonMatrix wm_i_rho = sm_yx[i] * gammas[rho];
        const WilsonMatrix wm_rho_i = gammas[rho] * sm_xy[i];
        for (int sigma = 0; sigma < 4; ++sigma) {
          const WilsonMatrix wm_sigma_i_rho = gammas[sigma] * wm_i_rho;
          const WilsonMatrix wm_rho_i_sigma = wm_rho_i * gammas[sigma];
          for (int lambda = 0; lambda < 4; ++lambda) {
            v_loop[64 * i + 16 * rho + 4 * sigma + lambda] +=
                final_coef * (matrix_trace(wm_sigma_i_rho, vc_xy[lambda]) +
                              matrix_trace(wm_rho_i_sigma, vc_yx[lambda]));
          }
        }
      }
    }
  });
}

inline void contract_four_combine(
    SlTable& t, SlTable& t_pi, const Complex& coef, const Geometry& geo,
    const Coordinate& xg_x, const Coordinate& xg_y,
    const SelectedField<Complex>& f_loop_i_rho_sigma_lambda,
    const SelectedField<ManyMagneticMoments>& smf, const FieldSelection& fsel)
// x and y are global coordinates
{
  TIMER("contract_four_combine");
  qassert(f_loop_i_rho_sigma_lambda.geo().multiplicity == 3 * 4 * 4 * 4);
  const Coordinate total_site = geo.total_site();
  SelectedField<Complex> fsum;
  fsum.init(fsel, 2);
  set_zero(fsum);
  qacc_for(idx, fsel.n_elems, {
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Vector<Complex> v_loop =
        f_loop_i_rho_sigma_lambda.get_elems_const(idx);
    const ManyMagneticMoments& mmm = smf.get_elem(idx);
    const ManyMagneticMoments pion_proj =
        pion_projection(xg_x, xg_y, xg, total_site, true);
    Vector<Complex> sums = fsum.get_elems(idx);
    Complex sum = 0;
    Complex pi_sum = 0;
    Complex pi_proj_sum = 0;
    for (int i = 0; i < 3; ++i) {
      for (int rho = 0; rho < 4; ++rho) {
        for (int sigma = 0; sigma < 4; ++sigma) {
          for (int lambda = 0; lambda < 4; ++lambda) {
            const Complex val =
                coef * v_loop[64 * i + 16 * rho + 4 * sigma + lambda];
            sum += val * get_m_comp(mmm, i, rho, sigma, lambda);
            pi_sum += get_m_comp(pion_proj, i, rho, sigma, lambda) *
                      get_m_comp(mmm, i, rho, sigma, lambda);
            pi_proj_sum += val * get_m_comp(pion_proj, i, rho, sigma, lambda);
          }
        }
      }
    }
    sums[0] += sum;                   // total contribution
    sums[1] += pi_sum * pi_proj_sum;  // total pion pole (and other pseudo
                                      // scalar exchange type) contribution
  });
  t.init(total_site);
  t_pi.init(total_site);
  qfor(idx, fsel.n_elems, {
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Vector<Complex> sums = fsum.get_elems_const(idx);
    const Complex& sum = sums[0];
    const Complex& sum_pi = sums[1];
    add_to_sl_table(t, sum, xg_x, xg_y, xg, total_site);
    add_to_sl_table(t_pi, sum_pi, xg_x, xg_y, xg, total_site);
  });
  acc_sl_table(t);
  acc_sl_table(t_pi);
  glb_sum_double_vec(get_data(t.table));
  glb_sum_double_vec(get_data(t_pi.table));
}

inline std::vector<std::string> get_clbl_inf_ref_tags(
    const std::string& job_tag)
{
  std::vector<std::string> tags;
  tags.push_back("ref-far");
  tags.push_back("ref-close");
  tags.push_back("ref-center");
  return tags;
}

// ------------------------------------------------------------------------

inline std::vector<std::string> contract_four_pair_labels(
    const std::vector<std::string>& tags)
{
  TIMER_VERBOSE("contract_four_pair_labels");
  std::vector<std::string> labels;
  for (int i = 0; i < (int)tags.size(); ++i) {
    std::string label = tags[i] + " proj-all";
    std::string label_pi = tags[i] + " proj-pi";
    labels.push_back(label);
    labels.push_back(label_pi);
  }
  return labels;
}

inline std::vector<SlTable> contract_four_pair(
    const Complex& coef, const FieldSelection& fsel,
    const SelectedField<RealD>& smf_d,
    const SelectedField<WilsonMatrix>& sprop_x,
    const SelectedField<WilsonMatrix>& sprop_y, const Coordinate& xg_x,
    const Coordinate& xg_y, const int type, const double weight_pair,
    const std::vector<std::string>& tags, const double muon_mass,
    const double z_v)
// default coef = 1.0
// type = 0 : light quark
// type = 1 : strange quark
// tags can include "ref-far", "ref-center", "ref-close"
{
  TIMER_VERBOSE("contract_four_pair");
  const SelectedField<ManyMagneticMoments>& smf =
      qcast_const<ManyMagneticMoments, RealD>(smf_d);
  const Geometry& geo = fsel.f_rank.geo();
  SelectedField<WilsonMatrix> sc_xy, sc_yx;
  set_local_current_from_props(sc_xy, sprop_y, sprop_x, fsel);
  set_local_current_from_props(sc_yx, sprop_x, sprop_y, fsel);
  CurrentMoments<WilsonMatrix> cm_xy, cm_yx;
  set_current_moments_from_current(cm_xy, sc_xy, fsel);
  set_current_moments_from_current(cm_yx, sc_yx, fsel);
  const double alpha_inv = 137.035999139;
  const double e_charge = std::sqrt(4 * qlat::PI / alpha_inv);
  const Complex coef0 = 1.0E10 * 2.0 * muon_mass * std::pow(e_charge, 6);
  const Complex coef1 =
      (type == 0 ? 16.0 + 1.0 : 1.0) / 81.0 * (-3.0) * std::pow(z_v, 4);
  const Complex coef_all = coef * coef0 * coef1 / 3.0 * weight_pair;
  std::vector<SlTable> ts;
  for (int i = 0; i < (int)tags.size(); ++i) {
    SelectedField<Complex> f_loop_i_rho_sigma_lambda;
    contract_four_loop(f_loop_i_rho_sigma_lambda, 1.0, xg_x, xg_y, sc_xy, sc_yx,
                       cm_xy, cm_yx, fsel, tags[i]);
    SlTable t, t_pi;
    contract_four_combine(t, t_pi, coef_all, geo, xg_x, xg_y,
                          f_loop_i_rho_sigma_lambda, smf, fsel);
    ts.push_back(t);
    ts.push_back(t_pi);
  }
  qcast_const<RealD, ManyMagneticMoments>(smf);
  return ts;
}

// ------------------------------------------------------------------------

inline std::vector<std::string> contract_two_plus_two_pair_labels()
{
  TIMER("contract_two_plus_two_pair_labels");
  std::vector<std::string> tags;
  tags.push_back("sub");
  tags.push_back("dsub");
  tags.push_back("sub pisl");
  tags.push_back("dsub pisl");
  std::vector<std::string> labels;
  for (int i = 0; i < (int)tags.size(); ++i) {
    std::string label = tags[i] + " proj-all";
    std::string label_pi = tags[i] + " proj-pi";
    labels.push_back(label);
    labels.push_back(label_pi);
  }
  return labels;
}

inline std::vector<SlTable> contract_two_plus_two_pair_no_glb_sum(
    long& n_points_in_r_sq_limit, long& n_points_computed, const Complex& coef,
    const Field<RealD>& rand_prob_sel_field,
    const Field<Complex>& hvp_x,
    const SelectedPoints<Complex>& edl_list_c, const Coordinate& xg_x,
    const PointsSelection& psel_edl, const long r_sq_limit,
    const double hvp_sel_threshold, const double weight_pair_0,
    const double muon_mass, const double z_v)
// hvp point source at x (rho)
// hvp point sink at y (sigma)
// hvp with external loop source at z (lambda)
// hvp with external loop sink summed over (i)
// psel_edl[k] = xg_z
// edl_list[k][i * 4 + lambda]
// -hvp_x.get_elem(xl_y, sigma * 4 + rho)
//
// glb_sum for SlTable not yet performed
{
  TIMER_VERBOSE("contract_two_plus_two_pair_no_glb_sum");
  qassert(rand_prob_sel_field.geo().multiplicity == 1);
  qassert(hvp_x.geo().multiplicity == 16);
  n_points_in_r_sq_limit = 0;
  n_points_computed = 0;
  const Geometry& geo = geo_reform(hvp_x.geo());
  const Coordinate total_site = geo.total_site();
  const double alpha_inv = 137.035999139;
  const double e_charge = std::sqrt(4 * qlat::PI / alpha_inv);
  const Complex coef0 = 1.0E10 * 2.0 * muon_mass * std::pow(e_charge, 6);
  const Complex coef1 = 25.0 / 81.0 * 3.0 * std::pow(z_v, 4);
  const Complex coef2 = coef * coef0 * coef1 / 3.0;
  const int sub_tag = 0;  // use subtracted muon line weighting function
  const long n_labels = 8;
  const long n_points = psel_edl.size();
  const SelectedPoints<array<ComplexD, 3 * 4>>& edl_list =
      qcast_const<array<ComplexD, 3 * 4>, ComplexD>(edl_list_c);
  qassert(n_points == edl_list.n_points);
  qassert(1 == edl_list.multiplicity);
  bool has_same_x_z = false;
  qfor(k, n_points, {
    const Coordinate& xg_z = psel_edl[k];
    if (xg_z == xg_x) {
      has_same_x_z = true;
      break;
    }
  });
  const long total_volume = geo.total_volume();
  double weight_pair = 0.0;
  if (has_same_x_z) {
    weight_pair = (total_volume - 1) / (n_points - 1);
  } else {
    weight_pair = (total_volume - 1) / n_points;
  }
  vector_acc<ManyMagneticMoments> mmm_0_list(n_points);
  qfor(k, n_points, {
    const Coordinate& xg_z = psel_edl[k];
    mmm_0_list[k] = get_muon_line_m_extra_lat(xg_x, xg_x, xg_z, total_site,
                                              muon_mass, sub_tag);
  });
  qassert(n_labels == (long)contract_two_plus_two_pair_labels().size());
  std::vector<SlTable> ts(n_labels);
  for (long i = 0; i < n_labels; ++i) {
    ts[i].init(total_site);
  }
  // original
  vector_acc<Complex> sums_sub(n_points, 0.0);
  vector_acc<Complex> sums_dsub(n_points, 0.0);
  // pion projection based on location same as connected diagram
  vector_acc<Complex> sums_sub_pi(n_points, 0.0);
  vector_acc<Complex> sums_dsub_pi(n_points, 0.0);
  // pion projection based on loop
  vector_acc<Complex> sums_sub_pi_pisl(n_points, 0.0);
  vector_acc<Complex> sums_dsub_pi_pisl(n_points, 0.0);
  qfor(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg_y = geo.coordinate_g_from_l(xl);
    if (sqr(smod(xg_y - xg_x, total_site)) > r_sq_limit) {
      continue;
    }
    const Vector<Complex> vhvp = hvp_x.get_elems_const(xl);
    qassert(vhvp.size() == 16);
    const double hvp_sel_ratio = std::sqrt(qnorm(vhvp)) / hvp_sel_threshold;
    const double prob = hvp_sel_ratio > 1.0 ? 1.0 : hvp_sel_ratio;
    qassert(prob > 0.0 and 1.0 >= prob);
    n_points_in_r_sq_limit += 1;
    if (prob < rand_prob_sel_field.get_elem(xl)) {
      continue;
    }
    n_points_computed += 1;
    // displayln_info(
    //     ssprintf("compute point with index=%ld prob=%.8lf", index, prob));
    const double weight1 = 1.0 / prob;
    set_zero(sums_sub);
    set_zero(sums_dsub);
    set_zero(sums_sub_pi);
    set_zero(sums_dsub_pi);
    set_zero(sums_sub_pi_pisl);
    set_zero(sums_dsub_pi_pisl);
    qthread_for(k, n_points, {
      const array<Complex, 3 * 4>& edl = edl_list.get_elem(k);
      const Coordinate& xg_z = psel_edl[k];
      if (sqr(smod(xg_z - xg_x, total_site)) > r_sq_limit or
          sqr(smod(xg_z - xg_y, total_site)) > r_sq_limit) {
        continue;
      }
      const double weight2 = xg_z == xg_x ? weight_pair_0 : weight_pair;
      const Complex coef_all = coef2 * weight1 * weight2;
      const ManyMagneticMoments mmm = get_muon_line_m_extra_lat(
          xg_x, xg_y, xg_z, total_site, muon_mass, sub_tag);
      const ManyMagneticMoments mmm_dsub = mmm - mmm_0_list[k];
      const ManyMagneticMoments pion_proj =
          pion_projection(xg_x, xg_y, xg_z, total_site, true);
      const ManyMagneticMoments pion_proj_pisl = permute_rho_sigma_nu(
          pion_projection(xg_z, xg_x, xg_y, total_site, false), 1, 2, 0);
      Complex sub_sum = 0;
      Complex dsub_sum = 0;
      Complex sub_pi_sum = 0;
      Complex dsub_pi_sum = 0;
      Complex pi_proj_sum = 0;
      Complex sub_pi_pisl_sum = 0;
      Complex dsub_pi_pisl_sum = 0;
      Complex pi_proj_pisl_sum = 0;
      for (int i = 0; i < 3; ++i) {
        for (int rho = 0; rho < 4; ++rho) {
          for (int sigma = 0; sigma < 4; ++sigma) {
            for (int lambda = 0; lambda < 4; ++lambda) {
              const Complex val =
                  coef_all * edl[i * 4 + lambda] * (-vhvp[sigma * 4 + rho]);
              sub_sum += val * get_m_comp(mmm, i, rho, sigma, lambda);
              dsub_sum += val * get_m_comp(mmm_dsub, i, rho, sigma, lambda);
              sub_pi_sum += get_m_comp(pion_proj, i, rho, sigma, lambda) *
                            get_m_comp(mmm, i, rho, sigma, lambda);
              dsub_pi_sum += get_m_comp(pion_proj, i, rho, sigma, lambda) *
                             get_m_comp(mmm_dsub, i, rho, sigma, lambda);
              pi_proj_sum += val * get_m_comp(pion_proj, i, rho, sigma, lambda);
              sub_pi_pisl_sum +=
                  get_m_comp(pion_proj_pisl, i, rho, sigma, lambda) *
                  get_m_comp(mmm, i, rho, sigma, lambda);
              dsub_pi_pisl_sum +=
                  get_m_comp(pion_proj_pisl, i, rho, sigma, lambda) *
                  get_m_comp(mmm_dsub, i, rho, sigma, lambda);
              pi_proj_pisl_sum +=
                  val * get_m_comp(pion_proj_pisl, i, rho, sigma, lambda);
            }
          }
        }
      }
      // original
      sums_sub[k] += sub_sum;
      sums_dsub[k] += dsub_sum;
      // pion polarization projection (based on vertex location, same as
      // connected diagram)
      sums_sub_pi[k] += sub_pi_sum * pi_proj_sum;
      sums_dsub_pi[k] += dsub_pi_sum * pi_proj_sum;
      // pion polarization projection based on loop topology (only possible for
      // disconnected diagram)
      sums_sub_pi_pisl[k] += sub_pi_pisl_sum * pi_proj_pisl_sum;
      sums_dsub_pi_pisl[k] += dsub_pi_pisl_sum * pi_proj_pisl_sum;
    });
    qcast_const<ComplexD, array<ComplexD, 3 * 4>>(edl_list);
    qfor(k, n_points, {
      const Coordinate& xg_z = psel_edl[k];
      add_to_sl_table(ts[0], sums_sub[k], xg_x, xg_y, xg_z, total_site);
      add_to_sl_table(ts[1], sums_sub_pi[k], xg_x, xg_y, xg_z, total_site);
      add_to_sl_table(ts[2], sums_dsub[k], xg_x, xg_y, xg_z, total_site);
      add_to_sl_table(ts[3], sums_dsub_pi[k], xg_x, xg_y, xg_z, total_site);
      add_to_sl_table(ts[4], sums_sub[k], sqr(smod(xg_y - xg_x, total_site)),
                      sqr(smod(xg_z - xg_x, total_site)));
      add_to_sl_table(ts[5], sums_sub_pi_pisl[k],
                      sqr(smod(xg_y - xg_x, total_site)),
                      sqr(smod(xg_z - xg_x, total_site)));
      add_to_sl_table(ts[6], sums_dsub[k], sqr(smod(xg_y - xg_x, total_site)),
                      sqr(smod(xg_z - xg_x, total_site)));
      add_to_sl_table(ts[7], sums_dsub_pi_pisl[k],
                      sqr(smod(xg_y - xg_x, total_site)),
                      sqr(smod(xg_z - xg_x, total_site)));
    });
  });
  // no glb sum performed
  for (long i = 0; i < n_labels; ++i) {
    acc_sl_table(ts[i]);
  }
  displayln_info(
      fname +
      ssprintf(
          ": n_points_in_r_sq_limit=%ld n_points_computed=%ld ratio=%.10lf",
          n_points_in_r_sq_limit, n_points_computed,
          (double)n_points_computed / (double)n_points_in_r_sq_limit));
  return ts;
}

}  // namespace qlat
