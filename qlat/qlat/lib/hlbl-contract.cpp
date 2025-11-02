#include <qlat/hlbl-contract.h>

namespace qlat
{  //

#define qacc_for_debug qthread_for

void set_m_z_field_tag(SelectedPoints<RealD>& smf_d,
                       const PointsSelection& psel_d, const Geometry& geo,
                       const Coordinate& xg_x, const Coordinate& xg_y,
                       const double a, const Int tag)
// interface
// tag = 0 sub
// tag = 1 nosub
{
  TIMER_VERBOSE("set_m_z_field_tag(smf_d,psel_d,xg_x,xg_y,a,tag)");
  const Int multiplicity = sizeof(ManyMagneticMoments) / sizeof(RealD);
  Qassert(multiplicity * (int)sizeof(RealD) ==
          (int)sizeof(ManyMagneticMoments));
  smf_d.init(psel_d, multiplicity);
  SelectedPoints<ManyMagneticMoments> smf;
  smf.set_view_cast(smf_d);
  const Coordinate total_site = geo.total_site();
  qthread_for(idx, psel_d.size(), {
    const Coordinate xg_z = psel_d[idx];  // z location
    ManyMagneticMoments& mmm = smf.get_elem(idx);
    mmm = get_muon_line_m_extra_lat(xg_x, xg_y, xg_z, total_site, a, tag);
  });
}

void set_local_current_from_props(SelectedPoints<WilsonMatrix>& scf,
                                  const SelectedPoints<WilsonMatrix>& sprop1,
                                  const SelectedPoints<WilsonMatrix>& sprop2,
                                  const PointsSelection& psel_d,
                                  const Geometry& geo)
// -<- gamma5 sprop2^+ gamma5 -<- gamma_mu -<- sprop1 -<-
{
  (void)geo;
  TIMER_VERBOSE("set_local_current_from_props");
  const array<SpinMatrix, 4>& gammas = SpinMatrixConstants::get_cps_gammas();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  scf.init(psel_d, 4);
  set_zero(scf);
  qacc_for(idx, psel_d.size(), {
    const WilsonMatrix& m1 = sprop1.get_elem(idx);
    const WilsonMatrix& m2 = sprop2.get_elem(idx);
    Vector<WilsonMatrix> v = scf.get_elems(idx);
    const WilsonMatrix m2rev =
        gamma5 * (WilsonMatrix)matrix_adjoint(m2) * gamma5;
    for (Int m = 0; m < 4; ++m) {
      v[m] = m2rev * gammas[m] * m1;
    }
  });
}

RealD set_psel_d_prob_xy(SelectedPoints<RealD>& psel_d_prob_xy,
                         const PointsSelection& psel,
                         SelectedPoints<RealD>& psel_prob,
                         const PointsSelection& psel_d,
                         const SelectedPoints<RealD>& psel_d_prob,
                         const Long idx_xg_x, const Long idx_xg_y)
// return prob_pair;
{
  TIMER("set_psel_d_prob_xy(psel_d_prob_xy,...)");
  psel_d_prob_xy = psel_d_prob;
  Qassert(0 <= idx_xg_x and idx_xg_x < (Long)psel.size());
  Qassert(0 <= idx_xg_y and idx_xg_y < (Long)psel.size());
  const RealD prob_xg_x = psel_prob.get_elem(idx_xg_x);
  const RealD prob_xg_y = psel_prob.get_elem(idx_xg_y);
  const Coordinate& xg_x = psel[idx_xg_x];
  const Coordinate& xg_y = psel[idx_xg_y];
  RealD prob_pair = 0.0;
  if (xg_x == xg_y) {
    Qassert(idx_xg_x == idx_xg_y);
    Qassert(prob_xg_x == prob_xg_y);
    prob_pair = prob_xg_x;
    qthread_for(idx, psel_d.size(), {
      const Coordinate xg_z = psel_d[idx];
      if (xg_z == xg_x) {
        Qassert(xg_z == xg_y);
        const RealD prob = psel_d_prob_xy.get_elem(idx);
        psel_d_prob_xy.get_elem(idx) = std::min(1.0, prob / prob_xg_x);
      }
    });
  } else {
    prob_pair = prob_xg_x * prob_xg_y;
    qthread_for(idx, psel_d.size(), {
      const Coordinate xg_z = psel_d[idx];
      if (xg_z == xg_x) {
        const RealD prob = psel_d_prob_xy.get_elem(idx);
        psel_d_prob_xy.get_elem(idx) = std::min(1.0, prob / prob_xg_x);
      } else if (xg_z == xg_y) {
        const RealD prob = psel_d_prob_xy.get_elem(idx);
        psel_d_prob_xy.get_elem(idx) = std::min(1.0, prob / prob_xg_y);
      }
    });
  }
  return prob_pair;
}

void set_current_moments_from_current_par(
    CurrentMoments<WilsonMatrix>& cm,
    const SelectedPoints<WilsonMatrix>& current, const PointsSelection& psel_d,
    const SelectedPoints<RealD>& psel_d_prob_xy, const Geometry& geo)
{
  TIMER("set_current_moments_from_current_par");
  Qassert(current.multiplicity == 4);
  const Coordinate total_site = geo.total_site();
  const Int lsize =
      std::max(total_site[0], std::max(total_site[1], total_site[2]));
  cm.init(lsize);
  qthread_for(xg_op_j, lsize, {
    qfor(idx, psel_d.size(), {
      const Coordinate& xg_op = psel_d[idx];  // z location
      for (Int j = 0; j < 3; ++j) {
        if (xg_op[j] == xg_op_j) {
          const Vector<WilsonMatrix> v = current.get_elems_const(idx);
          const RealD prob = psel_d_prob_xy.get_elem(idx);
          const ComplexD weight = 1.0 / prob;
          for (Int k = 0; k < 3; ++k) {
            cm.d[xg_op_j][3 * j + k] += weight * v[k];
          }
        }
      }
    });
  });
}

void set_current_moments_from_current_nopar(
    CurrentMoments<WilsonMatrix>& cm,
    const SelectedPoints<WilsonMatrix>& current, const PointsSelection& psel_d,
    const SelectedPoints<RealD>& psel_d_prob_xy, const Geometry& geo)
{
  TIMER("set_current_moments_from_current_nopar");
  Qassert(current.multiplicity == 4);
  const Coordinate total_site = geo.total_site();
  const Int lsize =
      std::max(total_site[0], std::max(total_site[1], total_site[2]));
  cm.init(lsize);
  qfor(idx, psel_d.size(), {
    const Coordinate& xg_op = psel_d[idx];  // z location
    const Vector<WilsonMatrix> v = current.get_elems_const(idx);
    const RealD prob = psel_d_prob_xy.get_elem(idx);
    const ComplexD weight = 1.0 / prob;
    for (Int j = 0; j < 3; ++j) {
      const Int xg_op_j = xg_op[j];
      if (xg_op[j] == xg_op_j) {
        for (Int k = 0; k < 3; ++k) {
          cm.d[xg_op_j][3 * j + k] += weight * v[k];
        }
      }
    }
  });
}

void set_current_moments_from_current(
    CurrentMoments<WilsonMatrix>& cm,
    const SelectedPoints<WilsonMatrix>& current, const PointsSelection& psel_d,
    const SelectedPoints<RealD>& psel_d_prob_xy, const Geometry& geo)
{
  // set_current_moments_from_current_nopar(cm, current, psel_d, psel_d_prob_xy,
  //                                        geo);
  set_current_moments_from_current_par(cm, current, psel_d, psel_d_prob_xy,
                                       geo);
}

void glb_sum_current_moments(CurrentMoments<WilsonMatrix>& cm)
{
  TIMER("glb_sum_current_moments");
  glb_sum(cm.d);
}

void contract_four_loop(SelectedPoints<Complex>& f_loop_i_rho_sigma_lambda,
                        const Complex& coef, const Coordinate& xg_x,
                        const Coordinate& xg_y,
                        const SelectedPoints<WilsonMatrix>& c_xy,
                        const SelectedPoints<WilsonMatrix>& c_yx,
                        const CurrentMoments<WilsonMatrix>& cm_xy,
                        const CurrentMoments<WilsonMatrix>& cm_yx,
                        const PointsSelection& psel_d,
                        const SelectedPoints<RealD>& psel_d_prob_xy,
                        const Geometry& geo, const Long r_sq_limit,
                        const std::string& label)
{
  TIMER_VERBOSE("contract_four_loop");
  const box<SpinMatrixConstantsT<>>& smc = get_spin_matrix_constants();
  const Coordinate total_site = geo.total_site();
  f_loop_i_rho_sigma_lambda.init(psel_d, 3 * 4 * 4 * 4);
  set_zero(f_loop_i_rho_sigma_lambda);
  const ChooseReferenceLabel cr_label = choose_reference_label(label);
  SelectedPoints<WilsonMatrix> f_sm_yx_g, f_sm_xy_g, f_vc_yx_g, f_vc_xy_g;
  f_sm_yx_g.init(psel_d, 3 * 4);
  f_sm_xy_g.init(psel_d, 3 * 4);
  f_vc_yx_g.init(psel_d, 4 * 4);
  f_vc_xy_g.init(psel_d, 4 * 4);
  qacc_for_debug(idx, psel_d.size(), {
    const array<SpinMatrix, 4>& gammas = smc().cps_gammas;
    const RealD prob = psel_d_prob_xy.get_elem(idx);
    const RealD weight = 1.0 / prob;
    const ComplexD final_coef = coef * weight;
    const Coordinate xg_z = psel_d[idx];
    if (sqr(smod(xg_z - xg_x, total_site)) <= r_sq_limit and
        sqr(smod(xg_z - xg_y, total_site)) <= r_sq_limit) {
      const CoordinateD xgref =
          choose_reference(xg_x, xg_y, xg_z, total_site, cr_label);
      const array<WilsonMatrix, 3> sm_yx = simple_moment_with_contact_subtract(
          cm_yx, xgref, total_site, c_yx, psel_d, psel_d_prob_xy, idx);
      const array<WilsonMatrix, 3> sm_xy = simple_moment_with_contact_subtract(
          cm_xy, xgref, total_site, c_xy, psel_d, psel_d_prob_xy, idx);
      const Vector<WilsonMatrix> vc_xy = c_xy.get_elems_const(idx);
      const Vector<WilsonMatrix> vc_yx = c_yx.get_elems_const(idx);
      Vector<Complex> v_loop = f_loop_i_rho_sigma_lambda.get_elems(idx);
      // matrix_trace(sm_yx[i] * gammas[rho], vc_xy[lambda] * gammas[sigma])
      // matrix_trace(sm_xy[i] * gammas[sigma], vc_yx[lambda] * gammas[rho])
      Vector<WilsonMatrix> sm_yx_g = f_sm_yx_g.get_elems(idx);
      Vector<WilsonMatrix> sm_xy_g = f_sm_xy_g.get_elems(idx);
      // array<WilsonMatrix, 3 * 4> sm_yx_g;
      // array<WilsonMatrix, 3 * 4> sm_xy_g;
      for (Int i = 0; i < 3; ++i) {
        for (Int rho = 0; rho < 4; ++rho) {
          sm_yx_g[i * 4 + rho] = sm_yx[i] * gammas[rho];
          sm_xy_g[i * 4 + rho] = sm_xy[i] * gammas[rho];
        }
      }
      Vector<WilsonMatrix> vc_yx_g = f_vc_yx_g.get_elems(idx);
      Vector<WilsonMatrix> vc_xy_g = f_vc_xy_g.get_elems(idx);
      // array<WilsonMatrix, 4 * 4> vc_yx_g;
      // array<WilsonMatrix, 4 * 4> vc_xy_g;
      for (Int lambda = 0; lambda < 4; ++lambda) {
        for (Int sigma = 0; sigma < 4; ++sigma) {
          vc_xy_g[lambda * 4 + sigma] = vc_xy[lambda] * gammas[sigma];
          vc_yx_g[lambda * 4 + sigma] = vc_yx[lambda] * gammas[sigma];
        }
      }
      for (Int i = 0; i < 3; ++i) {
        for (Int rho = 0; rho < 4; ++rho) {
          for (Int sigma = 0; sigma < 4; ++sigma) {
            for (Int lambda = 0; lambda < 4; ++lambda) {
              v_loop[64 * i + 16 * rho + 4 * sigma + lambda] +=
                  final_coef * (matrix_trace(sm_yx_g[i * 4 + rho],
                                             vc_xy_g[lambda * 4 + sigma]) +
                                matrix_trace(sm_xy_g[i * 4 + sigma],
                                             vc_yx_g[lambda * 4 + rho]));
            }
          }
        }
      }
    }
  });
}

void contract_four_combine(
    SlTable& t, SlTable& t_pi, const Complex& coef, const Geometry& geo,
    const Coordinate& xg_x, const Coordinate& xg_y,
    const SelectedPoints<Complex>& f_loop_i_rho_sigma_lambda,
    const SelectedPoints<ManyMagneticMoments>& smf,
    const PointsSelection& psel_d, const Long r_sq_limit)
// x and y are global coordinates
{
  TIMER("contract_four_combine");
  Qassert(f_loop_i_rho_sigma_lambda.multiplicity == 3 * 4 * 4 * 4);
  const Coordinate total_site = geo.total_site();
  SelectedPoints<Complex> fsum;
  fsum.init(psel_d, 2);
  set_zero(fsum);
  qacc_for(idx, psel_d.size(), {
    const Coordinate xg_z = psel_d[idx];
    if (sqr(smod(xg_z - xg_x, total_site)) <= r_sq_limit and
        sqr(smod(xg_z - xg_y, total_site)) <= r_sq_limit) {
      const Vector<Complex> v_loop =
          f_loop_i_rho_sigma_lambda.get_elems_const(idx);
      const ManyMagneticMoments& mmm = smf.get_elem(idx);
      const ManyMagneticMoments pion_proj =
          pion_projection(xg_x, xg_y, xg_z, total_site, true);
      Vector<Complex> sums = fsum.get_elems(idx);
      Complex sum = 0;
      Complex pi_sum = 0;
      Complex pi_proj_sum = 0;
      for (Int i = 0; i < 3; ++i) {
        for (Int rho = 0; rho < 4; ++rho) {
          for (Int sigma = 0; sigma < 4; ++sigma) {
            for (Int lambda = 0; lambda < 4; ++lambda) {
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
    }
  });
  t.init(total_site);
  t_pi.init(total_site);
  qfor(idx, psel_d.size(), {
    const Coordinate xg_z = psel_d[idx];
    const Vector<Complex> sums = fsum.get_elems_const(idx);
    const Complex& sum = sums[0];
    const Complex& sum_pi = sums[1];
    add_to_sl_table(t, sum, xg_x, xg_y, xg_z, total_site);
    add_to_sl_table(t_pi, sum_pi, xg_x, xg_y, xg_z, total_site);
  });
  // no glb sum performed
  acc_sl_table(t);
  acc_sl_table(t_pi);
}

std::vector<SlTable> contract_four_pair_no_glb_sum(
    const ComplexD& coef, const PointsSelection& psel,
    const PointsSelection& psel_d, const SelectedPoints<RealD>& psel_d_prob_xy,
    const Geometry& geo, const Long idx_xg_x, const Long idx_xg_y,
    const SelectedPoints<RealD>& smf_d,
    const SelectedPoints<WilsonMatrix>& sc_xy,
    const SelectedPoints<WilsonMatrix>& sc_yx,
    const CurrentMoments<WilsonMatrix>& cm_xy,
    const CurrentMoments<WilsonMatrix>& cm_yx, const Int inv_type,
    const std::vector<std::string>& tags, const Long r_sq_limit,
    const RealD muon_mass, const RealD z_v)
// default coef = 1.0
// inv_type = 0 : light quark
// inv_type = 1 : strange quark
// tags can include "ref-far", "ref-center", "ref-close"
{
  TIMER_VERBOSE("contract_four_pair_no_glb_sum");
  Qassert(0 <= idx_xg_x and idx_xg_x < (Long)psel.size());
  Qassert(0 <= idx_xg_y and idx_xg_y < (Long)psel.size());
  const Coordinate& xg_x = psel[idx_xg_x];
  const Coordinate& xg_y = psel[idx_xg_y];
  const Long total_volume = geo.total_volume();
  const Coordinate total_site = geo.total_site();
  if (sqr(smod(xg_y - xg_x, total_site)) > r_sq_limit) {
    qerr(fname + ssprintf(": xg_x=%s xg_y=%s total_site=%s r_sq_limit=%ld.",
                          show(xg_x).c_str(), show(xg_y).c_str(),
                          show(total_site).c_str(), (long)r_sq_limit));
  }
  SelectedPoints<ManyMagneticMoments> smf;
  smf.set_view_cast(smf_d);
  Qassert(psel_d.size() == psel_d_prob_xy.n_points);
  Qassert(psel_d.size() == sc_xy.n_points);
  Qassert(psel_d.size() == sc_yx.n_points);
  const RealD alpha_inv = 137.035999139;
  const RealD e_charge = std::sqrt(4 * qlat::PI / alpha_inv);
  const Complex coef0 = 1.0E10 * 2.0 * muon_mass * std::pow(e_charge, 6);
  const Complex coef1 =
      (inv_type == 0 ? 16.0 + 1.0 : 1.0) / 81.0 * (-3.0) * std::pow(z_v, 4);
  const Complex coef_all = coef * coef0 * coef1 / 3.0 / (RealD)total_volume;
  std::vector<SlTable> ts;
  for (Int i = 0; i < (int)tags.size(); ++i) {
    SelectedPoints<Complex> f_loop_i_rho_sigma_lambda;
    contract_four_loop(f_loop_i_rho_sigma_lambda, 1.0, xg_x, xg_y, sc_xy, sc_yx,
                       cm_xy, cm_yx, psel_d, psel_d_prob_xy, geo, r_sq_limit,
                       tags[i]);
    SlTable t, t_pi;
    // no glb sum performed
    contract_four_combine(t, t_pi, coef_all, geo, xg_x, xg_y,
                          f_loop_i_rho_sigma_lambda, smf, psel_d, r_sq_limit);
    ts.push_back(t);
    ts.push_back(t_pi);
  }
  return ts;
}

std::vector<SlTable> contract_two_plus_two_pair_no_glb_sum(
    Long& n_points_selected, Long& n_points_computed, const ComplexD& coef,
    const Geometry& geo, const PointsSelection& psel,
    const SelectedPoints<RealD>& psel_prob, const PointsSelection& psel_lps,
    const SelectedPoints<RealD>& psel_lps_prob, const Long idx_xg_x,
    const SelectedPoints<ComplexD>& lps_hvp_x,
    const SelectedPoints<ComplexD>& edl_list_c, const Long r_sq_limit,
    const RealD muon_mass, const RealD z_v)
// hvp point source at x (rho)
// hvp point sink at y (sigma)
// hvp with external loop source at z (lambda)
// hvp with external loop sink summed over (i)
// psel[k] = xg_z
// edl_list_c[k][i * 4 + lambda]
// -lps_hvp_x.get_elem(idx_xl_y, sigma * 4 + rho)
//
// glb_sum for SlTable not yet performed
//
// `rand_prob_sel_field` be the same rand field as the rand field for sampling
// points.
{
  TIMER_VERBOSE("contract_two_plus_two_pair_no_glb_sum");
  Qassert(0 <= idx_xg_x and idx_xg_x < (Long)psel.size());
  Qassert(psel_prob.multiplicity == 1);
  Qassert(psel_prob.n_points == (Long)psel.size());
  Qassert(psel_lps_prob.multiplicity == 1);
  Qassert(psel_lps_prob.n_points == (Long)psel_lps.size());
  Qassert(lps_hvp_x.multiplicity == 16);
  Qassert(lps_hvp_x.n_points == (Long)psel_lps.size());
  Qassert(edl_list_c.multiplicity == 3 * 4);
  Qassert(edl_list_c.n_points == (Long)psel.size());
  const Coordinate& xg_x = psel[idx_xg_x];
  const RealD prob_xg_x = psel_prob.get_elem(idx_xg_x);
  n_points_selected = 0;
  n_points_computed = 0;
  const Long total_volume = geo.total_volume();
  const Coordinate total_site = geo.total_site();
  const double alpha_inv = 137.035999139;
  const double e_charge = std::sqrt(4 * qlat::PI / alpha_inv);
  const Complex coef0 = 1.0E10 * 2.0 * muon_mass * std::pow(e_charge, 6);
  const Complex coef1 = 25.0 / 81.0 * 3.0 * std::pow(z_v, 4);
  const Complex coef2 = coef * coef0 * coef1 / 3.0 / (RealD)total_volume;
  const Int sub_tag = 0;  // use subtracted muon line weighting function
  const long n_labels = 8;
  const long n_points = psel.size();
  bool has_same_x_z = false;
  qfor(k, n_points, {
    const Coordinate& xg_z = psel[k];
    if (xg_z == xg_x) {
      has_same_x_z = true;
      break;
    }
  });
  Qassert(has_same_x_z);
  vector<ManyMagneticMoments> mmm_0_list(n_points);
  qfor(k, n_points, {
    const Coordinate& xg_z = psel[k];
    mmm_0_list[k] = get_muon_line_m_extra_lat(xg_x, xg_x, xg_z, total_site,
                                              muon_mass, sub_tag);
  });
  Qassert(n_labels == (long)contract_two_plus_two_pair_labels().size());
  std::vector<SlTable> ts(n_labels);
  for (long i = 0; i < n_labels; ++i) {
    ts[i].init(total_site);
  }
  // original
  vector<Complex> sums_sub(n_points);
  vector<Complex> sums_dsub(n_points);
  // pion projection based on location same as connected diagram
  vector<Complex> sums_sub_pi(n_points);
  vector<Complex> sums_dsub_pi(n_points);
  // pion projection based on loop
  vector<Complex> sums_sub_pi_pisl(n_points);
  vector<Complex> sums_dsub_pi_pisl(n_points);
  //
  set_zero(sums_sub);
  set_zero(sums_dsub);
  set_zero(sums_sub_pi);
  set_zero(sums_dsub_pi);
  set_zero(sums_sub_pi_pisl);
  set_zero(sums_dsub_pi_pisl);
  //
  qfor(idx, (Long)psel_lps.size(), {
    const Coordinate xg_y = psel_lps[idx];
    n_points_selected += 1;
    if (sqr(smod(xg_y - xg_x, total_site)) > r_sq_limit) {
      continue;
    }
    n_points_computed += 1;
    const Vector<Complex> vhvp = lps_hvp_x.get_elems_const(idx);
    qassert(vhvp.size() == 16);
    const RealD prob_xg_y = psel_lps_prob.get_elem(idx);
    qassert(prob_xg_y > 0.0 and 1.0 >= prob_xg_y);
    set_zero(sums_sub);
    set_zero(sums_dsub);
    set_zero(sums_sub_pi);
    set_zero(sums_dsub_pi);
    set_zero(sums_sub_pi_pisl);
    set_zero(sums_dsub_pi_pisl);
    qthread_for(k, n_points, {
      const Coordinate& xg_z = psel[k];
      const RealD prob_xg_z = psel_prob.get_elem(k);
      const Vector<Complex> edl = edl_list_c.get_elems_const(k);
      qassert(edl.size() == 3 * 4);
      if (sqr(smod(xg_z - xg_x, total_site)) > r_sq_limit or
          sqr(smod(xg_z - xg_y, total_site)) > r_sq_limit) {
        continue;
      }
      RealD prob = 1.0;
      if (xg_x == xg_y) {
        prob = std::min(prob_xg_x, prob_xg_y);
        if (xg_x == xg_z) {
          prob = std::min(prob, prob_xg_z);
        } else {
          prob = prob * prob_xg_z;
        }
      } else {
        if (xg_x == xg_z) {
          prob = std::min(prob_xg_x, prob_xg_z) * prob_xg_y;
        } else if (xg_y == xg_z) {
          prob = prob_xg_x * std::min(prob_xg_y, prob_xg_z);
        } else {
          prob = prob_xg_x * prob_xg_y * prob_xg_z;
        }
      }
      const Complex coef_all = coef2 / prob;
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
      for (Int i = 0; i < 3; ++i) {
        for (Int rho = 0; rho < 4; ++rho) {
          for (Int sigma = 0; sigma < 4; ++sigma) {
            for (Int lambda = 0; lambda < 4; ++lambda) {
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
    qfor(k, n_points, {
      const Coordinate& xg_z = psel[k];
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
      ssprintf(": xg_x=%s ; n_points_selected=%ld ; n_points_computed=%ld "
               "; ratio=%.10lf .",
               show(xg_x).c_str(), n_points_selected, n_points_computed,
               (double)n_points_computed / (double)n_points_selected));
  return ts;
}

}  // namespace qlat
