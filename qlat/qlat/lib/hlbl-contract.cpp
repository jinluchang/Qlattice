#include <qlat/hlbl-contract.h>

namespace qlat
{  //

void set_m_z_field_tag(SelectedField<RealD>& smf_d, const FieldSelection& fsel,
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

void set_local_current_from_props(SelectedField<WilsonMatrix>& scf,
                                  const SelectedField<WilsonMatrix>& sprop1,
                                  const SelectedField<WilsonMatrix>& sprop2,
                                  const FieldSelection& fsel)
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

void contract_four_loop(SelectedField<Complex>& f_loop_i_rho_sigma_lambda,
                        const Complex& coef, const Coordinate& xg_x,
                        const Coordinate& xg_y,
                        const SelectedField<WilsonMatrix>& c_xy,
                        const SelectedField<WilsonMatrix>& c_yx,
                        const CurrentMoments<WilsonMatrix>& cm_xy,
                        const CurrentMoments<WilsonMatrix>& cm_yx,
                        const FieldSelection& fsel,
                        const SelectedField<RealD>& fsel_prob_xy,
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
  const ChooseReferenceLabel cr_label = choose_reference_label(label);
  qacc_for(idx, fsel.n_elems, {
    const array<SpinMatrix, 4>& gammas = smc().cps_gammas;
    const RealD prob = fsel_prob_xy.get_elem(idx);
    const RealD weight = 1.0 / prob;
    const ComplexD final_coef = coef * weight;
    const long index = fsel.indices[idx];
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const CoordinateD xgref =
        choose_reference(xg_x, xg_y, xg, total_site, cr_label);
    const array<WilsonMatrix, 3> sm_yx = simple_moment_with_contact_subtract(
        cm_yx, xgref, total_site, c_yx, fsel, fsel_prob_xy, xg);
    const array<WilsonMatrix, 3> sm_xy = simple_moment_with_contact_subtract(
        cm_xy, xgref, total_site, c_xy, fsel, fsel_prob_xy, xg);
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

void contract_four_combine(
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

std::vector<SlTable> contract_four_pair(
    const ComplexD& coef, const PointsSelection& psel,
    SelectedPoints<RealD>& psel_prob, const FieldSelection& fsel,
    const SelectedField<RealD>& fsel_prob, const Long idx_xg_x,
    const Long idx_xg_y, const SelectedField<RealD>& smf_d,
    const SelectedField<WilsonMatrix>& sprop_x,
    const SelectedField<WilsonMatrix>& sprop_y, const Int inv_type,
    const std::vector<std::string>& tags, const Long r_sq_limit,
    const RealD muon_mass, const RealD z_v)
// default coef = 1.0
// inv_type = 0 : light quark
// inv_type = 1 : strange quark
// tags can include "ref-far", "ref-center", "ref-close"
{
  // TODO: implement r_sq_limit
  TIMER_VERBOSE("contract_four_pair");
  qassert(0 <= idx_xg_x and idx_xg_x < (Long)psel.size());
  qassert(0 <= idx_xg_y and idx_xg_y < (Long)psel.size());
  const Coordinate& xg_x = psel[idx_xg_x];
  const Coordinate& xg_y = psel[idx_xg_y];
  const RealD prob_xg_x = psel_prob.get_elem(idx_xg_x);
  const RealD prob_xg_y = psel_prob.get_elem(idx_xg_y);
  const RealD prob_pair =
      xg_x == xg_y ? std::min(prob_xg_x, prob_xg_y) : prob_xg_x * prob_xg_y;
  const RealD weight_pair = 1.0 / prob_pair;
  const SelectedField<ManyMagneticMoments>& smf =
      qcast_const<ManyMagneticMoments, RealD>(smf_d);
  const Geometry& geo = fsel.f_rank.geo();
  qassert(is_matching_geo_mult(geo, fsel_prob.geo()));
  qassert(is_matching_geo_mult(geo, smf.geo()));
  qassert(is_matching_geo_mult(geo, sprop_x.geo()));
  qassert(is_matching_geo_mult(geo, sprop_y.geo()));
  qassert(fsel.n_elems == fsel_prob.n_elems);
  qassert(fsel.n_elems == sprop_x.n_elems);
  qassert(fsel.n_elems == sprop_y.n_elems);
  SelectedField<RealD> fsel_prob_xy;
  fsel_prob_xy = fsel_prob;
  const Long idx_x = idx_from_xg(xg_x, fsel);
  const Long idx_y = idx_from_xg(xg_y, fsel);
  if (idx_x == idx_y) {
    if (idx_x >= 0) {
      const RealD prob = fsel_prob_xy.get_elem(idx_x);
      fsel_prob_xy.get_elem(idx_x) =
          std::min(1.0, prob / std::min(prob_xg_x, prob_xg_y));
    }
  } else {
    if (idx_x >= 0) {
      const RealD prob = fsel_prob_xy.get_elem(idx_x);
      fsel_prob_xy.get_elem(idx_x) = std::min(1.0, prob / prob_xg_x);
    }
    if (idx_y >= 0) {
      const RealD prob = fsel_prob_xy.get_elem(idx_y);
      fsel_prob_xy.get_elem(idx_y) = std::min(1.0, prob / prob_xg_y);
    }
  }
  SelectedField<WilsonMatrix> sc_xy, sc_yx;
  set_local_current_from_props(sc_xy, sprop_y, sprop_x, fsel);
  set_local_current_from_props(sc_yx, sprop_x, sprop_y, fsel);
  CurrentMoments<WilsonMatrix> cm_xy, cm_yx;
  set_current_moments_from_current(cm_xy, sc_xy, fsel, fsel_prob_xy);
  set_current_moments_from_current(cm_yx, sc_yx, fsel, fsel_prob_xy);
  const RealD alpha_inv = 137.035999139;
  const RealD e_charge = std::sqrt(4 * qlat::PI / alpha_inv);
  const Complex coef0 = 1.0E10 * 2.0 * muon_mass * std::pow(e_charge, 6);
  const Complex coef1 =
      (inv_type == 0 ? 16.0 + 1.0 : 1.0) / 81.0 * (-3.0) * std::pow(z_v, 4);
  const Complex coef_all = coef * coef0 * coef1 / 3.0 * weight_pair;
  std::vector<SlTable> ts;
  for (int i = 0; i < (int)tags.size(); ++i) {
    SelectedField<Complex> f_loop_i_rho_sigma_lambda;
    contract_four_loop(f_loop_i_rho_sigma_lambda, 1.0, xg_x, xg_y, sc_xy, sc_yx,
                       cm_xy, cm_yx, fsel, fsel_prob_xy, tags[i]);
    SlTable t, t_pi;
    contract_four_combine(t, t_pi, coef_all, geo, xg_x, xg_y,
                          f_loop_i_rho_sigma_lambda, smf, fsel);
    ts.push_back(t);
    ts.push_back(t_pi);
  }
  qcast_const<RealD, ManyMagneticMoments>(smf);
  return ts;
}

std::vector<SlTable> contract_two_plus_two_pair_no_glb_sum(
    Long& n_points_in_r_sq_limit, Long& n_points_computed, const ComplexD& coef,
    const PointsSelection& psel, const SelectedPoints<RealD> psel_prob,
    const Field<RealD>& rand_prob_sel_field, const RealD hvp_sel_threshold,
    const Long idx_xg_x, const Field<ComplexD>& hvp_x,
    const SelectedPoints<ComplexD>& edl_list_c, const Long r_sq_limit,
    const RealD muon_mass, const RealD z_v)
// hvp point source at x (rho)
// hvp point sink at y (sigma)
// hvp with external loop source at z (lambda)
// hvp with external loop sink summed over (i)
// psel[k] = xg_z
// edl_list[k][i * 4 + lambda]
// -hvp_x.get_elem(xl_y, sigma * 4 + rho)
//
// glb_sum for SlTable not yet performed
//
// `rand_prob_sel_field` be the same rand field as the rand field for sampling
// points.
{
  TIMER_VERBOSE("contract_two_plus_two_pair_no_glb_sum");
  qassert(0 <= idx_xg_x and idx_xg_x < (Long)psel.size());
  const Coordinate& xg_x = psel[idx_xg_x];
  const RealD prob_xg_x = psel_prob.get_elem(idx_xg_x);
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
  const long n_points = psel.size();
  const SelectedPoints<array<ComplexD, 3 * 4>>& edl_list =
      qcast_const<array<ComplexD, 3 * 4>, ComplexD>(edl_list_c);
  qassert(n_points == edl_list.n_points);
  qassert(1 == edl_list.multiplicity);
  bool has_same_x_z = false;
  qfor(k, n_points, {
    const Coordinate& xg_z = psel[k];
    if (xg_z == xg_x) {
      has_same_x_z = true;
      break;
    }
  });
  qassert(has_same_x_z);
  vector_acc<ManyMagneticMoments> mmm_0_list(n_points);
  qfor(k, n_points, {
    const Coordinate& xg_z = psel[k];
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
    const RealD prob_xg_y = std::min(1.0, hvp_sel_ratio);
    qassert(prob_xg_y > 0.0 and 1.0 >= prob_xg_y);
    n_points_in_r_sq_limit += 1;
    if (prob_xg_y < rand_prob_sel_field.get_elem(xl)) {
      continue;
    }
    n_points_computed += 1;
    // displayln_info(
    //     ssprintf("compute point with index=%ld prob=%.8lf", index, prob));
    set_zero(sums_sub);
    set_zero(sums_dsub);
    set_zero(sums_sub_pi);
    set_zero(sums_dsub_pi);
    set_zero(sums_sub_pi_pisl);
    set_zero(sums_dsub_pi_pisl);
    qthread_for(k, n_points, {
      const Coordinate& xg_z = psel[k];
      const RealD prob_xg_z = psel_prob.get_elem(k);
      const array<Complex, 3 * 4>& edl = edl_list.get_elem(k);
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
      ssprintf(
          ": n_points_in_r_sq_limit=%ld n_points_computed=%ld ratio=%.10lf",
          n_points_in_r_sq_limit, n_points_computed,
          (double)n_points_computed / (double)n_points_in_r_sq_limit));
  return ts;
}

}  // namespace qlat
