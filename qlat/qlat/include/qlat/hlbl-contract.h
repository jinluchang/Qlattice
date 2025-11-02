#pragma once

#include "hlbl-sl-table.h"
#include "muon-line.h"

namespace qlat
{  //

void set_m_z_field_tag(SelectedPoints<RealD>& smf_d,
                       const PointsSelection& psel_d, const Geometry& geo,
                       const Coordinate& xg_x, const Coordinate& xg_y,
                       const RealD a, const Int tag);

// ------------------------------------------------------------------------

qacc ManyMagneticMoments simple_pion_projection(const CoordinateD& x,
                                                const CoordinateD& y,
                                                const CoordinateD& z)
{
  const CoordinateD mid_yz = (y + z) / 2;
  const CoordinateD x_mid_yz = x - mid_yz;
  const CoordinateD y_z = y - z;
  array<array<RealD, 4>, 3> eps_eps_xyz;
  set_zero(eps_eps_xyz);
  for (Int i = 0; i < 3; ++i) {
    for (Int rho = 0; rho < 4; ++rho) {
      for (Int j = 0; j < 3; ++j) {
        for (Int k = 0; k < 3; ++k) {
          for (Int m = 0; m < 4; ++m) {
            eps_eps_xyz[i][rho] += epsilon_tensor_acc(i, j, k) *
                                   epsilon_tensor_acc(rho, k, j, m) *
                                   x_mid_yz[m];
          }
        }
      }
    }
  }
  array<array<RealD, 4>, 4> eps_yz_xyz;
  set_zero(eps_yz_xyz);
  for (Int sigma = 0; sigma < 4; ++sigma) {
    for (Int lambda = 0; lambda < 4; ++lambda) {
      for (Int r = 0; r < 4; ++r) {
        for (Int s = 0; s < 4; ++s) {
          eps_yz_xyz[sigma][lambda] +=
              epsilon_tensor_acc(sigma, lambda, r, s) * y_z[r] * x_mid_yz[s];
        }
      }
    }
  }
  ManyMagneticMoments ret;
  set_zero(ret);
  for (Int i = 0; i < 3; ++i) {
    for (Int rho = 0; rho < 4; ++rho) {
      for (Int sigma = 0; sigma < 4; ++sigma) {
        for (Int lambda = 0; lambda < 4; ++lambda) {
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
  RealD sum = 0.0;
  for (Int i = 0; i < 3; ++i) {
    for (Int rho = 0; rho < 4; ++rho) {
      for (Int sigma = 0; sigma < 4; ++sigma) {
        for (Int lambda = 0; lambda < 4; ++lambda) {
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

void set_local_current_from_props(SelectedPoints<WilsonMatrix>& scf,
                                  const SelectedPoints<WilsonMatrix>& sprop1,
                                  const SelectedPoints<WilsonMatrix>& sprop2,
                                  const PointsSelection& psel_d,
                                  const Geometry& geo);

template <class M>
struct CurrentMoments {
  vector<array<M, 3 * 3>> d;
  // moment = 0.5 * epsilon_tensor(i, j, k) * smod_sym(xg[j] - ref[j],
  // total_site[j]) * d[ xg[j] ][3*j + k]
  //
  void init() { clear(d); }
  void init(const Int lsize)
  {
    init();
    d.resize(lsize);
    set_zero(d);
  }
};

template <class M>
qacc array<M, 3> simple_moment(const CurrentMoments<M>& cm,
                               const CoordinateD& ref,
                               const Coordinate& total_site)
{
  const Int lsize =
      std::max(total_site[0], std::max(total_site[1], total_site[2]));
  array<M, 3> ret;
  set_zero(ret);
  for (Int x = 0; x < lsize; ++x) {
    for (Int i = 0; i < 3; ++i) {
      for (Int j = 0; j < 3; ++j) {
        if (i == j) {
          continue;
        }
        for (Int k = 0; k < 3; ++k) {
          if (i == k or j == k) {
            continue;
          }
          ret[i] += (Complex)(0.5 * epsilon_tensor_acc(i, j, k) *
                              smod_sym(x - ref[j], (RealD)total_site[j])) *
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
    const Coordinate& total_site, const SelectedPoints<M>& current,
    const PointsSelection& psel_d, const SelectedPoints<RealD>& psel_d_prob_xy,
    const Long idx)
// subtract over weighting due to sparsening for both xg_z and xg_op.
{
  const Coordinate xg_op = psel_d[idx];
  const RealD prob = psel_d_prob_xy.get_elem(idx);
  const RealD weight = 1.0 / prob;
  const RealD sub_coef = 1.0 - weight;
  const Vector<M> cv = current.get_elems_const(idx);
  array<M, 3> ret = simple_moment(cm, ref, total_site);
  for (Int i = 0; i < 3; ++i) {
    for (Int j = 0; j < 3; ++j) {
      if (i == j) {
        continue;
      }
      for (Int k = 0; k < 3; ++k) {
        if (i == k or j == k) {
          continue;
        }
        ret[i] +=
            (sub_coef * 0.5 * (Complex)epsilon_tensor_acc(i, j, k) *
             (Complex)smod_sym(xg_op[j] - ref[j], (RealD)total_site[j])) *
            cv[k];
      }
    }
  }
  return ret;
}

enum struct ChooseReferenceLabel : Int {
  RefFar,
  RefClose,
  RefCenter,
};

inline ChooseReferenceLabel choose_reference_label(const std::string& label)
{
  if (does_string_have_tag(label, "ref-far")) {
    return ChooseReferenceLabel::RefFar;
  } else if (does_string_have_tag(label, "ref-close")) {
    return ChooseReferenceLabel::RefClose;
  } else if (does_string_have_tag(label, "ref-center")) {
    return ChooseReferenceLabel::RefCenter;
  } else {
    Qassert(false);
    return ChooseReferenceLabel::RefFar;
  }
  Qassert(false);
  return ChooseReferenceLabel::RefFar;
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
  if (ChooseReferenceLabel::RefFar == label) {
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
  } else if (ChooseReferenceLabel::RefClose == label) {
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
  } else if (ChooseReferenceLabel::RefCenter == label) {
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

RealD set_psel_d_prob_xy(SelectedPoints<RealD>& psel_d_prob_xy,
                         const PointsSelection& psel,
                         SelectedPoints<RealD>& psel_prob,
                         const PointsSelection& psel_d,
                         const SelectedPoints<RealD>& psel_d_prob,
                         const Long idx_xg_x, const Long idx_xg_y);

void set_current_moments_from_current(
    CurrentMoments<WilsonMatrix>& cm,
    const SelectedPoints<WilsonMatrix>& current, const PointsSelection& psel_d,
    const SelectedPoints<RealD>& psel_d_prob_xy, const Geometry& geo);

void glb_sum_current_moments(CurrentMoments<WilsonMatrix>& cm);

void contract_four_loop(SelectedPoints<Complex>& f_loop_i_rho_sigma_lambda,
                        const Complex& coef, const Coordinate& xg_x,
                        const Coordinate& xg_y,
                        const SelectedPoints<WilsonMatrix>& c_xy,
                        const SelectedPoints<WilsonMatrix>& c_yx,
                        const CurrentMoments<WilsonMatrix>& cm_xy,
                        const CurrentMoments<WilsonMatrix>& cm_yx,
                        const PointsSelection& psel_d,
                        const SelectedPoints<RealD>& psel_d_prob_xy,
                        const Geometry& geo,
                        const Long r_sq_limit, const std::string& label);

void contract_four_combine(
    SlTable& t, SlTable& t_pi, const Complex& coef, const Geometry& geo,
    const Coordinate& xg_x, const Coordinate& xg_y,
    const SelectedPoints<Complex>& f_loop_i_rho_sigma_lambda,
    const SelectedPoints<ManyMagneticMoments>& smf, const PointsSelection& psel_d,
    const Long r_sq_limit);

inline std::vector<std::string> get_clbl_inf_ref_tags(
    const std::string& job_tag)
{
  (void)job_tag;
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
  for (Int i = 0; i < (int)tags.size(); ++i) {
    std::string label = tags[i] + " proj-all";
    std::string label_pi = tags[i] + " proj-pi";
    labels.push_back(label);
    labels.push_back(label_pi);
  }
  return labels;
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
    const RealD muon_mass, const RealD z_v);

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
  for (Int i = 0; i < (int)tags.size(); ++i) {
    std::string label = tags[i] + " proj-all";
    std::string label_pi = tags[i] + " proj-pi";
    labels.push_back(label);
    labels.push_back(label_pi);
  }
  return labels;
}

std::vector<SlTable> contract_two_plus_two_pair_no_glb_sum(
    Long& n_points_selected, Long& n_points_computed, const ComplexD& coef,
    const Geometry& geo, const PointsSelection& psel,
    const SelectedPoints<RealD>& psel_prob, const PointsSelection& psel_lps,
    const SelectedPoints<RealD>& psel_lps_prob, const Long idx_xg_x,
    const SelectedPoints<ComplexD>& lps_hvp_x,
    const SelectedPoints<ComplexD>& edl_list_c, const Long r_sq_limit,
    const RealD muon_mass, const RealD z_v);

}  // namespace qlat
