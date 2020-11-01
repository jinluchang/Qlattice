#pragma once

#include "data-load.h"

namespace qlat
{  //

inline void set_local_va_current(Vector<WilsonMatrix> fc,
                                 const WilsonMatrix& prop_a,
                                 const WilsonMatrix& prop_b)
// ->- prop_a ->- op ->- inv_prop_b ->-
{
  const std::array<SpinMatrix, 8>& va_ms = get_va_matrices();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const WilsonMatrix inv_prop_b =
      gamma5 * (WilsonMatrix)matrix_adjoint(prop_b) * gamma5;
  for (int m = 0; m < 8; ++m) {
    fc[m] = inv_prop_b * va_ms[m] * prop_a;
  }
}

inline int tsep_op_wall_src(const std::string& job_tag)
// parameter
{
  if (job_tag == "24D" or job_tag == "32D" or job_tag == "24DH") {
    return 8;
  } else if (job_tag == "32Dfine") {
    return 10;
  } else if (job_tag == "48I") {
    return 12;
  } else if (job_tag == "64I") {
    return 18;
  } else {
    qassert(false);
  }
  return 8;
}

inline SelPropCache& get_prop_psrc_ama_cache()
{
  static SelPropCache cache("PropPsrcAmaCache", 8, 2);
  return cache;
}

inline const PselProp& get_psel_prop_psrc_ama(const std::string& job_tag,
                                              const int traj,
                                              const Coordinate& xg,
                                              const int type)
{
  TIMER_VERBOSE("get_psel_prop_psrc_ama");
  const std::vector<PointInfo>& pis_xgt =
      get_point_src_info(job_tag, traj, xg, type);
  if (pis_xgt.size() == 1) {
    qassert(pis_xgt[0].accuracy == 0);
    return get_psel_prop_psrc(job_tag, traj, xg, type, 0);
  }
  const std::string key = ssprintf("%s,%d,%s,%d,psrc-ama", show(xg).c_str(),
                                   type, job_tag.c_str(), traj);
  PselPropCache& cache = get_psel_prop_cache();
  if (not cache.has(key)) {
    const TypeAccuracyTable& tat = get_type_accuracy_table(job_tag, traj);
    qassert(get_accuracy_weight(tat, type, 0) == 1.0);
    const int num_acc = pis_xgt.size();
    qassert(num_acc > 1);
    std::vector<double> coefs(num_acc, 0.0);
    coefs[0] = 1.0;
    for (int acc = 1; acc < num_acc; ++acc) {
      const double weight = get_accuracy_weight(tat, type, acc);
      coefs[acc] += weight;
      coefs[acc - 1] -= weight;
    }
    PselProp& ps_prop = cache[key];
    for (int acc = 1; acc < num_acc; ++acc) {
      PselProp ps_prop_acc = get_psel_prop_psrc(job_tag, traj, xg, type, acc);
      ps_prop_acc *= coefs[acc];
      ps_prop += ps_prop_acc;
    }
  }
  return cache[key];
}

inline const SelProp& get_prop_psrc_ama(const std::string& job_tag,
                                        const int traj, const Coordinate& xg,
                                        const int type)
{
  TIMER_VERBOSE("get_prop_psrc_ama");
  const std::vector<PointInfo>& pis_xgt =
      get_point_src_info(job_tag, traj, xg, type);
  if (pis_xgt.size() == 1) {
    qassert(pis_xgt[0].accuracy == 0);
    return get_prop_psrc(job_tag, traj, xg, type, 0);
  }
  const std::string key = ssprintf("%s,%d,%s,%d,psrc-ama", show(xg).c_str(),
                                   type, job_tag.c_str(), traj);
  SelPropCache& cache = get_prop_psrc_ama_cache();
  if (not cache.has(key)) {
    const TypeAccuracyTable& tat = get_type_accuracy_table(job_tag, traj);
    qassert(get_accuracy_weight(tat, type, 0) == 1.0);
    const int num_acc = pis_xgt.size();
    qassert(num_acc > 1);
    std::vector<double> coefs(num_acc, 0.0);
    coefs[0] = 1.0;
    for (int acc = 1; acc < num_acc; ++acc) {
      const double weight = get_accuracy_weight(tat, type, acc);
      coefs[acc] += weight;
      coefs[acc - 1] -= weight;
    }
    SelProp& s_prop = cache[key];
    for (int acc = 1; acc < num_acc; ++acc) {
      SelProp s_prop_acc = get_prop_psrc(job_tag, traj, xg, type, acc);
      s_prop_acc *= coefs[acc];
      s_prop += s_prop_acc;
    }
    if (is_check_prop_consistency()) {
      TIMER_VERBOSE("check_prop_consistency_ama");
      const PselProp& ps_prop = get_psel_prop_psrc_ama(job_tag, traj, xg, type);
      const PointSelection& psel = get_point_selection(job_tag, traj);
      const FieldSelection& fsel = get_field_selection(job_tag, traj);
      const Geometry& geo = fsel.f_rank.geo;
      qassert(ps_prop.n_points == (long)psel.size());
      for (long ps_idx = 0; ps_idx < ps_prop.n_points; ++ps_idx) {
        const Coordinate& xg = psel[ps_idx];
        const Coordinate xl = geo.coordinate_l_from_g(xg);
        const long s_idx = fsel.f_local_idx.get_elem(xl);
        double qnorm_comm = 0.0;
        double qnorm_diff = 0.0;
        if (s_idx >= 0) {
          qnorm_comm += qnorm(s_prop.get_elem(s_idx));
          qnorm_diff +=
              qnorm(ps_prop.get_elem(ps_idx) - s_prop.get_elem(s_idx));
        }
        glb_sum(qnorm_comm);
        glb_sum(qnorm_diff);
        displayln_info(
            fname +
            ssprintf(": ps_prop diff qnorm = %24.17E. ps_prop qnorm = %24.17E",
                     qnorm_diff, qnorm_comm));
        qassert(qnorm_diff == 0.0);
      }
    }
  }
  return cache[key];
}

}  // namespace qlat
