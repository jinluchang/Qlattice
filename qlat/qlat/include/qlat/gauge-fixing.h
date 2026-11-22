#pragma once

#include <qlat/qlat.h>

namespace qlat
{

// Plain local steepest-descent gauge fixing (Landau and Coulomb) on a
// distributed gauge field.
//
// The gauge functional is
//
//   F = sum_{x, mu in dirs} Re Tr U_mu(x) / (V * len(dirs) * 3)
//
// and the local update is
//
//   g(x) = exp(alpha * D(x))
//
// where D(x) is the anti-hermitian traceless part of the staple-sum matrix
// K(x) = sum_{mu in dirs} (U_mu(x) + U_mu(x-mu)^dagger), i.e. the local
// ascent direction of the functional.  The gauge field is then updated by
//
//   U_mu(x) -> g(x) U_mu(x) g(x+mu)^dagger.
//
// One sweep is a Jacobi sweep over all sites.  The loop stops when the RMS
// divergence theta (which vanishes exactly in the gauge) is at most eps or
// after max_iter sweeps.  There is no Fourier acceleration and no FFT.

// Convergence history and final values of a gauge-fixing run.
struct GaugeFixingInfo {
  bool converged;
  Int n_iter;
  Int max_iter;
  RealD eps;
  RealD alpha;
  RealD functional;
  RealD divergence;
  std::vector<RealD> functional_list;
  std::vector<RealD> divergence_list;
};

// Normalized gauge functional:
//   F = sum_{x, mu in dirs} Re Tr U_mu(x) / (V * len(dirs) * 3)
RealD gf_gauge_functional(const GaugeField& gf, const std::vector<Int>& dirs);

// Local gradient field D(x) of the gauge functional (anti-hermitian,
// traceless).  The gauge condition is D(x) = 0.
void gf_gauge_gradient(Field<ColorMatrix>& D, const GaugeField& gf,
                       const std::vector<Int>& dirs);

// RMS divergence theta = sqrt(sum_x Tr[D(x) D(x)^dagger] / (V * 3))
RealD gf_gauge_divergence_from_gradient(const Field<ColorMatrix>& D);

RealD gf_gauge_divergence(const GaugeField& gf, const std::vector<Int>& dirs);

// Maximal |(U U^dagger - 1)_{ij}| over all links and sites.
RealD gf_gauge_unitarity_error(const GaugeField& gf);

GaugeFixingInfo gf_fix_gauge(GaugeField& gf, const std::vector<Int>& dirs,
                             const Int max_iter, const RealD eps,
                             const RealD alpha, const bool verbose);

GaugeFixingInfo gf_fix_landau_gauge(GaugeField& gf, const Int max_iter = 10000,
                                    const RealD eps = 1e-10,
                                    const RealD alpha = 0.1,
                                    const bool verbose = false);

GaugeFixingInfo gf_fix_coulomb_gauge(GaugeField& gf, const Int max_iter = 10000,
                                     const RealD eps = 1e-10,
                                     const RealD alpha = 0.1,
                                     const bool verbose = false);

}  // namespace qlat
