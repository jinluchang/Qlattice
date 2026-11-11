#pragma once

#include <qlat-utils/matrix.h>
#include <qlat/qcd-utils.h>
#include <qlat/qcd.h>

namespace qlat
{  //

void gf_ape_smear(GaugeField& gf, const GaugeField& gf0, const double alpha,
                  const Long num_step = 1);

// APE smearing without the SU(3) projection, using the convex kernel
// `U' = (1 - alpha) U + (alpha / 6) * S` and multiplying the whole field by
// the real scalar `s = sqrt(3 / <N^2>)` after every step, where `<N^2>` is
// the global (MPI-reduced) average of `Tr(U_mu(x)^dag U_mu(x))` over all
// links, so that `<N^2> = 3` (the value for a unitary link) after every step.
// The rescale is skipped if the average is not positive.
// `alpha = 0` removes the staple mixing, but the rescale still applies.
// `num_step = 0` (or negative) returns an unchanged copy of `gf0`.
void gf_ape_smear_no_proj_rescaling(GaugeField& gf, const GaugeField& gf0,
                                    const RealD alpha, const Long num_step = 1);

// Same as `gf_ape_smear` but without the SU(3) projection step, using the
// step kernel `U' = U + (a / 6) * S` (no `(1 - a)` factor).
// `accumulated_alpha` is the target total effective smearing (input only).
// Each step uses `a = alpha * factor` with `factor = 3 / M`, where `M` is the
// global maximum of `Tr(U_mu(x)^dag U_mu(x))` over all links, so that `alpha`
// is the effective smearing step of the largest link. The field grows under
// this kernel (`U' ~ (1 + a) U` for near-unitary links), hence `M > 3` and
// `factor < 1`, so `a <= alpha`. The weight decays geometrically as the
// amplitude grows, so the accumulated smearing converges; a target above that
// limit is never reached and the loop stops on the no-progress guard,
// leaving a very large, non-unitary field. The field is never rescaled. The
// last step is truncated so that `sum(a)` does not exceed
// `accumulated_alpha`. The loop stops when `M` is degenerate, when `a` is
// not positive, or when `sum(a)` makes no progress in floating point.
// Progress is printed at verbose level (1) every 4 steps and on the last
// step. `alpha = 0` or `accumulated_alpha = 0` is the identity.
void gf_ape_smear_no_proj(GaugeField& gf, const GaugeField& gf0,
                          const RealD alpha, const RealD& accumulated_alpha);

void gf_ape_smear_spatial(GaugeField& gf, const GaugeField& gf0,
                          const double alpha, const Long num_step = 1);

void gf_hyp_smear(GaugeField& gf, const GaugeField& gf0, const double alpha1,
                  const double alpha2, const double alpha3,
                  const Long num_step = 1);

void prop_smear(Propagator4dT<RealD>& prop, const GaugeFieldT<RealD>& gf1,
                const double coef, const Int num_step,
                const CoordinateD& mom = CoordinateD(),
                const bool smear_in_time_dir = false);

void prop_smear_qlat_convention(Propagator4dT<RealD>& prop,
                                const GaugeFieldT<RealD>& gf, const double coef,
                                const Int num_step,
                                const CoordinateD& mom = CoordinateD(),
                                const bool smear_in_time_dir = false,
                                const Int mode = 1);

void prop_smear_qlat_convention(Propagator4dT<RealF>& prop,
                                const GaugeFieldT<RealF>& gf, const double coef,
                                const Int num_step,
                                const CoordinateD& mom = CoordinateD(),
                                const bool smear_in_time_dir = false,
                                const Int mode = 1);

void prop_spatial_smear_no_comm(std::vector<FermionField4d>& ff_vec,
                                const GaugeField& gf, const RealD coef,
                                const Long step,
                                const CoordinateD& mom = CoordinateD());

void gf_reduce_half(GaugeField& hgf, const GaugeField& gf);

}  // namespace qlat
