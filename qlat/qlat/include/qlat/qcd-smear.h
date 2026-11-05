#pragma once

#include <qlat-utils/matrix.h>
#include <qlat/qcd-utils.h>
#include <qlat/qcd.h>

namespace qlat
{  //

void gf_ape_smear(GaugeField& gf, const GaugeField& gf0, const double alpha,
                  const Long num_step = 1);

void gf_ape_smear_no_proj(GaugeField& gf, const GaugeField& gf0,
                          const RealD alpha, const Long num_step = 1);

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
