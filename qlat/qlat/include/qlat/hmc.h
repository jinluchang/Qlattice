#pragma once

#include <qlat/gauge-action.h>
#include <qlat/scalar-action.h>
#include <qlat/qcd-smear.h>
#include <qlat/qcd-topology.h>
#include <qlat/qcd-utils.h>
#include <qlat/qcd.h>

#include <cmath>
#include <sstream>
#include <string>

namespace qlat
{  //

struct API GaugeMomentum : FieldM<ColorMatrix, 4> {
};

bool metropolis_accept(RealD& accept_prob, const RealD delta_h,
                       const Int traj, const RngState& rs_);

void set_rand_gauge_momentum(GaugeMomentum& gm, const RealD sigma,
                             const RngState& rs);

void set_rand_gauge_momentum(GaugeMomentum& gm, const Field<RealD>& mf,
                             const RngState& rs);

RealD gm_hamilton_node(const GaugeMomentum& gm);

RealD gm_hamilton_node(const GaugeMomentum& gm, const Field<RealD>& mf);

RealD gf_sum_re_tr_plaq_node_no_comm(const GaugeField& gf);

RealD gf_sum_re_tr_rect_node_no_comm(const GaugeField& gf);

RealD gf_hamilton_node_no_comm(const GaugeField& gf, const GaugeAction& ga);

RealD gf_hamilton_node(const GaugeField& gf, const GaugeAction& ga);

RealD gf_hamilton(const GaugeField& gf, const GaugeAction& ga);

void gf_evolve(GaugeField& gf, const GaugeMomentum& gm, const RealD step_size);

void gf_evolve_dual(GaugeField& gf, const GaugeMomentum& gm_dual,
                    const RealD step_size);

void gf_evolve(GaugeField& gf, const GaugeMomentum& gm, const Field<RealD>& mf,
               const RealD step_size);

void gf_evolve_dual(GaugeField& gf, const GaugeMomentum& gm_dual,
                    const Field<RealD>& mf_dual, const RealD step_size);

void set_gm_force_no_comm(GaugeMomentum& gm_force, const GaugeField& gf,
                          const GaugeAction& ga);

void field_color_matrix_exp(Field<ColorMatrix>& fc,
                            const Field<ColorMatrix>& fc1, const ComplexD coef);

void field_color_matrix_mul(Field<ColorMatrix>& fc,
                            const Field<ColorMatrix>& fc1,
                            const Field<ColorMatrix>& fc2);

void set_gm_force(GaugeMomentum& gm_force, const GaugeField& gf,
                  const GaugeAction& ga);

void set_gm_force_dual(GaugeMomentum& gm_force_dual, const GaugeField& gf,
                       const GaugeMomentum& gm_force);

RealD project_gauge_transform(GaugeMomentum& gm, GaugeMomentum& gm_dual,
                              const Field<RealD>& mf,
                              const Field<RealD>& mf_dual);

void dot_gauge_momentum(Field<RealD>& f, const GaugeMomentum& gm1,
                        const GaugeMomentum& gm2);

void set_anti_hermitian_matrix_from_basis(Field<ColorMatrix>& fc,
                                          const Field<RealD>& basis);

void set_basis_from_anti_hermitian_matrix(Field<RealD>& basis,
                                          const Field<ColorMatrix>& fc);

// -------------------

}  // namespace qlat
