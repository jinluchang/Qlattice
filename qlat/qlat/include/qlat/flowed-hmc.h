#pragma once

#include <qlat/field-expand.h>
#include <qlat/hmc.h>

namespace qlat
{  //

// mask -> mask
// flow_size -> flow_type

// Masking scheme: mask

// get_mask_node_site: node_site -> mask_node_site
//
// get_num_mask : flow_type -> num_mask
//
// coordinate_from_mask_coordinate : mxg, mask, flow_type -> xg
//
// mask_from_coordinate : xg, flow_type -> mask
//
// mask_coordinate_from_coordinate : xg, flow_type -> mxg

struct FlowStepInfo {
  Int mask;
  Int mu;
  double epsilon;
  Int flow_size;
  //
  FlowStepInfo() {}
  FlowStepInfo(const Int mask_, const Int mu_, const double epsilon_,
               const Int flow_size_ = 1)
  {
    mask = mask_;
    mu = mu_;
    epsilon = epsilon_;
    flow_size = flow_size_;
  }
};

struct FlowInfo {
  std::vector<FlowStepInfo> v;
};

std::string show(const FlowInfo& fi);

FlowInfo mk_flow_info_step(const RngState& rs, const double epsilon);

FlowInfo mk_flow_info_step(const RngState& rs, const double epsilon,
                           const double epsilon2);

void gf_flow(GaugeField& gf, const GaugeField& gf0, const FlowInfo& fi);

void gf_flow_inv(GaugeField& gf, const GaugeField& gf1, const FlowInfo& fi);

double gf_hamilton_flowed_node(const GaugeField& gf0, const GaugeAction& ga,
                               const FlowInfo& fi);

void set_gm_force_flowed(GaugeMomentum& gm_force, const GaugeField& gf0,
                         const GaugeAction& ga, const FlowInfo& fi);

void set_gm_force_flowed_no_det(GaugeMomentum& gm_force,
                                const GaugeMomentum& gm_force_pre,
                                const GaugeField& gf0, const FlowInfo& fi);

double gf_flow_and_ln_det_node(GaugeField& gf, const GaugeField& gf0,
                               const FlowInfo& fi);

void set_flowed_gauge_fields(std::vector<GaugeField>& gf_ext_vec,
                             const GaugeField& gf0, const FlowInfo& fi);

void set_gm_force_propagated_det_from_flow(
    GaugeMomentum& gm_force, const GaugeMomentum& gm_force_pre,
    const std::vector<GaugeField>& gf_ext_vec, const FlowInfo& fi);

void set_gm_force_propagated_and_gm_force_det_from_flow_step(
    GaugeMomentum& gm_force_propagated, GaugeMomentum& gm_force_det,
    const GaugeMomentum& gm_force_pre, const GaugeField& gf0_ext,
    const FlowStepInfo& fsi);

void set_gm_force_propagated_no_det_from_flow(
    GaugeMomentum& gm_force, const GaugeMomentum& gm_force_pre,
    const std::vector<GaugeField>& gf_ext_vec, const FlowInfo& fi);

void set_gm_force_propagated_from_flow_step(GaugeMomentum& gm_force_propagated,
                                            const GaugeMomentum& gm_force_pre,
                                            const GaugeField& gf0_ext,
                                            const FlowStepInfo& fsi);

}  // namespace qlat
