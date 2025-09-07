#pragma once

#include <qlat/hmc.h>
#include <qlat/qcd-topology.h>

namespace qlat
{  //

void set_wilson_flow_z(GaugeMomentum& z, const GaugeField& gf,
                       const double c1 = 0.0);

void gf_wilson_flow_step_euler(GaugeField& gf, const double epsilon,
                               const double c1 = 0.0);

void gf_wilson_flow_step(GaugeField& gf, const double epsilon,
                         const double c1 = 0.0);

void gf_energy_density_dir_field(Field<RealD>& fd, const GaugeField& gf);

void gf_energy_density_field(Field<RealD>& fd, const GaugeField& gf);

RealD gf_energy_density(const GaugeField& gf);

std::vector<double> gf_wilson_flow(GaugeField& gf,
                                   const double existing_flow_time,
                                   const double flow_time, const int steps,
                                   const double c1 = 0.0);

void set_plaq_flow_z(GaugeMomentum& z, const GaugeField& gf,
                     const Field<RealD>& plaq_factor);

}  // namespace qlat
