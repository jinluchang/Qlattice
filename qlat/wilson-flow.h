#pragma once

#include <qlat/hmc.h>

namespace qlat
{  //

inline void set_wilson_flow_z(GaugeMomentum& z, const GaugeField& gf,
                              const double c1 = 0.0)
{
  TIMER("set_wilson_flow_z");
  const GaugeAction ga(3.0, c1);
  set_gm_force(z, gf, ga);
}

inline void gf_wilson_flow_step(GaugeField& gf, const double epsilon,
                                const double c1 = 0.0)
{
  TIMER("gf_wilson_flow_step");
  GaugeField& w = gf;
  GaugeMomentum z, zp;
  set_wilson_flow_z(z, w, c1);
  z *= 1.0 / 4.0;
  gf_evolve(w, z, epsilon);
  qswap(z, zp);
  zp *= 17.0 / 9.0;
  set_wilson_flow_z(z, w, c1);
  z *= 8.0 / 9.0;
  z -= zp;
  gf_evolve(w, z, epsilon);
  qswap(z, zp);
  set_wilson_flow_z(z, w, c1);
  z *= 3.0 / 4.0;
  z -= zp;
  gf_evolve(w, z, epsilon);
}

inline void gf_wilson_flow(GaugeField& gf, const double flow_time,
                           const int steps, const double c1 = 0.0)
{
  TIMER("gf_wilson_flow");
  const double epsilon = flow_time / (double)steps;
  for (int i = 0; i < steps; ++i) {
    gf_wilson_flow_step(gf, epsilon, c1);
  }
}

}  // namespace qlat
