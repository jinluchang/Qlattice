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
// Runge–Kutta scheme
// http://arxiv.org/abs/1006.4518v3
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

inline double gf_energy_density_no_comm(const GaugeField& gf)
{
  TIMER("gf_energy_density_no_comm");
  const Geometry geo = geo_reform(gf.geo());
  FieldM<double, 1> fd;
  fd.init(geo);
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    double s = 0.0;
    for (int mu = 0; mu < 3; ++mu) {
      for (int nu = mu + 1; nu < 4; ++nu) {
        const ColorMatrix g_mu_nu = make_tr_less_anti_herm_matrix(
            gf_clover_leaf_no_comm(gf, xl, mu, nu));
        s += -matrix_trace(g_mu_nu, g_mu_nu).real();
      }
    }
    fd.get_elem(index) = s;
  });
  double sum = 0.0;
  for (long index = 0; index < geo.local_volume(); ++index) {
    sum += fd.get_elem(index);
  }
  glb_sum(sum);
  sum *= 1.0 / (double)geo.total_volume();
  return sum;
}

inline double gf_energy_density(const GaugeField& gf)
// https://arxiv.org/pdf/1006.4518.pdf Eq. (2.1) (Fig. 1) (approximate Eq. (3.1))
// https://arxiv.org/pdf/1203.4469.pdf
{
  TIMER("gf_energy_density");
  GaugeField gf1;
  gf1.init(geo_resize(gf.geo(), 1));
  gf1 = gf;
  refresh_expanded(gf1);
  return gf_energy_density_no_comm(gf1);
}

inline std::vector<double> gf_wilson_flow(GaugeField& gf,
                                          const double existing_flow_time,
                                          const double flow_time,
                                          const int steps,
                                          const double c1 = 0.0)
{
  TIMER("gf_wilson_flow");
  std::vector<double> energy_density_list(steps, 0.0);
  const double epsilon = flow_time / (double)steps;
  for (int i = 0; i < steps; ++i) {
    gf_wilson_flow_step(gf, epsilon, c1);
    const double t = (i + 1) * epsilon + existing_flow_time;
    const double energy_density = gf_energy_density(gf);
    energy_density_list[i] = energy_density;
    displayln_info(fname +
                   ssprintf(": t = %24.17E ; E = %24.17E ; t^2 E = %24.17E.", t,
                            energy_density, sqr(t) * energy_density));
  }
  return energy_density_list;
}

}  // namespace qlat
