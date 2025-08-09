#include <qlat/qcd-acc.h>
#include <qlat/wilson-flow.h>

namespace qlat
{  //

void set_wilson_flow_z(GaugeMomentum& z, const GaugeField& gf, const double c1)
{
  TIMER("set_wilson_flow_z");
  const GaugeAction ga(3.0, c1);
  set_gm_force(z, gf, ga);
}

void gf_wilson_flow_step_euler(GaugeField& gf, const double epsilon,
                               const double c1)
{
  TIMER("gf_wilson_flow_step_euler");
  GaugeField& w = gf;
  GaugeMomentum z;
  set_wilson_flow_z(z, w, c1);
  gf_evolve(w, z, epsilon);
}

void gf_wilson_flow_step(GaugeField& gf, const double epsilon, const double c1)
// Runge-Kutta scheme
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

static void gf_energy_density_field_no_comm(Field<RealD>& fd,
                                            const GaugeField& gf)
{
  TIMER("gf_energy_density_field_no_comm");
  const Geometry geo = geo_resize(gf.geo());
  fd.init(geo, 1);
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
}

void gf_energy_density_field(Field<RealD>& fd, const GaugeField& gf)
// https://arxiv.org/pdf/1006.4518.pdf Eq. (2.1) (Fig. 1) (approximate Eq. (3.1))
// https://arxiv.org/pdf/1203.4469.pdf
{
  TIMER("gf_energy_density_field");
  GaugeField gf1;
  gf1.init(geo_resize(gf.geo(), 1));
  gf1 = gf;
  refresh_expanded(gf1);
  gf_energy_density_field_no_comm(fd, gf1);
}

RealD gf_energy_density(const GaugeField& gf)
// https://arxiv.org/pdf/1006.4518.pdf Eq. (2.1) (Fig. 1) (approximate Eq. (3.1))
// https://arxiv.org/pdf/1203.4469.pdf
{
  TIMER("gf_energy_density");
  const Geometry& geo = gf.geo();
  FieldM<RealD, 1> fd;
  gf_energy_density_field(fd, gf);
  return field_glb_sum(fd)[0] / (RealD)geo.total_volume();
}

std::vector<double> gf_wilson_flow(GaugeField& gf,
                                   const double existing_flow_time,
                                   const double flow_time, const int steps,
                                   const double c1)
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
