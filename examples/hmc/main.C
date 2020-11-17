#include "qlat-setup.h"

namespace qlat
{  //

inline void gm_evolve_fg(GaugeMomentum& gm, const GaugeField& gf_init,
                  const GaugeAction& ga, const double fg_dt, const double dt)
{
  TIMER("gm_evolve_fg");
  const Geometry& geo = gf_init.geo;
  GaugeField gf;
  gf.init(geo);
  gf = gf_init;
  //
  GaugeMomentum gm_force;
  gm_force.init(gm.geo);
  set_zero(gm_force);
  //
  set_gm_force(gm_force, gf, ga);
  //
  gf_evolve(gf, gm_force, fg_dt);
  //
  set_gm_force(gm_force, gf, ga);
  //
  // do the actual evolution
  gm_force *= dt;
  gm += gm_force;
}

inline void run_hmc(GaugeField& gf, const GaugeAction& ga, const int traj, const RngState& rs)
{
  TIMER_VERBOSE("run_hmc");
  const Geometry& geo = gf.geo;
  //
  GaugeField gf_init;
  gf_init = gf;
  //
  GaugeMomentum gm;
  gm.init(geo);
  set_rand_gauge_momentum(gm, 1.0, rs.split("set_rand_gauge_momentum"));
  //
  double energy = gm_hamilton_node(gm) + gf_hamilton_node(gf, ga);
  //
  int steps = 6;
  //
  double dt = 1.0 / steps;
  double lambda = 0.5 * (1.0 - 1.0 / sqrt(3.0));
  // double xi = 0.0;
  // double chi = 0.0;
  double theta = (2.0 - sqrt(3.0)) / 48.0;
  double ttheta = theta * dt * dt * dt;
  // double cchi = chi * dt * dt * dt;
  //
  gf_evolve(gf, gm, lambda * dt);
  for (int i = 0; i < steps; i++) {
    gm_evolve_fg(gm, gf, ga, 4.0 * ttheta / dt, 0.5 * dt);
    //
    gf_evolve(gf, gm, (1.0 - 2.0 * lambda) * dt);
    //
    gm_evolve_fg(gm, gf, ga, 4.0 * ttheta / dt, 0.5 * dt);
    //
    if (i < steps - 1) {
      gf_evolve(gf, gm, 2.0 * lambda * dt);
    } else {
      gf_evolve(gf, gm, lambda * dt);
    }
  }
  //
  unitarize(gf);
  //
  double delta_h = gm_hamilton_node(gm) + gf_hamilton_node(gf, ga) - energy;
  glb_sum(delta_h);
  //
  double accept_prob;
  bool flag =
      metropolis_accept(accept_prob, delta_h, rs.split("metropolis_accept"));
  displayln_info(
      fname +
      ssprintf(": accept flag = %d with prob accept = %.1f%% deltaH = %.16f traj = %d",
               flag, accept_prob * 100, delta_h, traj));
  //
  if (not flag and traj > 10) {
    displayln_info(fname + ssprintf(": restore gf (traj=%d).", traj));
    gf = gf_init;
  }
}

inline void test_hmc(const Coordinate& total_site)
{
  TIMER_VERBOSE("test_hmc");
  //
  Geometry geo;
  geo.init(total_site, 1);
  //
  GaugeAction ga(2.13, -0.331);
  // GaugeAction ga(5.5, 0.0);
  //
  const RngState rs =
      RngState(ssprintf("test_hmc-%s", show(total_site).c_str()));
  //
  GaugeField gf;
  gf.init(geo);
  set_unit(gf);
  //
  int traj = 0;
  //
  for (int i = 0; i < 30; i++) {
    traj++;
    run_hmc(gf, ga, traj, rs.split(ssprintf("hmc-%d", traj)));
    //
    const double plaq_avg = gf_avg_plaq(gf);
    const double plaq_sum = product(total_site) * 6 * (1.0 - plaq_avg);
    displayln_info(fname +
                   ssprintf(": traj=%d ; plaq_avg=%24.17E ; plaq_sum=%24.17E.",
                            plaq_avg, plaq_sum, traj));
  }
}

}  // namespace qlat

int main(int argc, char* argv[])
{
  using namespace qlat;
  Coordinate total_site(4, 4, 4, 8);
  // Coordinate total_site(8, 8, 8, 16);
  // Coordinate total_site(16, 16, 16, 16);
  begin(&argc, &argv);
  test_hmc(total_site);
  Timer::display();
  end();
  return 0;
}
