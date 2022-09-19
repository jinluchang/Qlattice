#include "qlat-setup.h"

namespace qlat
{  //

inline void gm_evolve_fg(GaugeMomentum& gm, const GaugeField& gf_init,
                  const GaugeAction& ga, const double fg_dt, const double dt)
{
  TIMER("gm_evolve_fg");
  const Geometry& geo = gf_init.geo();
  GaugeField gf;
  gf.init(geo);
  gf = gf_init;
  //
  GaugeMomentum gm_force;
  gm_force.init(gm.geo());
  //
  set_gm_force(gm_force, gf, ga);
  //
  gf_evolve(gf, gm_force, fg_dt);
  //
  set_gm_force(gm_force, gf, ga);
  //
  display_gm_force_magnitudes(gm_force, 5);
  //
  // do the actual evolution
  gm_force *= dt;
  gm += gm_force;
}

inline double run_hmc_evolve(GaugeMomentum& gm, GaugeField& gf,
                             const GaugeAction& ga, const RngState& rs,
                             const int steps, const double md_time = 1.0)
{
  TIMER_VERBOSE("run_hmc_evolve");
  //
  const double energy = gm_hamilton_node(gm) + gf_hamilton_node(gf, ga);
  //
  const double dt = md_time / steps;
  const double lambda = 0.5 * (1.0 - 1.0 / sqrt(3.0));
  const double theta = (2.0 - sqrt(3.0)) / 48.0;
  const double ttheta = theta * dt * dt * dt;
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
  return delta_h;
}

inline void run_hmc(GaugeField& gf, const GaugeAction& ga, const int traj, const RngState& rs)
{
  TIMER_VERBOSE("run_hmc");
  const Geometry& geo = gf.geo();
  //
  const bool is_reverse_test = traj < 3;
  //
  GaugeField gf0;
  gf0 = gf;
  //
  GaugeMomentum gm;
  gm.init(geo);
  set_rand_gauge_momentum(gm, 1.0, rs.split("set_rand_gauge_momentum"));
  //
  // ADJUST ME
  const int steps = 6;
  const double md_time = 1.0;
  //
  const double delta_h = run_hmc_evolve(gm, gf0, ga, rs, steps, md_time);
  //
  if (is_reverse_test) {
    GaugeMomentum gm_r;
    gm_r = gm;
    GaugeField gf0_r;
    gf0_r = gf0;
    const double delta_h_rev =
        run_hmc_evolve(gm_r, gf0_r, ga, rs, steps, -md_time);
    gf0_r -= gf;
    displayln_info(
        ssprintf("run_hmc_evolve reversed delta_diff: %24.17E / %24.17E",
                 delta_h + delta_h_rev, delta_h));
    displayln_info(
        ssprintf("run_hmc_evolve reversed gf_diff: %24.17E / %24.17E",
                 qnorm(gf0_r), qnorm(gf0)));
  }
  //
  double accept_prob;
  const bool flag = metropolis_accept(accept_prob, delta_h, traj,
                                      rs.split("metropolis_accept"));
  //
  if (flag or traj <= 20) {
    displayln_info(fname + ssprintf(": update gf (traj=%d).", traj));
    gf = gf0;
  }
}

inline void test_hmc(const Coordinate& total_site, const GaugeAction& ga)
{
  TIMER_VERBOSE("test_hmc");
  //
  qmkdir_info("results");
  qmkdir_info("results/gf_info");
  qmkdir_info("results/wilson_flow_energy_info");
  qmkdir_info("results/gm_force_info");
  //
  Geometry geo;
  geo.init(total_site, 1);
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
  for (int i = 0; i < 4; i++) {
    traj++;
    run_hmc(gf, ga, traj, rs.split(ssprintf("hmc-%d", traj)));
    //
    const double plaq_avg = gf_avg_plaq(gf);
    const double plaq_sum = product(total_site) * 6 * (1.0 - plaq_avg);
    displayln_info("CHECK: " + fname +
                   ssprintf(": traj=%d ; plaq_avg=%24.13E ; plaq_sum=%24.13E.",
                            traj, plaq_avg, plaq_sum));
    //
    if (traj % 2 == 0) {
      display_gauge_field_info_table_with_wilson_flow(
          ssprintf("results/gf_info/traj=%d.lat", traj),
          ssprintf("results/wilson_flow_energy_info/traj=%d.lat", traj), gf,
          0.1, 5, 2);
      save_gm_force_magnitudes_list(
          ssprintf("results/gm_force_info/traj=%d.lat", traj));
    }
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
  ColorMatrixConstants::get_instance().check();
  const GaugeAction ga1(2.13, -0.331);
  const GaugeAction ga2(5.5, 0.0);
  test_hmc(total_site, ga1);
  test_hmc(total_site, ga2);
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
