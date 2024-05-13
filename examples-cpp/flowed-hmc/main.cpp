#include "qlat/flowed-hmc.h"
#include "qlat-setup.h"

namespace qlat
{  //

inline void gm_evolve_fg(GaugeMomentum& gm, const GaugeField& gf_init,
                         const GaugeAction& ga, const FlowInfo& fi,
                         const double fg_dt, const double dt)
{
  TIMER("gm_evolve_fg");
  const Geometry& geo = gf_init.geo();
  GaugeField gf;
  gf.init(geo);
  gf = gf_init;
  //
  GaugeMomentum gm_force;
  //
  set_gm_force_flowed(gm_force, gf, ga, fi);
  //
  gf_evolve(gf, gm_force, fg_dt);
  //
  set_gm_force_flowed(gm_force, gf, ga, fi);
  //
  display_gm_force_magnitudes(gm_force, 5);
  //
  // do the actual evolution
  gm_force *= dt;
  gm += gm_force;
}

inline double run_hmc_evolve_flowed(GaugeMomentum& gm, GaugeField& gf,
                                    const GaugeAction& ga, const FlowInfo& fi,
                                    const RngState& rs, const int steps,
                                    const double md_time = 1.0)
{
  TIMER_VERBOSE("run_hmc_evolve_flowed");
  //
  const double energy =
      gm_hamilton_node(gm) + gf_hamilton_flowed_node(gf, ga, fi);
  //
  const double dt = md_time / steps;
  const double lambda = 0.5 * (1.0 - 1.0 / sqrt(3.0));
  // double xi = 0.0;
  // double chi = 0.0;
  const double theta = (2.0 - sqrt(3.0)) / 48.0;
  const double ttheta = theta * dt * dt * dt;
  // double cchi = chi * dt * dt * dt;
  //
  gf_evolve(gf, gm, lambda * dt);
  for (int i = 0; i < steps; i++) {
    gm_evolve_fg(gm, gf, ga, fi, 4.0 * ttheta / dt, 0.5 * dt);
    //
    gf_evolve(gf, gm, (1.0 - 2.0 * lambda) * dt);
    //
    gm_evolve_fg(gm, gf, ga, fi, 4.0 * ttheta / dt, 0.5 * dt);
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
  double delta_h =
      gm_hamilton_node(gm) + gf_hamilton_flowed_node(gf, ga, fi) - energy;
  glb_sum(delta_h);
  //
  return delta_h;
}

inline FlowInfo mk_flow_info(const RngState& rs)
{
  FlowInfo fi;
  // ADJUST ME
  vector_append(fi.v, mk_flow_info_step(rs, 0.1, -0.01).v);
  //
  return fi;
}

inline void run_hmc(GaugeField& gf, const GaugeAction& ga, const int traj,
                    const RngState& rs)
{
  TIMER_VERBOSE("run_hmc");
  const Geometry& geo = gf.geo();
  //
  const bool is_reverse_test = traj < 3;
  //
  const FlowInfo fi = mk_flow_info(rs.split("mk_flow_info"));
  displayln_info(show(fi));
  //
  GaugeField gf0;
  gf_flow_inv(gf0, gf, fi);
  //
  if (is_reverse_test) {
    GaugeField gf_r;
    gf_flow(gf_r, gf0, fi);
    gf_r -= gf;
    displayln_info(ssprintf("gf_flow_inv gf_diff: %24.17E / %24.17E",
                            qnorm(gf_r), qnorm(gf)));
  }
  //
  GaugeField gf0_init;
  if (is_reverse_test) {
    gf0_init = gf0;
  }
  //
  GaugeMomentum gm;
  gm.init(geo);
  set_rand_gauge_momentum(gm, 1.0, rs.split("set_rand_gauge_momentum"));
  //
  // ADJUST ME
  const int steps = 6;
  const double md_time = 1.0;
  //
  const double delta_h =
      run_hmc_evolve_flowed(gm, gf0, ga, fi, rs, steps, md_time);
  //
  const double plaq_avg = gf_avg_plaq(gf0);
  const double plaq_sum = product(geo.total_site()) * 6 * (1.0 - plaq_avg);
  displayln_info(
      fname + ssprintf(": U_0: traj=%d ; plaq_avg=%24.17E ; plaq_sum=%24.17E .",
                       traj, plaq_avg, plaq_sum));
  //
  if (is_reverse_test) {
    GaugeMomentum gm_r;
    gm_r = gm;
    GaugeField gf0_r;
    gf0_r = gf0;
    const double delta_h_rev =
        run_hmc_evolve_flowed(gm_r, gf0_r, ga, fi, rs, steps, -md_time);
    gf0_r -= gf0_init;
    displayln_info(
        ssprintf("run_hmc_evolve_flowed reversed delta_diff: %24.17E / %24.17E",
                 delta_h + delta_h_rev, delta_h));
    displayln_info(
        ssprintf("run_hmc_evolve_flowed reversed gf_diff: %24.17E / %24.17E",
                 qnorm(gf0_r), qnorm(gf0)));
  }
  //
  double accept_prob;
  const bool flag = metropolis_accept(accept_prob, delta_h, traj,
                                      rs.split("metropolis_accept"));
  //
  if (flag or traj <= 20) {
    displayln_info(fname + ssprintf(": update gf (traj=%d).", traj));
    gf_flow(gf, gf0, fi);
  }
}

inline std::string get_config_fn(const int traj)
{
  return ssprintf("results/configs/ckpoint_lat.%d", traj);
}

inline void test_hmc(const Coordinate& total_site)
{
  TIMER_VERBOSE("test_hmc");
  //
  qmkdir_info("results");
  qmkdir_info("results/gf_info");
  qmkdir_info("results/wilson_flow_energy_info");
  qmkdir_info("results/gm_force_info");
  qmkdir_info("results/configs");
  //
  // ADJUST ME
  // const int traj_max = 100000;
  const int traj_max = 3;
  //
  Geometry geo;
  geo.init(total_site);
  //
  // ADJUST ME
  // GaugeAction ga(1.70, -0.331);
  // GaugeAction ga(2.13, -0.331); // 1.730
  // GaugeAction ga(2.25, -0.331); // 2.359
  GaugeAction ga(2.31, -0.331);  // 2.774
  // GaugeAction ga(5.5, 0.0);
  // GaugeAction ga(0.7796, -1.4088);  // 0.2000 fm
  // GaugeAction ga(1.0038, -1.4088);  // 0.1000 fm
  // GaugeAction ga(0.91082, -1.4088);  // 0.1250 fm
  // GaugeAction ga(0.89, -1.4088);  // 0.1324 fm
  //
  const RngState rs =
      RngState(ssprintf("test_hmc-%s", show(total_site).c_str()));
  //
  GaugeField gf;
  gf.init(geo);
  set_unit(gf);
  // ADJUST ME
  set_g_rand_color_matrix_field(gf, RngState(), 1.0);
  //
  int traj = 0;
  //
  // ADJUST ME
  const int traj_save_skip = 10;
  //
  for (int i = 0; i <= traj_max; i += traj_save_skip) {
    if (does_file_exist_sync_node(get_config_fn(i))) {
      traj = i;
    }
  }
  //
  if (traj != 0) {
    load_gauge_field(gf, get_config_fn(traj));
  }
  //
  for (int i = traj; i < traj_max; ++i) {
    check_sigterm();
    check_stop();
    //
    traj++;
    run_hmc(gf, ga, traj, rs.split(ssprintf("hmc-%d", traj)));
    //
    const double plaq_avg = gf_avg_plaq(gf);
    const double plaq_sum = product(total_site) * 6 * (1.0 - plaq_avg);
    displayln_info("CHECK: " + fname +
                   ssprintf(": traj=%d ; plaq_avg=%24.12E ; plaq_sum=%24.12E.",
                            traj, plaq_avg, plaq_sum));
    //
    // ADJUST ME
    if (traj % 1 == 0) {
      display_gauge_field_info_table_with_wilson_flow(
          ssprintf("results/gf_info/traj=%d.lat", traj),
          ssprintf("results/wilson_flow_energy_info/traj=%d.lat", traj), gf,
          1.0, 100, 2);
    }
    //
    if (traj % traj_save_skip == 0) {
      save_gm_force_magnitudes_list(
          ssprintf("results/gm_force_info/traj=%d.lat", traj));
      save_gauge_field(gf, get_config_fn(traj));
      check_time_limit();
    }
    Timer::autodisplay();
  }
}

inline void test_color_matrix()
{
  TIMER_VERBOSE("test_color_matrix");
  ColorMatrixConstants::get_instance().check();
}

}  // namespace qlat

int main(int argc, char* argv[])
{
  using namespace qlat;
  // ADJUST ME
  Coordinate total_site(4, 4, 4, 8);
  // Coordinate total_site(8, 8, 8, 8);
  // Coordinate total_site(8, 8, 8, 16);
  // Coordinate total_site(16, 16, 16, 32);
  //
  begin(&argc, &argv);
  setup();
  test_color_matrix();
  test_hmc(total_site);
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
