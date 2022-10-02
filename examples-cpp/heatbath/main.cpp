#include <qlat/qlat.h>

#include "qlat-setup.h"

namespace qlat
{  //

struct CorrParams {
  int t1, t2, dt;
  //
  CorrParams(const int t1_, const int t2_, const int dt_) {
    t1 = t1_;
    t2 = t2_;
    dt = dt_;
  }
};

struct CorrFuncs {
  double phi2;
  double c2_0;
  double c2_t1, c2_t2;
  double c4_t1, c4_t2;
  //
  void init()
  {
    phi2 = 0.0;
    c2_0 = 0.0;
    c2_t1 = 0.0;
    c2_t2 = 0.0;
    c4_t1 = 0.0;
    c4_t2 = 0.0;
  }
  //
  CorrFuncs() { init(); }
  //
  const CorrFuncs& operator+=(const CorrFuncs& cf) {
    phi2 += cf.phi2;
    c2_0 += cf.c2_0;
    c2_t1 += cf.c2_t1;
    c2_t2 += cf.c2_t2;
    c4_t1 += cf.c4_t1;
    c4_t2 += cf.c4_t2;
    return *this;
  }
  //
  const CorrFuncs& operator*=(const double coef) {
    phi2 *= coef;
    c2_0 *= coef;
    c2_t1 *= coef;
    c2_t2 *= coef;
    c4_t1 *= coef;
    c4_t2 *= coef;
    return *this;
  }
};

struct Observables
{
  double phi2, m_eff, v_eff;
  //
  void init()
  {
    phi2 = 0.0;
    m_eff = 0.0;
    v_eff = 0.0;
  }
  //
  Observables() { init(); }
  //
  const Observables& operator+=(const Observables& obs) {
    phi2 += obs.phi2;
    m_eff += obs.m_eff;
    v_eff += obs.v_eff;
    return *this;
  }
  //
  const Observables& operator*=(const double coef) {
    phi2 *= coef;
    m_eff *= coef;
    v_eff *= coef;
    return *this;
  }
};

typedef FieldM<double, 1> ScalarField;

inline void set_rng_field(RngField& rf, const Coordinate& total_site,
                          const std::string& seed, const long traj)
{
  TIMER("set_rng_field");
  Geometry geo;
  geo.init(total_site, 1);
  rf.init();
  rf.init(geo, RngState(seed).split(traj));
}

inline void set_scalar_field(ScalarField& sf, const Coordinate& total_site)
{
  TIMER("set_scalar_field");
  Geometry geo;
  geo.init(total_site, 1);
  geo.resize(1);
  sf.init();
  sf.init(geo);
  set_zero(sf);
}

inline void refresh_scalar_field(ScalarField& sf)
{
  TIMER("refresh_scalar_field");
  refresh_expanded_1(sf);
}

inline double action_scalar_field_site(const double val,
                                       const double nearby_sum, const double k1,
                                       const double k2)
{
  const double val2 = sqr(val);
  return -nearby_sum * val + k1 * val2 + k2 * sqr(val2);
}

inline void update_scalar_field_site(double& val, const double nearby_sum,
                                     RngState& rs, const double k1,
                                     const double k2)
{
  // ADJUST ME
  const int max_iter = 5;
  const double sigma = 1.0;
  //
  double action = action_scalar_field_site(val, nearby_sum, k1, k2);
  for (int i = 0; i < max_iter; ++i) {
    const double val_new = g_rand_gen(rs, val, sigma);
    const double action_new =
        action_scalar_field_site(val_new, nearby_sum, k1, k2);
    const double action_diff = action_new - action;
    if (action_diff <= 0 or std::exp(-action_diff) >= u_rand_gen(rs)) {
      val = val_new;
      action = action_new;
    }
  }
}

inline void sweep_scalar_field(ScalarField& sf, RngField& rf,
                               const int eo, const double mass_sqr,
                               const double lambda)
{
  TIMER("sweep_scalar_field(sf,rf,eo,m2,lam)");
  const Geometry& geo = sf.geo();
  const double k1 = 4.0 + 0.5 * mass_sqr;
  const double k2 = 1.0 / 24.0 * lambda;
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    if ((xg[0] + xg[1] + xg[2] + xg[3]) % 2 == 2 - eo) {
      RngState& rs = rf.get_elem(xl);
      double nearby_sum = 0.0;
      for (int dir = -4; dir < 4; ++dir) {
        const Coordinate xl_shifted = coordinate_shifts(xl, dir);
        nearby_sum += sf.get_elem(xl_shifted);
      }
      update_scalar_field_site(sf.get_elem(xl), nearby_sum, rs, k1, k2);
    }
  }
}

inline void sweep_scalar_field(ScalarField& sf, RngField& rf,
                               const double mass_sqr, const double lambda)
// interface function
{
  TIMER("sweep_scalar_field");
  refresh_scalar_field(sf);
  sweep_scalar_field(sf, rf, 1, mass_sqr, lambda);
  refresh_scalar_field(sf);
  sweep_scalar_field(sf, rf, 2, mass_sqr, lambda);
}

inline CorrFuncs measure_corr_funcs(const ScalarField& sf, const CorrParams& cp)
{
  TIMER("measure_corr_funcs");
  const Geometry& geo = sf.geo();
  const Coordinate total_site = geo.total_site();
  std::vector<double> phi_ts(total_site[3], 0.0);
  std::vector<double> phi_ts2(total_site[3], 0.0);
  double phi2 = 0.0;
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const double val = sf.get_elem(xl);
    phi_ts[xg[3]] += val;
    phi2 += sqr(val);
  }
  glb_sum(phi2);
  glb_sum(get_data(phi_ts));
  phi_ts2 = phi_ts;
  for (int dt = 1; dt < cp.dt; ++dt) {
    for (int t = 0; t < total_site[3]; ++t) {
      phi_ts[t] += phi_ts2[mod(t + dt, total_site[3])];
    }
  }
  for (int t = 0; t < total_site[3]; ++t) {
    phi_ts2[t] = sqr(phi_ts[t]);
  }
  CorrFuncs cf;
  for (int t = 0; t < total_site[3]; ++t) {
    cf.c2_0 += sqr(phi_ts[t]);
    cf.c2_t1 += phi_ts[t] * phi_ts[mod(t + cp.t1, total_site[3])];
    cf.c2_t2 += phi_ts[t] * phi_ts[mod(t + cp.t2, total_site[3])];
    cf.c4_t1 += phi_ts2[t] * phi_ts2[mod(t + cp.t1, total_site[3])];
    cf.c4_t2 += phi_ts2[t] * phi_ts2[mod(t + cp.t2, total_site[3])];
  }
  cf *= 1.0 / (double)total_site[3];
  cf.phi2 = phi2 / (double)product(total_site);
  // displayln_info(show(cf.phi2));
  return cf;
}

inline Observables get_observables(const CorrFuncs& cf, const CorrParams& cp)
{
  TIMER("get_observables");
  const double r4_t1 = (cf.c4_t1 - sqr(cf.c2_0)) / sqr(cf.c2_t1);
  const double r4_t2 = (cf.c4_t2 - sqr(cf.c2_0)) / sqr(cf.c2_t2);
  // displayln_info(fname + ssprintf(": r4_t1=%.8lf ; r4_t2=%.8lf.", r4_t1, r4_t2));
  Observables obs;
  obs.phi2 = cf.phi2;
  obs.m_eff = std::log(cf.c2_t1 / cf.c2_t2) / (double)(cp.t2 - cp.t1);
  obs.v_eff = std::log(r4_t1 / r4_t2) / (double)(cp.t2 - cp.t1);
  return obs;
}

inline std::string show_result(const long traj, const CorrFuncs& cf)
{
  TIMER("show_result");
  std::ostringstream out;
  out << ssprintf("((traj %s)", show(traj).c_str());
  out << std::endl;
  out << ssprintf(" (phi2 %s)", show(cf.phi2).c_str());
  out << std::endl;
  out << ssprintf(" (c2_0 %s)", show(cf.c2_0).c_str());
  out << std::endl;
  out << ssprintf(" (c2_t1 %s)", show(cf.c2_t1).c_str());
  out << std::endl;
  out << ssprintf(" (c2_t2 %s)", show(cf.c2_t2).c_str());
  out << std::endl;
  out << ssprintf(" (c4_t1 %s)", show(cf.c4_t1).c_str());
  out << std::endl;
  out << ssprintf(" (c4_t2 %s)", show(cf.c4_t2).c_str());
  out << ")";
  out << std::endl;
  return out.str();
}

inline std::string show_results(const std::vector<CorrFuncs>& cfs,
                                const Coordinate& total_site,
                                const CorrParams& cp, const double mass_sqr,
                                const double lambda)
{
  TIMER("show_results");
  std::ostringstream out;
  out << ssprintf("((total-site (%d %d %d %d))", total_site[0], total_site[1], total_site[2], total_site[3]);
  out << std::endl;
  out << ssprintf(" (mass-sqr %.10lf)", mass_sqr);
  out << std::endl;
  out << ssprintf(" (lambda   %.10lf)", lambda);
  const long skip_traj = (long)cfs.size() / 3;
  std::vector<CorrFuncs> cfs_d = vector_drop(cfs, skip_traj);
  out << std::endl;
  out << ssprintf(" (n-traj      %ld)", cfs.size());
  out << std::endl;
  out << ssprintf(" (n-traj-used %ld)", cfs_d.size());
  std::vector<long> n_block_list;
  n_block_list.push_back(1024 * 1024);
  n_block_list.push_back(128);
  n_block_list.push_back(64);
  n_block_list.push_back(32);
  for (int k = 0; k < (int)n_block_list.size(); ++k) {
    const long n_block = n_block_list[k];
    std::vector<CorrFuncs> cfs_db = vector_block(cfs_d, n_block);
    std::vector<CorrFuncs> jcfs = jackknife(cfs_db);
    const long jsize = jcfs.size();
    std::vector<Observables> jobs(jsize);
    for (int i = 0; i < jsize; ++i) {
      jobs[i] = get_observables(jcfs[i], cp);
    }
    out << std::endl;
    out << ssprintf(" (n-block-%ld", n_block);
    std::vector<double> vals(jsize);
    for (int i = 0; i < jsize; ++i) {
      vals[i] = jobs[i].phi2;
    }
    out << std::endl;
    out << ssprintf("  (phi2  %.15lf %.15f)", vals[0], jackknife_sigma(vals));
    for (int i = 0; i < jsize; ++i) {
      vals[i] = jobs[i].m_eff;
    }
    out << std::endl;
    out << ssprintf("  (m_eff %.15lf %.15f)", vals[0], jackknife_sigma(vals));
    for (int i = 0; i < jsize; ++i) {
      vals[i] = jobs[i].v_eff;
    }
    out << std::endl;
    out << ssprintf("  (v_eff %.15lf %.15f)", vals[0], jackknife_sigma(vals));
    out << ssprintf(")");
  }
  out << ssprintf(")");
  out << std::endl;
  return out.str();
}

inline void evolution(const Coordinate& total_site, const CorrParams& cp,
                      const double mass_sqr, const double lambda)
// interface function
{
  const std::string fn =
      ssprintf("results/total_site=%s/lambda=%.10lf/mass_sqr=%.10lf.txt",
               show(total_site).c_str(), lambda, mass_sqr);
  if (does_file_exist_sync_node(fn)) {
    return;
  }
  qmkdir_info(ssprintf("results"));
  qmkdir_info(ssprintf("results/total_site=%s", show(total_site).c_str()));
  qmkdir_info(ssprintf("results/total_site=%s/lambda=%.10lf",
                       show(total_site).c_str(), lambda));
  if (not obtain_lock(fn + "-lock")) {
    return;
  }
  {
    TIMER_VERBOSE("evolution");
    std::vector<CorrFuncs> cfs;
    // ADJUST ME
    const long max_traj = 2 / 2 * 3;
    const long n_steps = 1;
    //
    ScalarField sf;
    set_scalar_field(sf, total_site);
    for (long traj = 0; traj < max_traj; ++traj) {
      {
        TIMER_VERBOSE("evolution-traj");
        RngField rf;
        set_rng_field(rf, total_site, "seed", traj);
        CorrFuncs cf;
        for (long i = 0; i < n_steps; ++i) {
          sweep_scalar_field(sf, rf, mass_sqr, lambda);
          cf += measure_corr_funcs(sf, cp);
        }
        cf *= 1.0 / (double)n_steps;
        cfs.push_back(cf);
        displayln_info("CHECK: " + ssprintf("%24.13E", cf.phi2));
        display_info(show_results(cfs, total_site, cp, mass_sqr, lambda));
      }
      Timer::autodisplay();
      const double budget = 600.0;
      displayln_info(
          fname +
          ssprintf(": ( get_actual_total_time() + budget ) / get_time_limit() "
                   "= ( %.2lf + %.2lf ) / %.2lf hours.",
                   get_actual_total_time() / 3600.0, budget / 3600.0,
                   get_time_limit() / 3600.0));
      if (budget + get_actual_total_time() > get_time_limit()) {
        displayln_info(fname +
                       ssprintf(": quit because too little time left."));
        break;
      }
    }
    qtouch_info(fn, show_results(cfs, total_site, cp, mass_sqr, lambda));
    release_lock();
  }
  Timer::display();
  // ADJUST ME
  check_time_limit();
  // qquit();
  //
}

inline void compute_several_mass(const Coordinate& total_site,
                                 const CorrParams& cp,
                                 const double mass_sqr_start,
                                 const double mass_sqr_step,
                                 const double lambda)
{
  // ADJUST ME
  for (int i = 0; i < 2; ++i) {
    evolution(total_site, cp, mass_sqr_start + mass_sqr_step * i, lambda);
  }
  //
}

inline void compute(const Coordinate& total_site,
                    const CorrParams& cp)
{
  // ADJUST ME
  // compute_several_mass(total_site, cp, +0.01, 0.01, 0.0);
  // compute_several_mass(total_site, cp, -0.01, 0.01, 0.1);
  // compute_several_mass(total_site, cp, -0.03, 0.01, 0.4);
  // compute_several_mass(total_site, cp, -0.16, 0.01, 1.6);
  // compute_several_mass(total_site, cp, -0.37, 0.01, 4.0);
  // compute_several_mass(total_site, cp, -0.70, 0.02, 8.0);
  // compute_several_mass(total_site, cp, -1.28, 0.02, 16.0);
  // compute_several_mass(total_site, cp, -1.86, 0.02, 24.0);
  // compute_several_mass(total_site, cp, -2.36, 0.02, 32.0);
  if (total_site == Coordinate(4,4,4,256)) {
    compute_several_mass(total_site, cp, -1.86, 0.04, 24.0);
    compute_several_mass(total_site, cp, -2.36, 0.04, 32.0);
  } else if (total_site == Coordinate(8,8,8,512)) {
    compute_several_mass(total_site, cp, -1.69, 0.01, 24.0);
    compute_several_mass(total_site, cp, -2.18, 0.01, 32.0);
  }
  //
}

}  // namespace qlat

int main(int argc, char* argv[])
{
  using namespace qlat;
  std::vector<Coordinate> size_node_list;
  size_node_list.push_back(Coordinate(1, 1, 1, 1));
  size_node_list.push_back(Coordinate(1, 1, 1, 2));
  size_node_list.push_back(Coordinate(1, 1, 1, 4));
  size_node_list.push_back(Coordinate(1, 1, 1, 8));
  size_node_list.push_back(Coordinate(1, 1, 1, 16));
  size_node_list.push_back(Coordinate(1, 1, 1, 32));
  size_node_list.push_back(Coordinate(1, 1, 1, 64));
  size_node_list.push_back(Coordinate(1, 1, 1, 128));
  begin(&argc, &argv, size_node_list);
  setup();
  // ADJUST ME
  {
    evolution(Coordinate(8, 8, 8, 32), CorrParams(2, 4, 1), 0.01, 0.0);
    evolution(Coordinate(8, 8, 8, 32), CorrParams(2, 4, 1), 0.00, 0.1);
    evolution(Coordinate(8, 8, 8, 32), CorrParams(2, 4, 1), 0.01, 0.1);
    evolution(Coordinate(16, 16, 16, 64), CorrParams(2, 4, 1), 0.01, 0.1);
    evolution(Coordinate(16, 16, 16, 64), CorrParams(2, 4, 1), -0.1, 0.5);
  }
  {
    const Coordinate total_site = Coordinate(4, 4, 4, 256);
    const int t1 = 2;
    const int t2 = 4;
    const int dt = 1;
    const CorrParams cp(t1, t2, dt);
    compute(total_site, cp);
  }
  {
    const Coordinate total_site = Coordinate(4, 4, 4, 512);
    const int t1 = 4;
    const int t2 = 8;
    const int dt = 2;
    const CorrParams cp(t1, t2, dt);
    compute(total_site, cp);
  }
  {
    const Coordinate total_site = Coordinate(8, 8, 8, 512);
    const int t1 = 4;
    const int t2 = 8;
    const int dt = 2;
    const CorrParams cp(t1, t2, dt);
    compute(total_site, cp);
  }
  //
  displayln_info("CHECK: finished successfully.");
  Timer::display();
  end();
  return 0;
}
