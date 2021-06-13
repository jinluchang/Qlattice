#pragma once

#include <qlat/qlat.h>

namespace qlat
{  //

inline long& get_log_idx()
{
  static long idx = -1;
  return idx;
}

inline void setup_log_idx(const std::string& path = "results")
{
  TIMER_VERBOSE("setup_log_idx");
  qmkdir_info(path);
  qmkdir_info(path + "/logs-lock");
  for (long i = 0; true; ++i) {
    const std::string fn = path + ssprintf("/logs-lock/%010d", i);
    if (!does_file_exist_sync_node(fn) && 0 == mkdir_lock(fn)) {
      get_log_idx() = i;
      return;
    }
  }
  qassert(false);
  return;
}

inline std::string get_job_path(const std::string& job_tag)
{
  return "results/" + job_tag;
}

inline std::string get_job_path(const std::string& job_tag, const int traj)
{
  return get_job_path(job_tag) + ssprintf("/results=%d", traj);
}

inline void setup()
{
  Timer::minimum_duration_for_show_info() = 1.0;
  displayln_info(
      ssprintf("get_time_limit()=%lf hours", get_time_limit() / 3600.0));
  displayln_info(ssprintf("get_default_budget()=%lf hours",
                          get_default_budget() / 3600.0));
  displayln_info(ssprintf("dist_read_par_limit()=%d", dist_read_par_limit()));
  displayln_info(ssprintf("dist_write_par_limit()=%d", dist_write_par_limit()));
}

inline void setup(const std::string& job_tag)
{
  qmkdir_info(get_job_path(job_tag));
  qmkdir_info(get_job_path(job_tag) + "/logs");
  switch_monitor_file_info(get_job_path(job_tag) +
                           ssprintf("/logs/%010ld.txt", get_log_idx()));
  setup();
  if (check_status()) {
    qquit("setup(job_tag)");
  }
}

inline void setup(const std::string& job_tag, const int traj)
{
  const std::string job_path = get_job_path(job_tag, traj);
  qmkdir_info("results");
  qmkdir_info(get_job_path(job_tag));
  qmkdir_info(job_path);
  qmkdir_sync_node(job_path + "/logs");
  switch_monitor_file_info(job_path +
                           ssprintf("/logs/%010ld.txt", get_log_idx()));
  TIMER_VERBOSE("setup(job_tag,traj)");
  setup();
  if (check_status()) {
    qquit("setup(job_tag,traj)");
  }
}

// -----------------------------------------------------------------------------------

inline std::vector<int> get_trajs(const std::string& job_tag)
{
  TIMER_VERBOSE("get_trajs");
  std::vector<int> ret;
  if (job_tag.substr(0, 5) == "free-") {
    for (int traj = 1000; traj < 1020; traj += 10) {
      ret.push_back(traj);
    }
  } else if (job_tag.substr(0, 5) == "test-") {
    for (int traj = 1000; traj < 1020; traj += 10) {
      ret.push_back(traj);
    }
  } else if (job_tag == "16I-0.01") {
    for (int traj = 1000; traj <= 4000; traj += 100) {
      ret.push_back(traj);
    }
  } else if (job_tag == "24I-0.01") {
    for (int traj = 2700; traj <= 8500; traj += 100) {
      ret.push_back(traj);
    }
  } else if (job_tag == "48I-0.00078") {
    for (int traj = 500; traj <= 3000; traj += 10) {
      ret.push_back(traj);
    }
  } else if (job_tag == "64I-0.000678") {
    for (int traj = 1000; traj <= 3000; traj += 10) {
      ret.push_back(traj);
    }
  } else if (job_tag == "24D-0.00107") {
    for (int traj = 1000; traj <= 4000; traj += 10) {
      ret.push_back(traj);
    }
  } else if (job_tag == "24D-0.0174") {
    for (int traj = 1000; traj >= 200; traj -= 10) {
      ret.push_back(traj);
    }
  } else if (job_tag == "32D-0.00107") {
    for (int traj = 680; traj <= 2000; traj += 10) {
      ret.push_back(traj);
    }
  } else if (job_tag == "32Dfine-0.0001") {
    for (int traj = 2000; traj >= 200; traj -= 10) {
      ret.push_back(traj);
    }
  } else {
    qassert(false);
  }
  return ret;
}

inline Coordinate get_total_site(const std::string& job_tag)
{
  if (job_tag.substr(0, 5) == "free-" or job_tag.substr(0, 5) == "test-") {
    // eg: job_tag = "free-4nt8"
    long size_l_dir = 0;
    long size_t_dir = 0;
    long cur = 5;
    qassert(parse_long(size_l_dir, cur, job_tag));
    char c = 0;
    qassert(parse_char(c, cur, job_tag) and c == 'n');
    qassert(parse_char(c, cur, job_tag) and c == 't');
    qassert(parse_long(size_t_dir, cur, job_tag));
    return Coordinate(size_l_dir, size_l_dir, size_l_dir, size_t_dir);
  } else if (job_tag == "16I-0.01") {
    return Coordinate(16, 16, 16, 32);
  } else if (job_tag == "24I-0.01" or job_tag == "24D-0.00107" or
             job_tag == "24D-0.0174" or job_tag == "24D" or job_tag == "24DH") {
    return Coordinate(24, 24, 24, 64);
  } else if (job_tag == "48I-0.00078" or job_tag == "48I") {
    return Coordinate(48, 48, 48, 96);
  } else if (job_tag == "64I-0.000678" or job_tag == "64I") {
    return Coordinate(64, 64, 64, 128);
  } else if (job_tag == "32D-0.00107" or job_tag == "32Dfine-0.0001" or
             job_tag == "32D" or job_tag == "32Dfine") {
    return Coordinate(32, 32, 32, 64);
  } else {
    qassert(false);
    return Coordinate();
  }
}

inline Geometry get_geo(const std::string& job_tag)
{
  Geometry geo;
  geo.init(get_total_site(job_tag), 1);
  return geo;
}

inline std::string get_config_fn(const std::string& job_tag, const int traj)
{
  std::string fn("");
  if (job_tag.substr(0, 5) == "free-" or job_tag.substr(0, 5) == "test-") {
    return job_tag + ssprintf("/%d", traj);
  } else if (job_tag == "16I-0.01") {
    fn = get_env("HOME") +
         ssprintf(
             "/qcdarchive/DWF_iwa_nf2p1/16c32/"
             "2plus1_16nt32_IWASAKI_b2p13_ls16_M1p8_ms0p04_mu0p01_rhmc_multi_"
             "timescale_ukqcd/ckpoint_lat.IEEE64BIG.%d",
             traj);
    if (does_file_exist_sync_node(fn)) {
      return fn;
    }
  } else if (job_tag == "24I-0.005") {
    fn = get_env("HOME") + ssprintf(
                               "/qcdarchive/DWF_iwa_nf2p1/24c64/"
                               "2plus1_24nt64_IWASAKI_b2p13_ls16_M1p8_ms0p04_"
                               "mu0p005_rhmc_H_R_G/ckpoint_lat.IEEE64BIG.%d",
                               traj);
    if (does_file_exist_sync_node(fn)) {
      return fn;
    }
  } else if (job_tag == "24I-0.01") {
    fn = get_env("HOME") +
         ssprintf(
             "/qcdarchive/DWF_iwa_nf2p1/24c64/"
             "2plus1_24nt64_IWASAKI_b2p13_ls16_M1p8_ms0p04_mu0p01_quo_"
             "hasenbusch_rhmc/ckpoint_lat.IEEE64BIG.%d",
             traj);
    if (does_file_exist_sync_node(fn)) {
      return fn;
    }
    fn = get_env("HOME") +
         ssprintf(
             "/qcdarchive/DWF_iwa_nf2p1/24c64/"
             "2plus1_24nt64_IWASAKI_b2p13_ls16_M1p8_ms0p04_mu0p01_quo_"
             "hasenbusch_rhmc_ukqcd/ckpoint_lat.IEEE64BIG.%d",
             traj);
    if (does_file_exist_sync_node(fn)) {
      return fn;
    }
    fn = get_env("HOME") + ssprintf(
                               "/qcdarchive/DWF_iwa_nf2p1/24c64/"
                               "2plus1_24nt64_IWASAKI_b2p13_ls16_M1p8_ms0p04_"
                               "mu0p01_rhmc_H_R_G/ckpoint_lat.IEEE64BIG.%d",
                               traj);
    if (does_file_exist_sync_node(fn)) {
      return fn;
    }
  } else if (job_tag == "48I-0.00078") {
    fn =
        get_env("HOME") +
        ssprintf(
            "/qcdarchive-ljin/48nt96-ainv1.73gev-mpi139mev-ls24/ckpoint_lat.%d",
            traj);
    if (does_file_exist_sync_node(fn)) {
      return fn;
    }
  } else if (job_tag == "64I-0.000678") {
    fn = get_env("HOME") +
         ssprintf(
             "/qcdarchive-ljin/64nt128-ainv2.359gev-mpi139mev-ls12-gfixed/"
             "ckpoint_lat.Coulomb.%d",
             traj);
    if (does_file_exist_sync_node(fn)) {
      return fn;
    }
  } else if (job_tag == "24D-0.00107") {
    fn = get_env("HOME") +
         ssprintf(
             "/qcddata-chulwoo/DWF/2+1f/24nt64/IWASAKI+DSDR/b1.633/ls24/M1.8/"
             "ms0.0850/ml0.00107/evol0/configurations/ckpoint_lat.%d",
             traj);
    if (does_file_exist_sync_node(fn)) {
      return fn;
    }
    fn = get_env("HOME") +
         ssprintf(
             "/qcddata-chulwoo/DWF/2+1f/24nt64/IWASAKI+DSDR/b1.633/ls24/M1.8/"
             "ms0.0850/ml0.00107/evol1/configurations/ckpoint_lat.%d",
             traj);
    if (does_file_exist_sync_node(fn)) {
      return fn;
    }
    fn = get_env("HOME") +
         ssprintf(
             "/qcddata-chulwoo/DWF/2+1f/24nt64/IWASAKI+DSDR/b1.633/ls24/M1.8/"
             "ms0.0850/ml0.00107/evol2/configurations/ckpoint_lat.%d",
             traj);
    if (does_file_exist_sync_node(fn)) {
      return fn;
    }
  } else if (job_tag == "24D-0.0174") {
    fn = get_env("HOME") +
         ssprintf(
             "/qcddata-chulwoo/DWF/2+1f/24nt64/IWASAKI+DSDR/b1.633/ls24/M1.8/"
             "ms0.0985/ml0.0174/configurations/ckpoint_lat.%d",
             traj);
    if (does_file_exist_sync_node(fn)) {
      return fn;
    }
  } else if (job_tag == "32D-0.00107") {
    fn = get_env("HOME") +
         ssprintf(
             "/qcddata-chulwoo/DWF/2+1f/32nt64/IWASAKI+DSDR/b1.633/ls24/M1.8/"
             "ms0.0850/ml0.00107/evol0/configurations/ckpoint_lat.%d",
             traj);
    if (does_file_exist_sync_node(fn)) {
      return fn;
    }
  } else if (job_tag == "32Dfine-0.0001") {
    fn = get_env("HOME") +
         ssprintf(
             "/qcddata-chulwoo/DWF/2+1f/32nt64/IWASAKI+DSDR/b1.75/ls32/ms0.045/"
             "mu0.0001/evol0/configurations/ckpoint_lat.%d",
             traj);
    if (does_file_exist_sync_node(fn)) {
      return fn;
    }
  } else {
    qassert(false);
  }
  return "";
}

inline long load_configuration(GaugeField& gf, const std::string& job_tag,
                               const int traj)
{
  TIMER_VERBOSE("load_configuration");
  long file_size = 0;
  const Geometry geo = get_geo(job_tag);
  gf.init(geo);
  if (job_tag.substr(0, 5) == "free-") {
    set_unit(gf);
    file_size += geo.geon.num_node * get_data_size(gf);
  } else if (job_tag.substr(0, 5) == "test-") {
    set_unit(gf);
    const RngState rs =
        RngState().split(job_tag).split(traj).split("load_configuration");
    set_g_rand_color_matrix_field(gf, rs, 1.0);
    for (int i = 0; i < 5; ++i) {
      gf_ape_smear(gf, gf, 0.5);
    }
    file_size += geo.geon.num_node * get_data_size(gf);
  } else {
    file_size += load_gauge_field(gf, get_config_fn(job_tag, traj));
    if (job_tag == "24D-0.00107" or job_tag == "32D-0.00107" or
        job_tag == "48D-0.00107" or job_tag == "32Dfine-0.0001" or
        job_tag == "24D-0.0174" or job_tag == "16I-0.01" or
        job_tag == "24I-0.01") {
      twist_boundary_at_boundary(gf, -0.5, 3);
    }
  }
  qassert(is_matching_geo(geo, gf.geo()));
  return file_size;
}

inline bool& is_partially_quench_for_48I_and_64I()
{
  static bool b = true;
  return b;
}

inline std::vector<FermionAction> get_fermion_actions(
    const std::string& job_tag)
// FermionAction(const double mass_, const int ls_, const double m5_,
//     const double mobius_scale_ = 1.0, const bool is_multiplying_dminus_ =
//     true, bool is_using_zmobius_ = false)
{
  std::vector<FermionAction> fas;
  if (job_tag.substr(0, 5) == "free-") {
    fas.push_back(FermionAction(0.1, 8, 1.0, 1.0, true, true));
    fas.push_back(FermionAction(0.3, 8, 1.0, 1.0, true, true));
  } else if (job_tag.substr(0, 5) == "test-") {
    fas.push_back(FermionAction(0.1, 8, 1.0, 1.0, true, true));
    fas.push_back(FermionAction(0.3, 8, 1.0, 1.0, true, true));
  } else if (job_tag == "24I-0.01" or job_tag == "16I-0.01") {
    fas.push_back(FermionAction(0.01, 16, 1.8, 1.0, true, true));
    fas.push_back(FermionAction(0.04, 16, 1.8, 1.0, true, true));
  } else if (job_tag == "48I-0.00078") {
    if (is_partially_quench_for_48I_and_64I()) {
      fas.push_back(
          FermionAction(0.0006979, 10, 1.8, (1.5 + 0.5) * 2.4, true, true));
      fas.push_back(
          FermionAction(0.03580, 10, 1.8, (1.5 + 0.5) * 2.4, true, true));
      fas.push_back(FermionAction(0.0006979, 24, 1.8, 1.5 + 0.5, true, false));
      fas.push_back(FermionAction(0.03580, 24, 1.8, 1.5 + 0.5, true, false));
    } else {
      fas.push_back(
          FermionAction(0.00078, 10, 1.8, (1.5 + 0.5) * 2.4, true, true));
      fas.push_back(
          FermionAction(0.0362, 10, 1.8, (1.5 + 0.5) * 2.4, true, true));
      fas.push_back(FermionAction(0.00078, 24, 1.8, 1.5 + 0.5, true, false));
      fas.push_back(FermionAction(0.0362, 24, 1.8, 1.5 + 0.5, true, false));
    }
    {
      std::vector<Complex> bs(10, 0);
      bs[0] = 8.4292038368159705e-01;
      bs[1] = 9.2289979238280184e-01;
      bs[2] = 1.1017200769981794e+00;
      bs[3] = 1.4219097980542994e+00;
      bs[4] = 1.9620523417564424e+00;
      bs[5] = 2.8654191667525488e+00;
      bs[6] = 4.4659153528626341e+00;
      bs[7] = 5.5498080139636414e+00;
      bs[8] = Complex(4.9320961582039766e+00, -3.5559998543638791e+00);
      bs[9] = Complex(4.9320961582039766e+00, 3.5559998543638791e+00);
      for (int i = 0; i < (int)bs.size(); i++) {
        fas[0].bs[i] = bs[i];
        fas[0].cs[i] = fas[0].bs[i] - 1.0;
        fas[1].bs[i] = bs[i];
        fas[1].cs[i] = fas[0].bs[i] - 1.0;
      }
    }
  } else if (job_tag == "64I-0.000678") {
    if (is_partially_quench_for_48I_and_64I()) {
      fas.push_back(FermionAction(0.0006203, 12, 1.8, 1.5 + 0.5, true, true));
      fas.push_back(FermionAction(0.02539, 12, 1.8, 1.5 + 0.5, true, false));
    } else {
      fas.push_back(FermionAction(0.000678, 12, 1.8, 1.5 + 0.5, true, true));
      fas.push_back(FermionAction(0.02661, 12, 1.8, 1.5 + 0.5, true, false));
    }
  } else if (job_tag == "24D-0.00107" or job_tag == "32D-0.00107" or
             job_tag == "48D-0.00107" or job_tag == "24D-0.0174") {
    fas.push_back(FermionAction(0.00107, 12, 1.8, (2.5 + 1.5) * 2, true, true));
    {
      std::vector<Complex> omega(12, 0);
      omega[0] = 1.0903256131299373;
      omega[1] = 0.9570283702230611;
      omega[2] = 0.7048886040934104;
      omega[3] = 0.48979921782791747;
      omega[4] = 0.328608311201356;
      omega[5] = 0.21664245377015995;
      omega[6] = 0.14121112711957107;
      omega[7] = 0.0907785101745156;
      omega[8] = Complex(0.05608303440064219, -0.007537158177840385);
      omega[9] = Complex(0.05608303440064219, 0.007537158177840385);
      omega[10] = Complex(0.0365221637144842, -0.03343945161367745);
      omega[11] = Complex(0.0365221637144842, 0.03343945161367745);
      for (int i = 0; i < (int)omega.size(); i++) {
        fas[0].bs[i] = 0.5 * (1.0 / omega[i] + 1.0);
        fas[0].cs[i] = fas[0].bs[i] - 1.0;
      }
    }
    fas.push_back(FermionAction(0.0850, 24, 1.8, 2.5 + 1.5, true, false));
    fas.push_back(FermionAction(0.00107, 24, 1.8, 2.5 + 1.5, true, false));
    if (job_tag == "24D-0.0174") {
      fas[0].mass = 0.0174;
      fas[1].mass = 0.0985;
      fas[2].mass = 0.0174;
    }
  } else if (job_tag == "32Dfine-0.0001") {
    fas.push_back(FermionAction(0.0001, 12, 1.8, 32.0 / 12.0, true, true));
    {
      std::vector<Complex> omega(12, 0.375);
      for (int i = 0; i < (int)omega.size(); i++) {
        fas[0].bs[i] = 0.5 * (1.0 / omega[i] + 1.0);
        fas[0].cs[i] = fas[0].bs[i] - 1.0;
      }
    }
    fas.push_back(FermionAction(0.045, 12, 1.8, 32.0 / 12.0, true, false));
    fas.push_back(FermionAction(0.0001, 32, 1.8, 1.0, true, false));
    fas.push_back(FermionAction(0.045, 32, 1.8, 1.0, true, false));
  } else {
    qassert(false);
  }
  return fas;
}

inline LancArg get_lanc_arg(const std::string& job_tag)
// LancArg(double ch_alpha_, double ch_beta_, long ch_ord_, long n_use_, long
// n_get_, long n_true_get_)
{
  LancArg la;
  if (job_tag.substr(0, 5) == "free-") {
    qassert(false);
  } else if (job_tag == "16I-0.01") {
    la = LancArg(5.5, 0.20, 200, 200, 110, 100);
    // la = LancArg(5.5, 0.50, 200, 600, 510, 500);
  } else if (job_tag == "24I-0.01") {
    la = LancArg(5.5, 0.18, 200, 1100, 700, 600);
  } else if (job_tag == "48I-0.00078") {
    la = LancArg(2.5, 0.022, 200, 2600, 2100, 2000);
  } else if (job_tag == "64I-0.000678") {
    la = LancArg(2.5, 0.026, 200, 2600, 2100, 2000);
  } else if (job_tag == "24D-0.00107") {
    la = LancArg(sqrt(5.5), sqrt(0.000684), 200, 2600, 2100, 2000);
  } else if (job_tag == "24D-0.0174") {
    la = LancArg(sqrt(5.5), 0.20, 200, 800, 550, 500);
  } else if (job_tag == "32D-0.00107") {
    la = LancArg(sqrt(5.5), sqrt(5.860158125869930E-04), 300, 4500, 4200, 4000);
  } else if (job_tag == "32Dfine-0.0001") {
    la = LancArg(sqrt(5.5), sqrt(0.0018), 200, 2600, 2100, 2000);
  } else {
    qassert(false);
  }
  return la;
}

inline std::string make_low_modes_path(const std::string& job_tag,
                                       const int traj)
{
  qmkdir_info("lancs");
  qmkdir_info(ssprintf("lancs/%s", job_tag.c_str()));
  qmkdir_info(ssprintf("lancs/%s/qcdtraj=%d", job_tag.c_str(), traj));
  return ssprintf("lancs/%s/qcdtraj=%d/huge-data-lanc", job_tag.c_str(), traj);
}

inline std::string get_low_modes_path(const std::string& job_tag,
                                      const int traj)
{
  TIMER_VERBOSE("get_low_modes_path");
  if (job_tag.substr(0, 5) == "free-") {
    return "";
  } else {
    return ssprintf("lancs/%s/qcdtraj=%d/huge-data-lanc", job_tag.c_str(),
                    traj);
  }
}

}  // namespace qlat
