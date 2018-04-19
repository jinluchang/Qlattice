#pragma once

// used setup measurement calculation
//
// initialize()
//
// get_project_root_rng_state() = RngState();
//
// Setup setup(path);
//
// csi = make_configurations_info_test(total_site);
//
// load_configuration(gf, csi.infos[k])

#include <qlat/qlat.h>

QLAT_START_NAMESPACE

inline void setup_params()
{
  TIMER_VERBOSE("setup_params");
  Timer::max_function_name_length_shown() = 50;
  Timer::max_call_times_for_always_show_info() = 3;
  Timer::minimum_duration_for_show_stop_info() = 60;
  Timer::minimum_autodisplay_interval() = 365 * 24 * 3600;
  get_lock_expiration_time_limit() = 24.0 * 60.0 * 60.0;
  set_lock_expiration_time_limit();
  get_time_limit() = get_lock_expiration_time_limit();
  get_default_budget() = 0;
  dist_write_par_limit() = 128;
  dist_read_par_limit() = 128;
  displayln_info(ssprintf("get_start_time()=%lf", get_start_time()));
  displayln_info(ssprintf("get_lock_expiration_time_limit()=%lf", get_lock_expiration_time_limit()));
  displayln_info(ssprintf("expiration_time=%lf", get_start_time() + get_lock_expiration_time_limit()));
}

inline long& get_log_idx()
{
  static long idx = -1;
  return idx;
}

inline void setup_log_idx(const std::string& path = ".")
{
  TIMER_VERBOSE("setup_log_idx");
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

inline void initialize(const std::string& path = ".")
{
  qset_line_buf(get_output_file());
  setup_params();
  setup_log_idx(path);
}

inline std::string& get_result_path() {
  static std::string path = "";
  return path;
}

inline RngState& get_project_root_rng_state()
{
  static RngState rs;
  return rs;
}

inline void update_log_rng()
{
  qclose(get_monitor_file());
  if (get_result_path() != "") {
    qassert(get_log_idx() >= 0);
    get_monitor_file() = qopen(get_result_path() + ssprintf("/logs/%010d.txt", get_log_idx()), "a");
    qset_line_buf(get_monitor_file());
  }
  get_global_rng_state() = get_project_root_rng_state().split(get_result_path());
}

struct Setup
{
  std::string pre_name;
  //
  Setup(const std::string& name)
  {
    pre_name = get_result_path();
    get_result_path() = name;
    qmkdir_sync_node(get_result_path());
    qmkdir_sync_node(get_result_path() + "/logs");
    update_log_rng();
  }
  //
  ~Setup()
  {
    get_result_path() = pre_name;
    update_log_rng();
  }
};

struct ConfigurationsInfo;

struct ConfigurationInfo
{
  ConstHandle<ConfigurationsInfo> csi;
  int traj;
  std::string path;
  std::string conf_format;
  std::string low_modes_path;
};

struct ConfigurationsInfo
{
  std::string tag;
  Coordinate total_site;
  std::vector<FermionAction> fas;
  LancArg la;
  std::vector<ConfigurationInfo> infos;
};

inline void load_gauge_field_artificial(GaugeField& gf, const std::string& path)
  // assuming gf already initialized and have correct size;
{
  TIMER_VERBOSE("load_gauge_field_artificial");
  set_g_rand_color_matrix_field(gf, RngState(RngState("load_gauge_field_artificial"), path), 1.0);
  for (int i = 0; i < 10000; ++i) {
    gf_ape_smear(gf, gf, 0.1);
    displayln_info(ssprintf("ape-smear 0.1 %d times", i + 1));
    if (gf_avg_plaq(gf) >= 0.6) {
      displayln_info(ssprintf("total ape-smear 0.1 %d times", i + 1));
      break;
    }
  }
}

inline void load_configuration(GaugeField& gf, const ConfigurationInfo& ci)
{
  Geometry geo;
  geo.init(ci.csi().total_site, 4);
  gf.init(geo);
  if (ci.conf_format == "milc") {
    load_gauge_field_milc(gf, ci.path);
  } else if (ci.conf_format == "cps") {
    load_gauge_field(gf, ci.path);
  } else if (ci.conf_format == "cps3x3") {
    load_gauge_field_cps3x3(gf, ci.path);
  } else if (ci.conf_format == "artificial") {
    load_gauge_field_artificial(gf, ci.path);
  } else if (ci.conf_format == "free") {
    set_unit(gf);
  } else {
    qassert(false);
  }
}

inline ConfigurationsInfo make_configurations_info_test(const Coordinate& total_site, const bool make_cis = true)
{
  ConfigurationsInfo csi;
  csi.total_site = total_site;
  csi.fas.push_back(FermionAction(0.1, 8, 1.8));
  csi.fas.push_back(FermionAction(0.4, 8, 1.8));
  if (make_cis) {
    for (int traj = 0; traj < 3; ++traj) {
      ConfigurationInfo ci;
      ci.csi.init(csi);
      ci.traj = traj;
      ci.path = ssprintf("artificial-traj=%d", traj);
      ci.conf_format = "artificial";
      // qmkdir_info("./results-lancs");
      // qmkdir_info("./results-lancs/qcdtraj=" + show(traj));
      // ci.low_modes_path = "./results-lancs"
      //   "/qcdtraj=" + show(traj) + "/huge-data-lanc";
      csi.infos.push_back(ci);
    }
  }
  return csi;
}

inline ConfigurationsInfo make_configurations_info_free(const Coordinate& total_site, const bool make_cis = true)
{
  ConfigurationsInfo csi;
  csi.tag = "free";
  csi.total_site = total_site;
  csi.fas.push_back(FermionAction(0.10, 16, 1.0));
  csi.fas.push_back(FermionAction(0.20, 16, 1.0));
  if (make_cis) {
    for (int traj = 0; traj < 3; ++traj) {
      ConfigurationInfo ci;
      ci.csi.init(csi);
      ci.traj = traj;
      ci.path = ssprintf("free-traj=%d", traj);
      ci.conf_format = "free";
      // qmkdir_info("./results-lancs");
      // qmkdir_info("./results-lancs/qcdtraj=" + show(traj));
      // ci.low_modes_path = "./results-lancs"
      //   "/qcdtraj=" + show(traj) + "/huge-data-lanc";
      csi.infos.push_back(ci);
    }
  }
  return csi;
}

inline ConfigurationsInfo make_configurations_info_16c32_mu0p01_ms0p04(const bool make_cis = true)
{
  ConfigurationsInfo csi;
  csi.tag = "16I";
  csi.total_site = Coordinate(16, 16, 16, 32);
  csi.fas.push_back(FermionAction(0.01, 16, 1.8));
  csi.fas.push_back(FermionAction(0.04, 16, 1.8));
  // csi.la = LancArg(5.5, 0.20, 200, 200, 110, 100);
  csi.la = LancArg(5.5, 0.50, 200, 600, 510, 500);
  if (make_cis) {
    for (int traj = 1000; traj <= 4000; traj += 100) {
      ConfigurationInfo ci;
      ci.csi.init(csi);
      ci.traj = traj;
      ci.path = get_env("HOME") +
        "/qcdarchive/DWF_iwa_nf2p1/16c32"
        "/2plus1_16nt32_IWASAKI_b2p13_ls16_M1p8_ms0p04_mu0p01_rhmc_multi_timescale_ukqcd"
        "/ckpoint_lat.IEEE64BIG." + show(traj);
      ci.conf_format = "cps";
      ci.low_modes_path = get_env("HOME") +
        ssprintf("/application/Public/Muon-GM2-cc/jobs/16I/lancs/qcdtraj=%d/huge-data-lanc", traj);
      if (does_file_exist_sync_node(ci.path)) {
        csi.infos.push_back(ci);
      }
    }
  }
  return csi;
}

inline ConfigurationsInfo make_configurations_info_16c32_mu0p01_ms0p032(const bool make_cis = true)
{
  ConfigurationsInfo csi;
  csi.tag = "16I_ms0p032";
  csi.total_site = Coordinate(16, 16, 16, 32);
  csi.fas.push_back(FermionAction(0.01, 16, 1.8));
  csi.fas.push_back(FermionAction(0.032, 16, 1.8));
  csi.la = LancArg(5.5, 0.20, 200, 200, 110, 100);
  if (make_cis) {
    for (int traj = 1000; traj <= 11000; traj += 5) {
      ConfigurationInfo ci;
      ci.csi.init(csi);
      ci.traj = traj;
      ci.path = get_env("HOME") +
        "/qcdarchive/DWF_iwa_nf2p1/16c32"
        "/2plus1_16nt32_IWASAKI_b2p13_ls16_M1p8_ms0p032_mu0p01_rhmc_H_R_G"
        "/ckpoint_lat.IEEE64BIG." + show(traj);
      ci.conf_format = "cps";
      ci.low_modes_path = "./results-lancs"
        "/qcdtraj=" + show(traj) + "/huge-data-lanc";
      if (does_file_exist_sync_node(ci.path)) {
        csi.infos.push_back(ci);
      }
    }
  }
  return csi;
}

inline std::string find_conf_24c64_mu0p01_ms0p04(const int traj)
{
  const std::string parent_path = get_env("HOME") + "/qcdarchive/DWF_iwa_nf2p1/24c64";
  const std::string fname = "/ckpoint_lat.IEEE64BIG.";
  std::string path;
  path = parent_path +
    "/2plus1_24nt64_IWASAKI_b2p13_ls16_M1p8_ms0p04_mu0p01_quo_hasenbusch_rhmc"
    + fname + show(traj);
  if (does_file_exist_sync_node(path)) {
    return path;
  }
  path = parent_path +
    "/2plus1_24nt64_IWASAKI_b2p13_ls16_M1p8_ms0p04_mu0p01_quo_hasenbusch_rhmc_ukqcd"
    + fname + show(traj);
  if (does_file_exist_sync_node(path)) {
    return path;
  }
  path = parent_path +
    "/2plus1_24nt64_IWASAKI_b2p13_ls16_M1p8_ms0p04_mu0p01_rhmc_H_R_G"
    + fname + show(traj);
  if (does_file_exist_sync_node(path)) {
    return path;
  }
  path = parent_path +
    "/2plus1_24nt64_IWASAKI_b2p13_ls16_M1p8_ms0p04_mu0p01_rhmc_multi_timescale_ukqcd"
    + fname + show(traj);
  if (does_file_exist_sync_node(path)) {
    return path;
  }
  return "";
}

inline ConfigurationsInfo make_configurations_info_24c64_mu0p01_ms0p04(const bool make_cis = true)
{
  ConfigurationsInfo csi;
  csi.tag = "24I";
  csi.total_site = Coordinate(24, 24, 24, 64);
  csi.fas.push_back(FermionAction(0.01, 16, 1.8));
  csi.fas.push_back(FermionAction(0.04, 16, 1.8));
  csi.la = LancArg(5.5, 0.18, 200, 1100, 700, 600);
  if (make_cis) {
    for (int traj = 1000; traj <= 10000; traj += 5) {
      ConfigurationInfo ci;
      ci.csi.init(csi);
      ci.traj = traj;
      ci.path = find_conf_24c64_mu0p01_ms0p04(traj);
      ci.conf_format = "cps";
      qmkdir_info("./results-lancs");
      qmkdir_info("./results-lancs/qcdtraj=" + show(traj));
      ci.low_modes_path = "./results-lancs"
        "/qcdtraj=" + show(traj) + "/huge-data-lanc";
      if (ci.path != "") {
        csi.infos.push_back(ci);
      }
    }
  }
  return csi;
}

inline std::string find_conf_24c64_mu0p005_ms0p04(const int traj)
{
  const std::string parent_path = get_env("HOME") + "/qcdarchive/DWF_iwa_nf2p1/24c64";
  const std::string fname = "/ckpoint_lat.IEEE64BIG.";
  std::string path;
  path = parent_path +
    "/2plus1_24nt64_IWASAKI_b2p13_ls16_M1p8_ms0p04_mu0p005_rhmc_H_R_G"
    + fname + show(traj);
  if (does_file_exist_sync_node(path)) {
    return path;
  }
  return "";
}

inline ConfigurationsInfo make_configurations_info_24c64_mu0p005_ms0p04(const bool make_cis = true)
{
  ConfigurationsInfo csi;
  csi.tag = "24IL";
  csi.total_site = Coordinate(24, 24, 24, 64);
  csi.fas.push_back(FermionAction(0.005, 16, 1.8));
  csi.fas.push_back(FermionAction(0.04, 16, 1.8));
  csi.la = LancArg(5.5, 0.16, 100, 600, 555, 550);
  if (make_cis) {
    for (int traj = 8545; traj >= 2000 ; traj -= 40) {
      ConfigurationInfo ci;
      ci.csi.init(csi);
      ci.traj = traj;
      ci.path = find_conf_24c64_mu0p005_ms0p04(traj);
      ci.conf_format = "cps";
      qmkdir_info("./results-lancs");
      qmkdir_info("./results-lancs/qcdtraj=" + show(traj));
      ci.low_modes_path = "./results-lancs"
        "/qcdtraj=" + show(traj) + "/huge-data-lanc";
      if (ci.path != "") {
        csi.infos.push_back(ci);
      }
    }
  }
  return csi;
}

inline ConfigurationsInfo make_configurations_info_32c64_dsdr_mu0p001_ms0p045(const bool make_cis = true)
{
  ConfigurationsInfo csi;
  csi.tag = "32ID";
  csi.total_site = Coordinate(32, 32, 32, 64);
  csi.fas.push_back(FermionAction(0.001, 12, 1.8, 32.0/12.0));
  csi.fas.push_back(FermionAction(0.045, 12, 1.8, 32.0/12.0));
  // csi.la = LancArg(15.0, 0.08, 200, 600, 555, 550);
  csi.la = LancArg(15.0, 0.22, 200, 2600, 2100, 2000);
  if (make_cis) {
    for (int traj = 520; traj <= 1280; traj += 40) {
      ConfigurationInfo ci;
      ci.csi.init(csi);
      ci.traj = traj;
      ci.path = get_env("HOME") +
        "/qcdarchive-zbai/32nt64-ainv-1.37gev-mpi171mev-ls32"
        "/ckpoint_lat.IEEE64BIG." + show(traj);
      ci.conf_format = "cps3x3";
        "/qcdtraj=" + show(traj) + "/huge-data-lanc";
      ci.low_modes_path = "/home/ljin/application/Public/Muon-GM2-cc/jobs/32ID/lancs"
        "/qcdtraj=" + show(traj) + "/huge-data-lanc";
      if (ci.path != "") {
        csi.infos.push_back(ci);
      }
    }
  }
  return csi;
}

inline ConfigurationsInfo make_configurations_info_32c64_dsdr_mu0p001_ms0p045_unitary(const bool make_cis = true)
{
  ConfigurationsInfo csi;
  csi.tag = "32ID_unitary";
  csi.total_site = Coordinate(32, 32, 32, 64);
  csi.fas.push_back(FermionAction(0.001, 32, 1.8));
  csi.fas.push_back(FermionAction(0.045, 32, 1.8));
  // csi.la = LancArg(15.0, 0.08, 200, 600, 555, 550);
  // csi.la = LancArg(15.0, 0.22, 200, 2600, 2100, 2000);
  if (make_cis) {
    for (int traj = 520; traj <= 1280; traj += 40) {
      ConfigurationInfo ci;
      ci.csi.init(csi);
      ci.traj = traj;
      ci.path = get_env("HOME") +
        "/qcdarchive-zbai/32nt64-ainv-1.37gev-mpi171mev-ls32"
        "/ckpoint_lat.IEEE64BIG." + show(traj);
      ci.conf_format = "cps3x3";
        "/qcdtraj=" + show(traj) + "/huge-data-lanc";
      qmkdir_info("./results-lancs");
      qmkdir_info("./results-lancs/qcdtraj=" + show(traj));
      ci.low_modes_path = "./results-lancs"
        "/qcdtraj=" + show(traj) + "/huge-data-lanc";
      if (ci.path != "") {
        csi.infos.push_back(ci);
      }
    }
  }
  return csi;
}

inline std::string find_conf_24c64_dsdr_mu0p0017_ms0p0850(const int traj)
{
  const std::string parent_path = get_env("HOME") + "/qcdarchive-jtu/24nt64/IWASAKI+DSDR/b1.633/ls24/M1.8/ms0.0850/ml0.00107";
  const std::string fname = "/ckpoint_lat.";
  std::string path;
  path = parent_path +
    "/evol0/configurations"
    + fname + show(traj);
  if (does_file_exist_sync_node(path)) {
    return path;
  }
  path = parent_path +
    "/evol1/configurations"
    + fname + show(traj);
  if (does_file_exist_sync_node(path)) {
    return path;
  }
  path = parent_path +
    "/evol2/configurations"
    + fname + show(traj);
  if (does_file_exist_sync_node(path)) {
    return path;
  }
  return "";
}

inline ConfigurationsInfo make_configurations_info_24c64_dsdr_mu0p00107_ms0p0850(const bool make_cis = true)
{
  ConfigurationsInfo csi;
  csi.tag = "24D";
  csi.total_site = Coordinate(24, 24, 24, 64);
  csi.fas.push_back(FermionAction(0.00107, 12, 1.8, 4.0));
  csi.fas.push_back(FermionAction(0.0850, 12, 1.8, 4.0));
  csi.la = LancArg(5.5, 0.02, 200, 200, 150, 50);
  for (int i = 0; i < csi.fas.size(); ++i) {
    FermionAction& fa = csi.fas[i];
    fa.is_using_zmobius = true;
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
    qassert(fa.bs.size() == fa.ls);
    qassert(fa.cs.size() == fa.ls);
    qassert(omega.size() == fa.ls);
    for (int i = 0; i < omega.size(); i++) {
      fa.bs[i] = 0.5 * (1.0 / omega[i] + 1.0);
      fa.cs[i] = fa.bs[i] - 1.0;
    }
  }
  if (make_cis) {
    for (int traj = 2280; traj >= 1000 ; traj -= 40) {
      ConfigurationInfo ci;
      ci.csi.init(csi);
      ci.traj = traj;
      ci.path = find_conf_24c64_dsdr_mu0p0017_ms0p0850(traj);
      ci.conf_format = "cps";
      ci.low_modes_path = get_env("HOME") + ssprintf("/application/Public/Muon-GM2-cc/jobs/24D/lancs/qcdtraj=%d/huge-data-clanc", traj);
      if (ci.path != "") {
        csi.infos.push_back(ci);
      }
    }
  }
  return csi;
}

inline ConfigurationsInfo make_configurations_info_24c64_dsdr_mu0p00107_ms0p0850_unitary(const bool make_cis = true)
{
  ConfigurationsInfo csi;
  csi.tag = "24D_unitary";
  csi.total_site = Coordinate(24, 24, 24, 64);
  csi.fas.push_back(FermionAction(0.00107, 24, 1.8, 4.0));
  csi.fas.push_back(FermionAction(0.0850, 24, 1.8, 4.0));
  if (make_cis) {
    for (int traj = 2280; traj >= 1000 ; traj -= 40) {
      ConfigurationInfo ci;
      ci.csi.init(csi);
      ci.traj = traj;
      ci.path = find_conf_24c64_dsdr_mu0p0017_ms0p0850(traj);
      ci.conf_format = "cps";
      if (ci.path != "") {
        csi.infos.push_back(ci);
      }
    }
  }
  return csi;
}

inline ConfigurationsInfo make_configurations_info_milc(const bool make_cis = true)
{
  ConfigurationsInfo csi;
  csi.total_site = Coordinate(24, 24, 24, 64);
  csi.fas.push_back(FermionAction(0.01, 8, 1.8));
  csi.fas.push_back(FermionAction(0.04, 8, 1.8));
  if (make_cis) {
    for (int traj = 500; traj <= 1900; traj += 100) {
      ConfigurationInfo ci;
      ci.csi.init(csi);
      ci.traj = traj;
      ci.path = get_env("HOME") +
        "/qcdarchive-milc/24c64"
        "/2plus1plus1"
        "/l2464f211b600m0102m0509m635a.hyp." + show(traj);
      ci.conf_format = "milc";
      if (does_file_exist_sync_node(ci.path)) {
        csi.infos.push_back(ci);
      }
    }
  }
  return csi;
}

QLAT_END_NAMESPACE
