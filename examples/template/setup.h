#include <qlat/qlat.h>

QLAT_START_NAMESPACE

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

inline void setup(const std::string& job_tag)
{
  qmkdir_info(get_job_path(job_tag));
  qmkdir_info(get_job_path(job_tag) + "/logs");
  switch_monitor_file_info(get_job_path(job_tag) + ssprintf("/logs/%010ld.txt", get_log_idx()));
  Timer::max_function_name_length_shown() = 50;
  Timer::max_call_times_for_always_show_info() = 3;
  Timer::minimum_duration_for_show_stop_info() = 60;
  Timer::minimum_autodisplay_interval() = 365 * 24 * 3600;
  get_lock_expiration_time_limit() = 1.0 * 60.0 * 60.0;
  set_lock_expiration_time_limit();
  get_time_limit() = get_lock_expiration_time_limit();
  get_default_budget() = 15.0 * 60.0;
  dist_write_par_limit() = 128;
  dist_read_par_limit() = 128;
  displayln_info(ssprintf("get_start_time()=%lf", get_start_time()));
  displayln_info(ssprintf("lock_expiration_time=%lf", get_start_time() + get_lock_expiration_time_limit()));
  displayln_info(ssprintf("get_time_limit()=%lf hours", get_time_limit() / 3600.0));
  displayln_info(ssprintf("get_default_budget()=%lf hours", get_default_budget() / 3600.0));
  displayln_info(ssprintf("dist_read_par_limit()=%d", dist_read_par_limit()));
  displayln_info(ssprintf("dist_write_par_limit()=%d", dist_write_par_limit()));
}

// -----------------------------------------------------------------------------------

inline std::vector<int> get_trajs(const std::string& job_tag)
{
  TIMER_VERBOSE("get_trajs");
  std::vector<int> ret;
  if (job_tag == "free-4nt8") {
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
  } else {
    qassert(false);
  }
  return ret;
}

inline Coordinate get_total_site(const std::string& job_tag)
{
  if (job_tag == "free-4nt8") {
    return Coordinate(4,4,4,8);
  } else if (job_tag == "16I-0.01") {
    return Coordinate(16,16,16,32);
  } else if (job_tag == "24I-0.01") {
    return Coordinate(24,24,24,64);
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
  if (job_tag == "free-4nt8") {
    qassert(false);
  } else if (job_tag == "16I-0.01") {
    fn = get_env("HOME") + ssprintf("/qcdarchive/DWF_iwa_nf2p1/16c32/2plus1_16nt32_IWASAKI_b2p13_ls16_M1p8_ms0p04_mu0p01_rhmc_multi_timescale_ukqcd/ckpoint_lat.IEEE64BIG.%d", traj);
    if (does_file_exist_sync_node(fn)) {
      return fn;
    }
    return "";
  } else if (job_tag == "24I-0.01") {
    fn = get_env("HOME") + ssprintf("/qcdarchive/DWF_iwa_nf2p1/24c64/2plus1_24nt64_IWASAKI_b2p13_ls16_M1p8_ms0p04_mu0p005_rhmc_H_R_G/ckpoint_lat.IEEE64BIG.%d", traj);
    if (does_file_exist_sync_node(fn)) {
      return fn;
    }
    fn = get_env("HOME") + ssprintf("/qcdarchive/DWF_iwa_nf2p1/24c64/2plus1_24nt64_IWASAKI_b2p13_ls16_M1p8_ms0p04_mu0p01_quo_hasenbusch_rhmc/ckpoint_lat.IEEE64BIG.%d", traj);
    if (does_file_exist_sync_node(fn)) {
      return fn;
    }
    fn = get_env("HOME") + ssprintf("/qcdarchive/DWF_iwa_nf2p1/24c64/2plus1_24nt64_IWASAKI_b2p13_ls16_M1p8_ms0p04_mu0p01_quo_hasenbusch_rhmc_ukqcd/ckpoint_lat.IEEE64BIG.%d", traj);
    if (does_file_exist_sync_node(fn)) {
      return fn;
    }
    fn = get_env("HOME") + ssprintf("/qcdarchive/DWF_iwa_nf2p1/24c64/2plus1_24nt64_IWASAKI_b2p13_ls16_M1p8_ms0p04_mu0p01_rhmc_H_R_G/ckpoint_lat.IEEE64BIG.%d", traj);
    if (does_file_exist_sync_node(fn)) {
      return fn;
    }
    return "";
  } else {
    qassert(false);
  }
  return "";
}

inline long load_configuration(GaugeField& gf, const std::string& job_tag, const int traj)
{
  TIMER_VERBOSE("load_configuration");
  long file_size = 0;
  const Geometry geo = get_geo(job_tag);
  gf.init(geo);
  if (job_tag == "free-4nt8") {
    set_unit(gf);
    file_size += geo.geon.num_node * get_data_size(gf);
  } else {
    file_size += load_gauge_field(gf, get_config_fn(job_tag, traj));
  }
  qassert(is_matching_geo(geo, gf.geo));
  return file_size;
}

inline std::vector<FermionAction> get_fermion_actions(const std::string& job_tag)
{
  std::vector<FermionAction> fas;
  if (job_tag == "free-4nt8") {
    fas.push_back(FermionAction(0.1, 8, 1.0, 1.0, true, true));
    fas.push_back(FermionAction(0.3, 8, 1.0, 1.0, true, true));
  } else if (job_tag == "24I-0.01" or job_tag == "16I-0.01") {
    fas.push_back(FermionAction(0.01, 16, 1.8));
    fas.push_back(FermionAction(0.04, 16, 1.8));
  } else {
    qassert(false);
  }
  return fas;
}

inline LancArg get_lanc_arg(const std::string& job_tag)
{
  LancArg la;
  if (job_tag == "free-4nt8") {
    qassert(false);
  } else if (job_tag == "16I-0.01") {
    la = LancArg(5.5, 0.20, 200, 200, 110, 100);
    // la = LancArg(5.5, 0.50, 200, 600, 510, 500);
  } else if (job_tag == "24I-0.01") {
    la = LancArg(5.5, 0.18, 200, 1100, 700, 600);
  } else {
    qassert(false);
  }
  return la;
}

inline std::string get_low_modes_path(const std::string& job_tag, const int traj)
{
  if (job_tag == "free-4nt8") {
    qassert(false);
    return "";
  } else {
    qmkdir_info("lancs");
    qmkdir_info(ssprintf("lancs/%s", job_tag.c_str()));
    qmkdir_sync_node(ssprintf("lancs/%s/qcdtraj=%d", job_tag.c_str(), traj));
    return ssprintf("lancs/%s/qcdtraj=%d/huge-data-lanc", job_tag.c_str(), traj);
  }
}

QLAT_END_NAMESPACE
