// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <qlat/setup.h>
#include <qlat/mpi.h>
#include <sys/types.h>
#include <unistd.h>

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace qlat
{  //

inline double get_time_limit_default()
{
  const double time_limit_default = 12.0 * 3600.0;
  if (get_env("q_end_time") == "") {
    return get_env_double_default("q_time_limit", time_limit_default);
  } else {
    return get_env_double_default(
               "q_end_time", get_actual_start_time() + time_limit_default) -
           get_actual_start_time();
  }
}

API inline double& get_time_limit()
// qlat parameter
{
  static double limit = get_time_limit_default();
  return limit;
}

API inline double& get_default_budget()
// qlat parameter
{
  static double budget = get_env_double_default("q_budget", 15.0 * 60.0);
  return budget;
}

inline double get_remaining_time()
{
  return get_time_limit() - get_actual_total_time();
}

API inline std::string& get_lock_location()
{
  static std::string path;
  return path;
}

// -------------------

bool obtain_lock(const std::string& path);

void release_lock();

void close_all_shuffled_fields_writer();

void qquit(const std::string& msg);

void check_time_limit(const double budget = get_default_budget());

void check_stop(const std::string& fn = "stop.txt");

bool does_file_exist_sync_node(const std::string& fn);

bool does_regular_file_exist_qar_sync_node(const std::string& fn);

bool does_file_exist_qar_sync_node(const std::string& fn);

bool is_directory_sync_node(const std::string& fn);

bool is_regular_file_sync_node(const std::string& fn);

int qmkdir_sync_node(const std::string& path,
                     const mode_t mode = default_dir_mode());

int mkdir_lock(const std::string& path, const mode_t mode = default_dir_mode());

int mkdir_lock_all_node(const std::string& path, const mode_t mode = default_dir_mode());

int rmdir_lock(const std::string& path);

int rmdir_lock_all_node(const std::string& path);

std::vector<std::string> qls_sync_node(const std::string& path);

std::vector<std::string> qls_all_sync_node(
    const std::string& path, const bool is_folder_before_files = false);

int qremove_info(const std::string& path);

int qremove_all_info(const std::string& path);

int qar_create_info(const std::string& path_qar,
                    const std::string& path_folder_,
                    const bool is_remove_folder_after = false);

int qar_extract_info(const std::string& path_qar,
                     const std::string& path_folder_,
                     const bool is_remove_qar_after = false);

int qcopy_file_info(const std::string& path_src, const std::string& path_dst);

std::string qcat_sync_node(const std::string& path);

DataTable qload_datatable_sync_node(const std::string& path, const bool is_par = false);

LatData lat_data_load_info(const std::string& path);

// ------------------------------

inline double& get_lock_expiration_time_limit()
// obsolete
{
  displayln_info(
      "WARNING: do not use this function. get_lock_expiration_time_limit");
  return get_time_limit();
}

inline void set_lock_expiration_time_limit()
// obsolete
{
  TIMER_VERBOSE("set_lock_expiration_time_limit");
  displayln_info(
      "WARNING: do not use this function. set_lock_expiration_time_limit");
}

inline void check_sigterm()
{
  if (is_sigterm_received() > 0) {
    qquit("because sigterm received.");
  }
}

inline bool check_status()
{
  TIMER_VERBOSE("check_status");
  displayln_info(fname + ssprintf(": ( get_actual_total_time() + "
                                  "get_default_budget() ) / get_time_limit() "
                                  "= ( %.2lf + %.2lf ) / %.2lf hours.",
                                  get_actual_total_time() / 3600.0,
                                  get_default_budget() / 3600.0,
                                  get_time_limit() / 3600.0));
  if (get_default_budget() + get_actual_total_time() > get_time_limit()) {
    displayln_info(fname + ssprintf(": too little time left."));
    return true;
  }
  if (is_sigterm_received() > 0) {
    displayln_info(fname + ssprintf(": sigterm received."));
    return true;
  }
  if (does_file_exist_sync_node("stop.txt")) {
    displayln_info(fname + ssprintf(": File 'stop.txt' detected."));
    return true;
  }
  return false;
}

inline bool obtain_lock_all_node(const std::string& path)
{
  TIMER_VERBOSE("obtain_lock_all_node");
  const std::string path_time = path + "/time.txt";
  const double expiration_time = get_actual_start_time() + get_time_limit();
  displayln_info(fname +
                 ssprintf(": Trying to obtain lock '%s'.", path.c_str()));
  qassert(get_lock_location() == "");
  if (0 == mkdir_lock_all_node(path)) {
    qtouch(path_time, show(expiration_time) + "\n");
    get_lock_location() = path;
    displayln_info(fname + ssprintf(": Lock obtained '%s'.", path.c_str()));
    return true;
  } else if (does_file_exist(path_time)) {
    long ret = 0;
    double time;
    reads(time, qcat(path_time));
    if (get_time() - time > 0.0 && 0 == qremove(path_time)) {
      ret = 1;
      qtouch(path_time, show(expiration_time) + "\n");
    }
    if (ret > 0) {
      get_lock_location() = path;
      displayln_info(
          fname +
          ssprintf(": Lock obtained '%s' (old lock expired).", path.c_str()));
      return true;
    } else {
      displayln_info(fname +
                     ssprintf(": Failed to obtain '%s'.", path.c_str()));
      return false;
    }
  } else {
    displayln_info(fname +
                   ssprintf(": Failed to obtain '%s' (no creation time info).",
                            path.c_str()));
    return false;
  }
}

inline void release_lock_all_node()
{
  TIMER_VERBOSE("release_lock_all_node");
  std::string& path = get_lock_location();
  const std::string path_time = path + "/time.txt";
  displayln_info(fname + ssprintf(": Release lock '%s'", path.c_str()));
  if (path != "") {
    qremove(path_time);
    rmdir_lock_all_node(path);
    path = "";
  }
}


}  // namespace qlat
