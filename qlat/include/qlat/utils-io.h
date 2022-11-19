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

inline void close_all_all_shuffled_fields_writer();

inline void release_lock();

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

inline void qquit(const std::string& msg)
// everything needed for gracefully quit and then quit.
{
  clear_all_caches();
  close_all_all_shuffled_fields_writer();
  release_lock();
  Timer::display();
  displayln_info("qquit: " + msg);
  Timer::display_stack();
  ssleep(1.0);
  end();
  ssleep(1.0);
  exit(0);
}

inline void check_time_limit(const double budget = get_default_budget())
{
  TIMER_VERBOSE("check_time_limit");
  displayln_info(
      fname +
      ssprintf(": ( get_actual_total_time() + budget ) / get_time_limit() "
               "= ( %.2lf + %.2lf ) / %.2lf hours.",
               get_actual_total_time() / 3600.0, budget / 3600.0,
               get_time_limit() / 3600.0));
  bool b = budget + get_actual_total_time() > get_time_limit();
  bcast(get_data_one_elem(b));
  if (b) {
    qquit("because too little time left.");
  }
}

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

inline bool does_file_exist_sync_node(const std::string& fn)
{
  long nfile = 0;
  if (0 == get_id_node()) {
    if (does_file_exist(fn)) {
      nfile = 1;
    }
  }
  glb_sum(nfile);
  return 0 != nfile;
}

inline bool does_regular_file_exist_qar_sync_node(const std::string& fn)
{
  long nfile = 0;
  if (0 == get_id_node()) {
    if (does_regular_file_exist_qar(fn)) {
      nfile = 1;
    }
  }
  glb_sum(nfile);
  return 0 != nfile;
}

inline bool does_file_exist_qar_sync_node(const std::string& fn)
{
  long nfile = 0;
  if (0 == get_id_node()) {
    if (does_file_exist_qar(fn)) {
      nfile = 1;
    }
  }
  glb_sum(nfile);
  return 0 != nfile;
}

inline bool is_directory_sync_node(const std::string& fn)
{
  long nfile = 0;
  if (0 == get_id_node()) {
    if (is_directory(fn)) {
      nfile = 1;
    }
  }
  glb_sum(nfile);
  return 0 != nfile;
}

inline bool is_regular_file_sync_node(const std::string& fn)
{
  long nfile = 0;
  if (0 == get_id_node()) {
    if (is_regular_file(fn)) {
      nfile = 1;
    }
  }
  glb_sum(nfile);
  return 0 != nfile;
}

inline void check_stop(const std::string& fn = "stop.txt")
{
  if (does_file_exist_sync_node(fn)) {
    qquit("File 'stop.txt' detected.");
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

inline int qmkdir_sync_node(const std::string& path,
                            const mode_t mode = default_dir_mode())
{
  TIMER("qmkdir_sync_node");
  if (0 == get_id_node()) {
    qmkdir(path, mode);
  }
  sync_node();
  return check_dir(path, mode);
}

inline int mkdir_lock(const std::string& path,
                      const mode_t mode = default_dir_mode())
{
  TIMER("mkdir_lock");
  long ret = 0;
  if (0 == get_id_node()) {
    ret = mkdir(path.c_str(), mode);
  }
  glb_sum(ret);
  return ret;
}

inline int mkdir_lock_all_node(const std::string& path,
                               const mode_t mode = default_dir_mode())
{
  TIMER("mkdir_lock_all_node");
  return mkdir(path.c_str(), mode);
}

inline int rmdir_lock(const std::string& path)
{
  TIMER("rmdir_lock");
  long ret = 0;
  if (0 == get_id_node()) {
    ret = rmdir(path.c_str());
  }
  glb_sum(ret);
  return ret;
}

inline int rmdir_lock_all_node(const std::string& path)
{
  TIMER("rmdir_lock_all_node");
  return rmdir(path.c_str());
}

inline std::vector<std::string> qls_sync_node(const std::string& path)
{
  std::vector<std::string> ret;
  if (0 == get_id_node()) {
    ret = qls(path);
  }
  bcast(ret);
  return ret;
}

inline std::vector<std::string> qls_all_sync_node(
    const std::string& path, const bool is_folder_before_files = false)
{
  std::vector<std::string> ret;
  if (0 == get_id_node()) {
    ret = qls_all(path, is_folder_before_files);
  }
  bcast(ret);
  return ret;
}

inline int qremove_info(const std::string& path)
{
  TIMER_VERBOSE("qremove_info");
  long ret = 0;
  if (0 == get_id_node()) {
    ret = qremove(path);
  }
  glb_sum(ret);
  return ret;
}

inline int qremove_all_info(const std::string& path)
{
  TIMER_VERBOSE("qremove_all_info");
  long ret = 0;
  if (0 == get_id_node()) {
    ret = qremove_all(path);
  }
  glb_sum(ret);
  return ret;
}

inline int qar_create_info(const std::string& path_qar,
                           const std::string& path_folder_,
                           const bool is_remove_folder_after = false)
{
  long ret = 0;
  if (0 == get_id_node()) {
    ret = qar_create(path_qar, path_folder_, is_remove_folder_after);
  }
  glb_sum(ret);
  return ret;
}

inline int qar_extract_info(const std::string& path_qar,
                            const std::string& path_folder_,
                            const bool is_remove_qar_after = false)
{
  long ret = 0;
  if (0 == get_id_node()) {
    ret = qar_extract(path_qar, path_folder_, is_remove_qar_after);
  }
  glb_sum(ret);
  return ret;
}

inline int qcopy_file_info(const std::string& path_src, const std::string& path_dst)
{
  long ret = 0;
  if (0 == get_id_node()) {
    ret = qcopy_file(path_src, path_dst);
  }
  glb_sum(ret);
  return ret;
}

inline std::string qcat_sync_node(const std::string& path)
{
  TIMER("qcat_sync_node");
  std::string ret;
  if (0 == get_id_node()) {
    ret = qcat(path);
  }
  bcast(ret);
  return ret;
}

API inline std::string& get_lock_location()
{
  static std::string path;
  return path;
}

inline bool obtain_lock(const std::string& path)
{
  TIMER_VERBOSE("obtain_lock");
  const std::string path_time = path + "/time.txt";
  const double expiration_time = get_actual_start_time() + get_time_limit();
  displayln_info(fname +
                 ssprintf(": Trying to obtain lock '%s'.", path.c_str()));
  qassert(get_lock_location() == "");
  if (0 == mkdir_lock(path)) {
    qtouch_info(path_time, show(expiration_time) + "\n");
    get_lock_location() = path;
    displayln_info(fname + ssprintf(": Lock obtained '%s'.", path.c_str()));
    return true;
  } else if (does_file_exist_sync_node(path_time)) {
    long ret = 0;
    if (0 == get_id_node()) {
      double time;
      reads(time, qcat(path_time));
      if (get_time() - time > 0.0 && 0 == qremove(path_time)) {
        ret = 1;
        qtouch(path_time, show(expiration_time) + "\n");
      }
    }
    glb_sum(ret);
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

inline void release_lock()
{
  TIMER_VERBOSE("release_lock");
  std::string& path = get_lock_location();
  const std::string path_time = path + "/time.txt";
  displayln_info(fname + ssprintf(": Release lock '%s'", path.c_str()));
  if (path != "") {
    qremove_info(path_time);
    rmdir_lock(path);
    path = "";
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

inline DataTable qload_datatable_sync_node(const std::string& path, const bool is_par = false)
{
  TIMER_VERBOSE("qload_datatable_sync_node");
  DataTable dt;
  if (0 == get_id_node()) {
    dt = qload_datatable(path, is_par);
  }
  bcast(dt);
  return dt;
}

inline LatData lat_data_load_info(const std::string& path)
{
  TIMER("lat_data_load_info");
  LatData ld;
  if (get_id_node() == 0) {
    ld.load(path);
  }
  bcast(ld);
  return ld;
}

}  // namespace qlat
