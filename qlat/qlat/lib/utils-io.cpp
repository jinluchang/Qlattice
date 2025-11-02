#include <qlat/utils-io.h>
#include <qlat/fields-io.h>

namespace qlat
{  //

bool obtain_lock(const std::string& path)
{
  TIMER_VERBOSE("obtain_lock");
  const std::string path_time = path + "/time.txt";
  const double expiration_time = get_actual_start_time() + get_time_limit();
  displayln_info(fname +
                 ssprintf(": Trying to obtain lock '%s'.", path.c_str()));
  qassert(get_lock_location() == "");
  const std::string path_dir = dirname(path);
  qmkdir_p_sync_node(path_dir);
  if (0 == mkdir_lock(path)) {
    qtouch_info(path_time, show(expiration_time) + "\n");
    get_lock_location() = path;
    displayln_info(fname + ssprintf(": Lock obtained '%s'.", path.c_str()));
    return true;
  } else if (does_file_exist_sync_node(path_time)) {
    Long ret = 0;
    if (0 == get_id_node()) {
      double time;
      QFile qfile = qfopen(path_time, QFileMode::Read);
      if (not qfile.null()) {
        reads(time, qcat(qfile));
        qfclose(qfile);
        if (get_time() - time > 0.0) {
          ret = 1;
        }
      }
    }
    glb_sum(ret);
    if (ret > 0) {
      displayln_info(
          fname +
          ssprintf(": Lock expired '%s'.", path.c_str()));
      return obtain_lock(path + "-");
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

void release_lock()
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

 Int  mkdir_lock(const std::string& path, const mode_t mode)
{
  TIMER("mkdir_lock");
  return qmkdir_sync_node(path, mode);
}

 Int  rmdir_lock(const std::string& path)
{
  TIMER("rmdir_lock");
  return qremove_sync_node(path);
}

void qquit(const std::string& msg)
// everything needed for gracefully quit and then quit.
{
  clear_all_caches();
  close_all_shuffled_fields_writer();
  close_all_shuffled_fields_reader();
  close_all_qar_file();
  release_lock();
  Timer::display();
  displayln_info("qquit: " + msg);
  Timer::display_stack();
  ssleep(1.0);
  end();
  ssleep(1.0);
  exit(0);
}

void check_time_limit(const double budget)
{
  TIMER_VERBOSE("check_time_limit");
  displayln_info(
      fname +
      ssprintf(": ( get_actual_total_time() + budget ) / get_time_limit() "
               "= ( %.2lf + %.2lf ) / %.2lf hours.",
               get_actual_total_time() / 3600.0, budget / 3600.0,
               get_time_limit() / 3600.0));
  double time_deficit = budget + get_actual_total_time() - get_time_limit();
  bcast(time_deficit);
  if (time_deficit > 0) {
    qquit("because too little time left.");
  }
}

void check_stop(const std::string& fn)
{
  if (does_file_exist_sync_node(fn)) {
    qquit("File 'stop.txt' detected.");
  }
}

void check_sigterm()
{
  if (is_sigterm_received() > 0) {
    qquit("because sigterm received.");
  }
}

bool check_status()
{
  TIMER_VERBOSE("check_status");
  displayln_info(fname + ssprintf(": ( get_actual_total_time() + "
                                  "get_time_budget() ) / get_time_limit() "
                                  "= ( %.2lf + %.2lf ) / %.2lf hours.",
                                  get_actual_total_time() / 3600.0,
                                  get_time_budget() / 3600.0,
                                  get_time_limit() / 3600.0));
  if (get_time_budget() + get_actual_total_time() > get_time_limit()) {
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

bool obtain_lock_all_node(const std::string& path)
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
    Long ret = 0;
    double time;
    QFile qfile = qfopen(path_time, QFileMode::Read);
    if (not qfile.null()) {
      reads(time, qcat(qfile));
      qfclose(qfile);
      if (get_time() - time > 0.0) {
        ret = 1;
      }
    }
    if (ret > 0) {
      displayln_info(
          fname +
          ssprintf(": Lock expired '%s'.", path.c_str()));
      return obtain_lock_all_node(path + "-");
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

void release_lock_all_node()
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

 Int  mkdir_lock_all_node(const std::string& path, const mode_t mode)
{
  TIMER("mkdir_lock_all_node");
  return qmkdir(path, mode);
}

 Int  rmdir_lock_all_node(const std::string& path)
{
  TIMER("rmdir_lock_all_node");
  return qremove(path);
}

// double& get_lock_expiration_time_limit()
// // obsolete
// {
//   displayln_info(
//       "WARNING: do not use this function. get_lock_expiration_time_limit");
//   return get_time_limit();
// }
//
// void set_lock_expiration_time_limit()
// // obsolete
// {
//   TIMER_VERBOSE("set_lock_expiration_time_limit");
//   displayln_info(
//       "WARNING: do not use this function. set_lock_expiration_time_limit");
// }

}  // namespace qlat
