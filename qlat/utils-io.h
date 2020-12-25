// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <dirent.h>
#include <qlat/config.h>
#include <qlat/mpi.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

namespace qlat
{  //

inline void set_time_limit_auto();

inline void set_default_budget_auto();

inline void release_lock();

inline double& get_actual_start_time()
// not affected by Timer::reset()
{
  static double time = get_start_time();
  return time;
}

inline double get_actual_total_time()
{
  return get_time() - get_actual_start_time();
}

inline int ssleep(const double seconds)
{
  return usleep((useconds_t)(seconds * 1.0e6));
}

inline double& get_time_limit()
// qlat parameter
{
  static double limit = 0.0;
  if (0.0 == limit) {
    TIMER_VERBOSE("get_time_limit");
    limit = 12.0 * 3600.0;
    set_time_limit_auto();
  }
  return limit;
}

inline double& get_default_budget()
// qlat parameter
{
  static double budget = 0.0;
  if (0.0 == budget) {
    TIMER_VERBOSE("get_default_budget");
    budget = 15.0 * 60.0;
    set_default_budget_auto();
  }
  return budget;
}

inline void set_default_budget_auto()
{
  TIMER_VERBOSE("set_default_budget_auto");
  std::string stime = get_env("Q_BUDGET");
  if (stime == "") {
    stime = get_env("q_budget");
  }
  if (stime != "") {
    double budget = 0.0;
    reads(budget, stime);
    get_default_budget() = budget;
    displayln_info(fname + ssprintf(": get_default_budget() = %.2lf hours.",
                                    get_default_budget() / 3600.0));
  }
}

inline void set_time_limit_auto()
{
  TIMER_VERBOSE("set_time_limit_auto");
  std::string setime = get_env("Q_END_TIME");
  if (setime == "") {
    setime = get_env("q_end_time");
  }
  std::string stime = get_env("Q_TIME_LIMIT");
  if (stime == "") {
    stime = get_env("q_time_limit");
  }
  const std::string ss = get_env("COBALT_STARTTIME");
  const std::string se = get_env("COBALT_ENDTIME");
  if (setime != "") {
    double etime = 0.0;
    reads(etime, setime);
    get_time_limit() = etime - get_actual_start_time();
    displayln_info(fname + ssprintf(": via Q_END_TIME."));
  } else if (stime != "") {
    double time = 0.0;
    reads(time, stime);
    get_time_limit() = time;
    displayln_info(fname + ssprintf(": via Q_TIME_LIMIT."));
  } else if (ss != "" and se != "") {
    double start_time = 0.0;
    double end_time = 0.0;
    reads(start_time, ss);
    reads(end_time, se);
    get_time_limit() = end_time - get_actual_start_time();
    displayln_info(fname + ssprintf(": via COBALT_ENDTIME."));
    displayln_info(fname + ssprintf(": job total time = %.2lf hours.",
                                    (end_time - start_time) / 3600.0));
    displayln_info(fname +
                   ssprintf(": job init time = %.2lf hours.",
                            (get_actual_start_time() - start_time) / 3600.0));
  }
  displayln_info(fname + ssprintf(": get_time_limit() = %.2lf hours.",
                                  get_time_limit() / 3600.0));
}

inline double get_remaining_time()
{
  return get_time_limit() - get_actual_total_time();
}

inline void qquit(const std::string& msg)
{
  release_lock();
  Timer::display();
  displayln_info("qquit: " + msg);
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
  if (budget + get_actual_total_time() > get_time_limit()) {
    qquit("because too little time left.");
  }
}

inline void check_time_limit(bool timer_display,
                             const double budget = get_default_budget())
// obsolete
{
  displayln_info(
      "WARNING: do not use this function. "
      "check_time_limit(timer_display,budget)");
  check_time_limit(budget, timer_display);
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
  set_time_limit_auto();
}

inline void check_sigint()
{
  if (is_sigint_received() > 0) {
    qquit("because sigint received.");
  }
}

inline bool does_file_exist(const std::string& fn)
{
  struct stat sb;
  return 0 == stat(fn.c_str(), &sb);
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

inline bool is_directory(const std::string& fn)
{
  struct stat sb;
  if (0 != stat(fn.c_str(), &sb)) {
    return false;
  }
  return S_ISDIR(sb.st_mode);
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

inline void check_stop()
{
  if (does_file_exist_sync_node("stop.txt")) {
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
  if (is_sigint_received() > 0) {
    displayln_info(fname + ssprintf(": sigint received."));
    return true;
  }
  if (does_file_exist_sync_node("stop.txt")) {
    displayln_info(fname + ssprintf(": File 'stop.txt' detected."));
    return true;
  }
  return false;
}

inline bool truncate(const std::string& evilFile)
{
  std::ofstream evil;
  evil.open(evilFile.c_str());
  bool does_exist = evil.good();
  if (does_exist) evil.close();
  return does_exist;
}

inline mode_t& default_dir_mode()
// qlat parameter
{
  static mode_t mode = 0775;
  return mode;
}

inline int check_dir(const std::string& path,
                     const mode_t mode = default_dir_mode())
{
  TIMER("check_dir");
  int ret = 0;
  while (!does_file_exist(path)) {
    ret = mkdir(path.c_str(), mode);
    ssleep(0.001);
  }
  return ret;
}

inline int qmkdir(const std::string& path,
                  const mode_t mode = default_dir_mode())
{
  TIMER("qmkdir");
  mkdir(path.c_str(), mode);
  return check_dir(path, mode);
}

inline int qmkdir_info(const std::string& path,
                       const mode_t mode = default_dir_mode())
{
  TIMER("qmkdir_info");
  if (0 == get_id_node()) {
    return qmkdir(path, mode);
  } else {
    return 0;
  }
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

inline std::string remove_trailing_slashes(const std::string& fn)
{
  long cur = fn.size() - 1;
  while (cur > 0 and fn[cur] == '/') {
    cur -= 1;
  }
  return std::string(fn, 0, cur + 1);
}

inline std::vector<std::string> qls_aux(const std::string& path)
{
  std::vector<std::string> contents;
  DIR* dir = opendir(path.c_str());
  if (dir == NULL) {
    return contents;
  }
  struct dirent* d;
  while ((d = readdir(dir)) != NULL) {
    if (!strcmp(d->d_name, ".") || !strcmp(d->d_name, "..")) {
      continue;
    }
    contents.push_back(path + "/" + d->d_name);
  }
  closedir(dir);
  return contents;
}

inline std::vector<std::string> qls(const std::string& path)
{
  return qls_aux(remove_trailing_slashes(path));
}

inline int qremove(const std::string& path)
{
  displayln(ssprintf("qremove: '%s'", path.c_str()));
  return std::remove(path.c_str());
}

inline int qremove_all_aux(const std::string& path)
{
  if (not is_directory(path)) {
    return qremove(path);
  } else {
    int ret = 0;
    const std::vector<std::string> paths = qls_aux(path);
    for (long i = 0; i < (long)paths.size(); ++i) {
      ret += qremove_all_aux(paths[i]);
    }
    return ret + qremove(path);
  }
}

inline int qremove_all(const std::string& path)
{
  return qremove_all_aux(remove_trailing_slashes(path));
}

inline int qremove_info(const std::string& path)
{
  TIMER_VERBOSE("qremove_info");
  if (0 == get_id_node()) {
    return qremove(path);
  } else {
    return 0;
  }
}

inline int qremove_all_info(const std::string& path)
{
  TIMER_VERBOSE("qremove_all_info");
  if (0 == get_id_node()) {
    return qremove_all(path);
  } else {
    return 0;
  }
}

inline std::string show_file_crc32(const std::pair<std::string, crc32_t>& fcrc)
{
  return ssprintf("%08X  fn='%s'", fcrc.second, fcrc.first.c_str());
}

inline std::string show_files_crc32(
    const std::vector<std::pair<std::string, crc32_t> >& fcrcs)
{
  std::ostringstream out;
  for (long i = 0; i < (long)fcrcs.size(); ++i) {
    out << ssprintf("%5ld ", i) << show_file_crc32(fcrcs[i]) << std::endl;
  }
  return out.str();
}

inline std::pair<std::string, crc32_t> check_file_crc32(const std::string& fn)
{
  TIMER_VERBOSE("check_file_crc32");
  std::pair<std::string, crc32_t> p;
  p.first = fn;
  p.second = compute_crc32(fn);
  displayln_info(show_file_crc32(p));
  return p;
}

inline void check_all_files_crc32_aux(
    std::vector<std::pair<std::string, crc32_t> >& acc, const std::string& path)
{
  if (not is_directory(path)) {
    acc.push_back(check_file_crc32(path));
  } else {
    const std::vector<std::string> paths = qls_aux(path);
    for (long i = 0; i < (long)paths.size(); ++i) {
      check_all_files_crc32_aux(acc, paths[i]);
    }
  }
}

inline std::vector<std::pair<std::string, crc32_t> > check_all_files_crc32(
    const std::string& path)
{
  TIMER_VERBOSE("check_all_files_crc32");
  std::vector<std::pair<std::string, crc32_t> > ret;
  check_all_files_crc32_aux(ret, remove_trailing_slashes(path));
  return ret;
}

inline void check_all_files_crc32_info(const std::string& path)
{
  TIMER_VERBOSE("check_all_files_crc32_info");
  if (0 == get_id_node()) {
    displayln(fname + ssprintf(": start checking path='%s'", path.c_str()));
    std::vector<std::pair<std::string, crc32_t> > fcrcs;
    fcrcs = check_all_files_crc32(path);
    displayln(fname + ssprintf(": summary for path='%s'", path.c_str()));
    display(show_files_crc32(fcrcs));
  }
}

inline void qset_fully_buf(FILE* f)
{
  TIMER("qset_fully_buf");
  std::setvbuf(f, NULL, _IOFBF, BUFSIZ);
}

inline FILE* qopen_info(const std::string& path, const std::string& mode)
{
  TIMER("qopen_info");
  if (0 == get_id_node()) {
    FILE* f = qopen(path, mode);
    qassert(f != NULL);
    qset_line_buf(f);
    return f;
  } else {
    return NULL;
  }
}

inline int qclose_info(FILE*& file)
{
  TIMER("qclose_info");
  return qclose(file);
}

inline int qrename_info(const std::string& old_path,
                        const std::string& new_path)
{
  TIMER("qrename_info");
  if (0 == get_id_node()) {
    return qrename(old_path, new_path);
  } else {
    return 0;
  }
}

inline int qtouch_info(const std::string& path)
{
  TIMER("qtouch_info");
  if (0 == get_id_node()) {
    return qtouch(path);
  } else {
    return 0;
  }
}

inline int qtouch_info(const std::string& path, const std::string& content)
{
  TIMER("qtouch_info");
  if (0 == get_id_node()) {
    return qtouch(path, content);
  } else {
    return 0;
  }
}

inline int qappend_info(const std::string& path, const std::string& content)
{
  TIMER("qappend_info");
  if (0 == get_id_node()) {
    return qappend(path, content);
  } else {
    return 0;
  }
}

inline std::string qcat_info(const std::string& path)
{
  TIMER("qcat_info");
  if (0 == get_id_node()) {
    return qcat(path);
  } else {
    return std::string();
  }
}

inline std::string qcat_sync_node(const std::string& path)
{
  TIMER("qcat_sync_node");
  std::string ret;
  long length = 0;
  if (0 == get_id_node()) {
    ret = qcat(path);
    length = ret.length();
  }
  glb_sum(length);
  ret.resize(length, 0);
  bcast(get_data(ret));
  return ret;
}

inline void switch_monitor_file_info(const std::string& path)
{
  if (0 == get_id_node()) {
    switch_monitor_file(path);
  }
}

inline std::string& get_lock_location()
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

typedef std::vector<std::vector<double> > DataTable;

inline DataTable qload_datatable(FILE* fp)
{
  TIMER("qload_datatable(fp)");
  DataTable ret;
  while (!feof(fp)) {
    const std::string line = qgetline(fp);
    if (line.length() > 0 && line[0] != '#') {
      const std::vector<double> xs = read_doubles(line);
      if (xs.size() > 0) {
        ret.push_back(xs);
      }
    }
  }
  return ret;
}

inline DataTable qload_datatable_par(FILE* fp)
{
  TIMER("qload_datatable(fp)");
  const size_t line_buf_size = 1024;
  DataTable ret;
  std::vector<std::string> lines;
  DataTable xss;
  while (!feof(fp)) {
    lines.clear();
    for (size_t i = 0; i < line_buf_size; ++i) {
      lines.push_back(qgetline(fp));
      if (feof(fp)) {
        break;
      }
    }
    xss.resize(lines.size());
#pragma omp parallel for
    for (size_t i = 0; i < lines.size(); ++i) {
      const std::string& line = lines[i];
      if (line.length() > 0 && line[0] != '#') {
        xss[i] = read_doubles(line);
      } else {
        clear(xss[i]);
      }
    }
    for (size_t i = 0; i < xss.size(); ++i) {
      if (xss[i].size() > 0) {
        ret.push_back(xss[i]);
      }
    }
  }
  return ret;
}

inline DataTable qload_datatable(const std::string& path)
{
  TIMER("qload_datatable(path)");
  if (!does_file_exist(path)) {
    return DataTable();
  }
  FILE* fp = qopen(path, "r");
  qassert(fp != NULL);
  DataTable ret = qload_datatable(fp);
  qclose(fp);
  return ret;
}

inline DataTable qload_datatable_par(const std::string& path)
{
  TIMER("qload_datatable(path)");
  if (!does_file_exist(path)) {
    return DataTable();
  }
  FILE* fp = qopen(path, "r");
  qassert(fp != NULL);
  DataTable ret = qload_datatable_par(fp);
  qclose(fp);
  return ret;
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

inline void lat_data_save_info(const std::string& path, const LatData& ld)
{
  TIMER("lat_data_save_info");
  if (get_id_node() == 0) {
    ld.save(path);
  }
}

}  // namespace qlat
