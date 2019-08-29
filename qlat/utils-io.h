// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <qlat/config.h>
#include <qlat/mpi.h>

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

QLAT_START_NAMESPACE

inline bool truncate(const std::string& evilFile)
{
  std::ofstream evil;
  evil.open(evilFile.c_str());
  bool does_exist = evil.good();
  if (does_exist) evil.close();
  return does_exist;
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

inline mode_t& default_dir_mode()
{
  static mode_t mode = 0775;
  return mode;
}

inline int ssleep(const double seconds)
{
  return usleep((useconds_t)(seconds * 1.0e6));
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

inline std::string get_env(const std::string& var_name)
{
  const char* value = getenv(var_name.c_str());
  if (value == NULL) {
    return std::string();
  } else {
    return std::string(value);
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

inline double& get_lock_expiration_time_limit()
{
  static double limit = 7.0 * 24 * 3600;
  return limit;
}

inline void set_lock_expiration_time_limit()
{
  TIMER_VERBOSE("set_lock_expiration_time_limit");
  const std::string ss = get_env("COBALT_STARTTIME");
  const std::string se = get_env("COBALT_ENDTIME");
  if (ss != "" and se != "") {
    double start_time, end_time;
    reads(start_time, ss);
    reads(end_time, se);
    get_lock_expiration_time_limit() = end_time - get_start_time();
    displayln_info(fname + ssprintf(": job total time = %.2lf hours.",
                                    (end_time - start_time) / 3600.0));
    displayln_info(fname + ssprintf(": job init time = %.2lf hours.",
                                    (get_start_time() - start_time) / 3600.0));
    displayln_info(fname +
                   ssprintf(": get_lock_expiration_time_limit() = %.2lf hours.",
                            get_lock_expiration_time_limit() / 3600.0));
  } else {
    displayln_info(
        fname +
        ssprintf(
            ": get_lock_expiration_time_limit() = %.2lf hours. (NOT CHANGED)",
            get_lock_expiration_time_limit() / 3600.0));
  }
}

inline bool obtain_lock(const std::string& path)
{
  TIMER_VERBOSE("obtain_lock");
  const std::string path_time = path + "/time.txt";
  const double expiration_time =
      get_start_time() + get_lock_expiration_time_limit();
  displayln_info(
      ssprintf("%s: Trying to obtain lock '%s'.", fname, path.c_str()));
  qassert(get_lock_location() == "");
  if (0 == mkdir_lock(path)) {
    qtouch_info(path_time, show(expiration_time) + "\n");
    get_lock_location() = path;
    displayln_info(ssprintf("%s: Lock obtained '%s'.", fname, path.c_str()));
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
      displayln_info(ssprintf("%s: Lock obtained '%s' (old lock expired).",
                              fname, path.c_str()));
      return true;
    } else {
      displayln_info(
          ssprintf("%s: Failed to obtained '%s'.", fname, path.c_str()));
      return false;
    }
  } else {
    displayln_info(
        ssprintf("%s: Failed to obtained '%s' (no creation time info).", fname,
                 path.c_str()));
    return false;
  }
}

inline bool obtain_lock_all_node(const std::string& path)
{
  TIMER_VERBOSE("obtain_lock_all_node");
  const std::string path_time = path + "/time.txt";
  const double expiration_time =
      get_start_time() + get_lock_expiration_time_limit();
  displayln_info(
      ssprintf("%s: Trying to obtain lock '%s'.", fname, path.c_str()));
  qassert(get_lock_location() == "");
  if (0 == mkdir_lock_all_node(path)) {
    qtouch(path_time, show(expiration_time) + "\n");
    get_lock_location() = path;
    displayln_info(ssprintf("%s: Lock obtained '%s'.", fname, path.c_str()));
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
      displayln_info(ssprintf("%s: Lock obtained '%s' (old lock expired).",
                              fname, path.c_str()));
      return true;
    } else {
      displayln_info(
          ssprintf("%s: Failed to obtained '%s'.", fname, path.c_str()));
      return false;
    }
  } else {
    displayln_info(
        ssprintf("%s: Failed to obtained '%s' (no creation time info).", fname,
                 path.c_str()));
    return false;
  }
}

inline void release_lock()
{
  TIMER_VERBOSE("release_lock");
  std::string& path = get_lock_location();
  const std::string path_time = path + "/time.txt";
  displayln_info(ssprintf("%s: Release lock '%s'", fname, path.c_str()));
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
  displayln_info(ssprintf("%s: Release lock '%s'", fname, path.c_str()));
  if (path != "") {
    qremove(path_time);
    rmdir_lock_all_node(path);
    path = "";
  }
}

inline double& get_time_limit()
{
  static double limit = 7.0 * 24 * 3600;
  return limit;
}

inline double& get_default_budget()
{
  static double budget = 600;
  return budget;
}

inline void check_time_limit(bool timer_display = false,
                             const double budget = get_default_budget())
{
  TIMER_VERBOSE("check_time_limit");
  if (timer_display) {
    Timer::display("check_time_limit");
  }
  if (budget + get_total_time() > get_time_limit()) {
    release_lock();
    const bool time_out = false;
    assert(time_out);
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

QLAT_END_NAMESPACE
