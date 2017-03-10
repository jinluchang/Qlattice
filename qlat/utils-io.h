// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <qlat/config.h>
#include <qlat/mpi.h>

#include <array>
#include <vector>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdarg>
#include <fstream>

#include <sys/stat.h>
#include <unistd.h>

QLAT_START_NAMESPACE

inline bool truncate(const std::string &evilFile) {
  std::ofstream evil;
  evil.open(evilFile.c_str());
  bool does_exist = evil.good();
  if(does_exist) evil.close();
  return does_exist;
}

inline bool does_file_exist(const std::string& fn) {
  struct stat sb;
  return 0 == stat(fn.c_str(), &sb);
}

inline bool does_file_exist_sync_node(const std::string& fn) {
  long nfile = 0;
  if (0 == get_id_node()) {
    if (does_file_exist(fn)) {
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

inline int ssleep(const double seconds) {
  return usleep((useconds_t)(seconds * 1.0e6));
}

inline int check_dir(const std::string& path, const mode_t mode = default_dir_mode())
{
  TIMER("check_dir");
  int ret = 0;
  while (!does_file_exist(path)) {
    ret = mkdir(path.c_str(), mode);
    ssleep(0.001);
  }
  return ret;
}

inline int qmkdir(const std::string& path, const mode_t mode = default_dir_mode())
{
  TIMER("qmkdir");
  mkdir(path.c_str(), mode);
  return check_dir(path, mode);
}

inline int qmkdir_info(const std::string& path, const mode_t mode = default_dir_mode())
{
  TIMER("qmkdir_info");
  if (0 == get_id_node()) {
    return qmkdir(path, mode);
  } else {
    return 0;
  }
}

inline int qmkdir_sync_node(const std::string& path, const mode_t mode = default_dir_mode())
{
  TIMER("qmkdir_sync_node");
  if (0 == get_id_node()) {
    qmkdir(path, mode);
  }
  sync_node();
  return check_dir(path, mode);
}

inline int mkdir_lock(const std::string& path, const mode_t mode = default_dir_mode())
{
  TIMER("mkdir_lock");
  long ret = 0;
  if (0 == get_id_node()) {
    ret = mkdir(path.c_str(), mode);
  }
  glb_sum(ret);
  return ret;
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

inline int qremove(const std::string& path)
{
  TIMER("qremove");
  return std::remove(path.c_str());
}

inline int qremove_info(const std::string& path)
{
  TIMER("qremove_info");
  if (0 == get_id_node()) {
    return qremove(path);
  } else {
    return 0;
  }
}

inline std::string get_env(const std::string& var_name) {
  const char* value = getenv(var_name.c_str());
  return std::string(value);
}

inline FILE* qopen(const std::string& path, const std::string& mode)
{
  TIMER("qopen");
  return std::fopen(path.c_str(), mode.c_str());
}

inline void qset_line_buf(FILE* f)
{
  TIMER("qset_line_buf");
  std::setvbuf(f, NULL, _IOLBF, 0);
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
    qset_line_buf(f);
    return f;
  } else {
    return NULL;
  }
}

inline int qclose(FILE*& file)
{
  TIMER("qclose");
  if (NULL != file) {
    FILE* tmp_file = file;
    file = NULL;
    return std::fclose(tmp_file);
  }
  return 0;
}

inline int qclose_info(FILE*& file)
{
  TIMER("qclose_info");
  return qclose(file);
}

inline int qrename(const std::string& old_path, const std::string& new_path)
{
  TIMER("qrename");
  return rename(old_path.c_str(), new_path.c_str());
}

inline int qrename_info(const std::string& old_path, const std::string& new_path)
{
  TIMER("qrename_info");
  if (0 == get_id_node()) {
    return qrename(old_path, new_path);
  } else {
    return 0;
  }
}

inline int qtouch(const std::string& path, const std::string& content = "")
{
  TIMER("qtouch");
  FILE* file = qopen(path, "w");
  display(content, file);
  return qclose(file);
}

inline int qtouch_info(const std::string& path, const std::string& content = "")
{
  TIMER("qtouch_info");
  if (0 == get_id_node()) {
    return qtouch(path, content);
  } else {
    return 0;
  }
}

inline std::string qcat(const std::string& path)
{
  TIMER("qcat");
  FILE* fp = qopen(path, "r");
  fseek(fp, 0, SEEK_END);
  const long length = ftell(fp);
  fseek(fp, 0, SEEK_SET);
  std::string ret(length, 0);
  fread(&ret[0], 1, length, fp);
  qclose(fp);
  return ret;
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

inline std::string& qcat_sync_node(const std::string& path)
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

inline bool obtain_lock(const std::string& path)
{
  TIMER_VERBOSE("obtain_lock");
  const std::string path_time = path + "/time.txt";
  displayln_info(ssprintf("%s: Trying to obtain lock '%s'.", fname, path.c_str()));
  qassert(get_lock_location() == "");
  if (0 == mkdir_lock(path)) {
    qtouch_info(path_time, show(get_time()) + "\n");
    get_lock_location() = path;
    displayln_info(ssprintf("%s: Lock obtained '%s'.", fname, path.c_str()));
    return true;
  } else if (does_file_exist_sync_node(path_time)) {
    long ret = 0;
    if (0 == get_id_node()) {
      double time;
      reads(time, qcat_info(path_time));
      if (get_time() - time > get_lock_expiration_time_limit() && 0 == qremove(path_time)) {
        ret = 1;
        qtouch(path_time, show(get_time()) + "\n");
      }
    }
    glb_sum(ret);
    if (ret > 0) {
      displayln_info(ssprintf("%s: Lock obtained '%s' (old lock expired).", fname, path.c_str()));
      return true;
    } else {
      displayln_info(ssprintf("%s: Failed to obtained '%s'.", fname, path.c_str()));
      return false;
    }
  } else {
    displayln_info(ssprintf("%s: Failed to obtained '%s' (no creation time info).", fname, path.c_str()));
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

inline void check_time_limit(const double budget = get_default_budget())
{
  if (budget + get_total_time() > get_time_limit()) {
    release_lock();
    const bool time_out = false;
    assert(time_out);
  }
}

inline std::string qgetline(FILE* fp)
{
  char* lineptr = NULL;
  size_t n = 0;
  if (getline(&lineptr, &n, fp) > 0) {
    std::string ret(lineptr);
    std::free(lineptr);
    return ret;
  } else {
    std::free(lineptr);
    return std::string();
  }
}

inline bool is_space(const char c)
{
  return c == ' ' || c == '\n' || c == '\r' || c == '\t';
}

inline std::vector<std::string> split_line_with_spaces(const std::string& str)
{
  const size_t len = str.length();
  std::vector<std::string> words;
  size_t start = 0;
  size_t stop = 0;
  while (stop < len) {
    while (start < len && is_space(str[start])) {
      start += 1;
    }
    stop = start;
    while (stop < len && !is_space(str[stop])) {
      stop += 1;
    }
    if (stop > start) {
      words.push_back(std::string(str, start, stop - start));
    }
    start = stop;
  }
  return words;
}

inline double read_double(const std::string& str)
{
  double x;
  sscanf(str.c_str(), "%lf", &x);
  return x;
}

inline std::vector<double> read_doubles(const std::string& str)
{
  const std::vector<std::string> strs = split_line_with_spaces(str);
  std::vector<double> ret(strs.size());
  for (size_t i = 0; i < strs.size(); ++i) {
    ret[i] = read_double(strs[i]);
  }
  return ret;
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
        xss[i].clear();
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
  DataTable ret = qload_datatable_par(fp);
  qclose(fp);
  return ret;
}

QLAT_END_NAMESPACE
