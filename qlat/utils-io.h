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
#include <stdio.h>

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

inline std::string get_env(const std::string& var_name) {
  const char* value = getenv(var_name.c_str());
  return std::string(value);
}

inline FILE* qopen(const std::string& path, const std::string& mode)
{
  TIMER("qopen");
  return fopen(path.c_str(), mode.c_str());
}

inline void qset_line_buf(FILE* f)
{
  TIMER("qset_line_buf");
  std::setvbuf(f, NULL, _IOLBF, BUFSIZ);
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
    return fclose(tmp_file);
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
  return reads(x, str);
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

inline std::vector<std::vector<double> > qload_datatable(FILE* fp)
{
  TIMER("qload_datatable(fp)");
  const size_t line_buf_size = 128;
  std::vector<std::vector<double> > ret;
  std::vector<std::string> lines;
  std::vector<std::vector<double> > xss;
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

inline std::vector<std::vector<double> > qload_datatable(const std::string& path)
{
  TIMER("qload_datatable(path)");
  FILE* fp = qopen(path, "r");
  std::vector<std::vector<double> > ret = qload_datatable(fp);
  qclose(fp);
  return ret;
}

QLAT_END_NAMESPACE
