#pragma once

#include <qutils/qutils-vec.h>
#include <qutils/qutils.h>
#include <qutils/timer.h>
#include <signal.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <dirent.h>
#include <fstream>
#include <iostream>

namespace qlat
{  //

inline bool does_file_exist(const std::string& fn)
{
  struct stat sb;
  return 0 == stat(fn.c_str(), &sb);
}

inline bool is_directory(const std::string& fn)
{
  struct stat sb;
  if (0 != stat(fn.c_str(), &sb)) {
    return false;
  }
  return S_ISDIR(sb.st_mode);
}

inline bool is_regular_file(const std::string& fn)
{
  struct stat sb;
  if (0 != stat(fn.c_str(), &sb)) {
    return false;
  }
  return S_ISREG(sb.st_mode);
}

inline bool qtruncate(const std::string& evilFile)
{
  std::ofstream evil;
  evil.open(evilFile.c_str());
  bool does_exist = evil.good();
  if (does_exist) {
    evil.close();
  }
  return does_exist;
}

inline bool qtruncate(const std::string& path, const long offset)
// return true if successful.
{
  const int ret = truncate(path.c_str(), offset);
  return ret == 0;
}

inline mode_t& default_dir_mode()
// qlat parameter
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

inline std::string remove_trailing_slashes(const std::string& fn)
{
  long cur = fn.size() - 1;
  while (cur > 0 and fn[cur] == '/') {
    cur -= 1;
  }
  return std::string(fn, 0, cur + 1);
}

inline std::string dirname(const std::string& fn)
// try to follow libgen.h version see man 3 dirname
{
  long cur = fn.size() - 1;
  // remove trailing '/'
  while (cur > 0 and fn[cur] == '/') {
    cur -= 1;
  }
  if (cur < 0) {
    return ".";
  } else if (cur == 0) {
    if (fn[cur] == '/') {
      return "/";
    } else {
      return ".";
    }
  } else {
    // remove last component
    while (cur >= 0 and fn[cur] != '/') {
      cur -= 1;
    }
    if (cur < 0) {
      return ".";
    } else {
      // remove trailing '/'
      while (cur > 0 and fn[cur] == '/') {
        cur -= 1;
      }
      return std::string(fn, 0, cur + 1);
    }
  }
  qassert(false);
  return std::string();
}

inline std::string basename(const std::string& fn)
// try to follow libgen.h version see man 3 basename
{
  long cur = fn.size() - 1;
  // remove trailing '/'
  while (cur > 0 and fn[cur] == '/') {
    cur -= 1;
  }
  if (cur < 0) {
    return "";
  } else if (cur == 0) {
    if (fn[cur] == '/') {
      return "/";
    } else {
      return std::string(fn, 0, cur + 1);
    }
  } else {
    const long pos_stop = cur + 1;
    // skip last component
    while (cur >= 0 and fn[cur] != '/') {
      cur -= 1;
    }
    return std::string(fn, cur + 1, pos_stop);
  }
  qassert(false);
  return std::string();
}

inline int qmkdir(const std::string& path,
                  const mode_t mode = default_dir_mode())
{
  TIMER("qmkdir");
  mkdir(path.c_str(), mode);
  return check_dir(path, mode);
}

inline int qmkdir_p(const std::string& path_,
                    const mode_t mode = default_dir_mode())
// return 0 if successful
// may try repeatedly in case of failure.
{
  TIMER("qmkdir_p");
  std::string path = remove_trailing_slashes(path_);
  if (is_directory(path)) {
    return 0;
  }
  std::vector<std::string> paths;
  int ret = 0;
  while (true) {
    if (0 == mkdir(path.c_str(), mode)) {
      break;
    } else {
      paths.push_back(path);
      path = dirname(path);
      if (does_file_exist(path)) {
        qwarn(fname + ssprintf(": '%s' failed.", path_.c_str()));
        return 1;
      }
    }
  }
  for (long i = paths.size() - 1; i >= 0; i -= 1) {
    if (not(0 == mkdir(paths[i].c_str(), mode))) {
      qwarn(fname + ssprintf(": '%s' failed.", path_.c_str()));
      return 2;
    }
  }
  return check_dir(path_, mode);
}

inline std::vector<std::string> qls_aux(const std::string& path,
                                        const bool is_sort = true)
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
  if (is_sort) {
    std::sort(contents.begin(), contents.end());
  }
  return contents;
}

inline std::vector<std::string> qls(const std::string& path,
                                    const bool is_sort = true)
{
  return qls_aux(remove_trailing_slashes(path), is_sort);
}

inline std::vector<std::string> qls_all_aux(
    const std::string& path, const bool is_folder_before_files = false,
    const bool is_sort = true)
// list all files and folder in path (not including it self)
{
  std::vector<std::string> all_contents;
  if (not is_directory(path)) {
    return all_contents;
  }
  const std::vector<std::string> contents = qls_aux(path, is_sort);
  for (long i = 0; i < (long)contents.size(); ++i) {
    const std::string& path_i = contents[i];
    if (not is_directory(path_i)) {
      all_contents.push_back(path_i);
    } else {
      if (is_folder_before_files) {
        all_contents.push_back(path_i);
        vector_append(all_contents,
                      qls_all_aux(path_i, is_folder_before_files));
      } else {
        // default behavior
        vector_append(all_contents,
                      qls_all_aux(path_i, is_folder_before_files));
        all_contents.push_back(path_i);
      }
    }
  }
  return all_contents;
}

inline std::vector<std::string> qls_all(
    const std::string& path, const bool is_folder_before_files = false,
    const bool is_sort = true)
// list files before its folder
{
  return qls_all_aux(remove_trailing_slashes(path), is_folder_before_files,
                     is_sort);
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

inline int qrename(const std::string& old_path, const std::string& new_path)
{
  TIMER("qrename");
  if (verbose_level() > 0) {
    displayln(
        ssprintf("qrename: '%s' '%s'", old_path.c_str(), new_path.c_str()));
  }
  return rename(old_path.c_str(), new_path.c_str());
}

inline int qtouch(const std::string& path)
{
  TIMER("qtouch");
  FILE* file = qopen(path, "w");
  qassert(file != NULL);
  return qclose(file);
}

inline int qtouch(const std::string& path, const std::string& content)
{
  TIMER("qtouch");
  FILE* file = qopen(path + ".partial", "w");
  qassert(file != NULL);
  display(content, file);
  qclose(file);
  return qrename(path + ".partial", path);
}

inline int qappend(const std::string& path, const std::string& content)
{
  TIMER("qappend");
  FILE* file = qopen(path, "a");
  qassert(file != NULL);
  display(content, file);
  qclose(file);
  return qrename(path, path);
}

inline std::string qcat(const std::string& path)
{
  TIMER("qcat");
  FILE* fp = qopen(path, "r");
  if (fp == NULL) {
    return "";
  }
  fseek(fp, 0, SEEK_END);
  const long length = ftell(fp);
  fseek(fp, 0, SEEK_SET);
  std::string ret(length, 0);
  const long length_actual = fread(&ret[0], 1, length, fp);
  qassert(length == length_actual);
  qclose(fp);
  return ret;
}

inline void switch_monitor_file(const std::string& path)
{
  qclose(get_monitor_file());
  get_monitor_file() = qopen(path, "a");
  qset_line_buf(get_monitor_file());
}

inline std::string qgetline(FILE* fp)
{
  char* lineptr = NULL;
  size_t n = 0;
  const long size = getline(&lineptr, &n, fp);
  if (size > 0) {
    std::string ret(lineptr, size);
    std::free(lineptr);
    return ret;
  } else {
    std::free(lineptr);
    return std::string();
  }
}

inline std::vector<std::string> qgetlines(FILE* fp)
{
  std::vector<std::string> ret;
  while (!feof(fp)) {
    ret.push_back(qgetline(fp));
  }
  return ret;
}

inline std::vector<std::string> qgetlines(const std::string& fn)
{
  FILE* fp = qopen(fn, "r");
  qassert(fp != NULL);
  std::vector<std::string> lines = qgetlines(fp);
  qclose(fp);
  return lines;
}

inline int& is_sigterm_received()
{
  static int n = 0;
  return n;
}

inline double& get_last_sigint_time()
{
  static double time = 0.0;
  return time;
}

inline void qhandler_sig(const int signum)
{
  if (signum == SIGTERM) {
    is_sigterm_received() += 1;
    displayln(
        ssprintf("qhandler_sig: sigterm triggered, current count is: %d / 10.",
                 is_sigterm_received()));
    Timer::display();
    Timer::display_stack();
    ssleep(3.0);
    if (is_sigterm_received() >= 10) {
      qassert(false);
    }
  } else if (signum == SIGINT) {
    displayln(ssprintf("qhandler_sig: sigint triggered."));
    Timer::display();
    Timer::display_stack();
    const double time = get_total_time();
    if (time - get_last_sigint_time() <= 3.0) {
      displayln(ssprintf(
          "qhandler_sig: sigint triggered interval = %.2f <= 3.0. Quit.",
          time - get_last_sigint_time()));
      qassert(false);
    } else {
      get_last_sigint_time() = time;
    }
  } else {
    displayln(ssprintf("qhandler_sig: cannot handle this signal: %d.", signum));
    Timer::display();
    Timer::display_stack();
  }
}

inline int install_qhandle_sig()
{
  TIMER_VERBOSE("install_qhandle_sig");
  struct sigaction act;
  act.sa_handler = qhandler_sig;
  return sigaction(SIGINT, &act, NULL) + sigaction(SIGTERM, &act, NULL);
}

template <class M>
long qwrite_data(const Vector<M>& v, FILE* fp)
{
  TIMER_FLOPS("qwrite_data");
  timer.flops += v.data_size();
  return sizeof(M) * std::fwrite((void*)v.p, sizeof(M), v.n, fp);
}

template <class M>
long qread_data(const Vector<M>& v, FILE* fp)
{
  TIMER_FLOPS("qread_data");
  timer.flops += v.data_size();
  return sizeof(M) * std::fread((void*)v.p, sizeof(M), v.n, fp);
}

inline crc32_t compute_crc32(const std::string& path)
{
  TIMER_VERBOSE_FLOPS("compute_crc32");
  const size_t chunk_size = 16 * 1024 * 1024;
  std::vector<char> data(chunk_size);
  crc32_t crc = 0;
  FILE* fp = qopen(path, "r");
  qassert(fp != NULL);
  while (true) {
    const long size = qread_data(get_data(data), fp);
    timer.flops += size;
    if (size == 0) {
      break;
    }
    crc = crc32_par(crc, Vector<char>(data.data(), size));
  }
  qclose(fp);
  return crc;
}

}  // namespace qlat
