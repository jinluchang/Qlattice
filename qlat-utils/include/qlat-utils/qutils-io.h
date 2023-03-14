#pragma once

#include <qlat-utils/env.h>
#include <qlat-utils/qutils-vec.h>
#include <qlat-utils/qutils.h>
#include <qlat-utils/timer.h>
#include <qlat-utils/cache.h>
#include <signal.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <fstream>
#include <iostream>

namespace qlat
{  //

std::string dirname(const std::string& fn);

std::string basename(const std::string& fn);

std::string remove_trailing_slashes(const std::string& fn);

int qrename(const std::string& old_path, const std::string& new_path);

int qrename_info(const std::string& old_path, const std::string& new_path);

std::vector<std::string> qls(const std::string& path,
                             const bool is_sort = true);

std::vector<std::string> qls_all(const std::string& path,
                                 const bool is_folder_before_files = false,
                                 const bool is_sort = true);

bool does_file_exist(const std::string& fn);

bool is_directory(const std::string& fn);

bool is_regular_file(const std::string& fn);

int qremove(const std::string& path);

int qremove_all(const std::string& path);

int qmkdir(const std::string& path, const mode_t mode = default_dir_mode());

int qmkdir_p(const std::string& path_, const mode_t mode = default_dir_mode());

int qmkdir_info(const std::string& path,
                const mode_t mode = default_dir_mode());

int qmkdir_p_info(const std::string& path,
                  const mode_t mode = default_dir_mode());

void flush();

// --------------------------

API inline Cache<std::string, bool>& get_is_directory_cache()
// Note: key should end with '/'.
{
  static Cache<std::string, bool> cache("IsDirectoryCache", 1024, 128);
  return cache;
}

inline void clear_is_directory_cache() { get_is_directory_cache().clear(); }

bool is_directory_cache(const std::string& dir_);

bool is_regular_file_cache(const std::string& fn);

bool does_file_exist_cache(const std::string& fn);

// --------------------------

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

inline int qfclose(FILE*& file)
{
  TIMER("qfclose");
  if (NULL != file) {
    FILE* tmp_file = file;
    file = NULL;
    return std::fclose(tmp_file);
  }
  return 0;
}

inline void switch_monitor_file(const std::string& path)
{
  qfclose(get_monitor_file());
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

API inline int& is_sigterm_received()
{
  static int n = 0;
  return n;
}

API inline double& get_last_sigint_time()
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

inline void switch_monitor_file_info(const std::string& path)
{
  if (0 == get_id_node()) {
    switch_monitor_file(path);
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

inline int qfclose_info(FILE*& file)
{
  TIMER("qfclose_info");
  return qfclose(file);
}

}  // namespace qlat
