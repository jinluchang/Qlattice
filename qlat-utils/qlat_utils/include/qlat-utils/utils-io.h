#pragma once

#include <qlat-utils/env.h>
#include <qlat-utils/utils-vec.h>
#include <qlat-utils/utils.h>
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

std::vector<std::string> all_dirname_vec(const std::string& fn);

std::string remove_trailing_slashes(const std::string& fn);

std::vector<std::string> qls(const std::string& path,
                             const bool is_sort = true);

std::vector<std::string> qls_all(const std::string& path,
                                 const bool is_folder_before_files = false,
                                 const bool is_sort = true);

bool does_file_exist(const std::string& fn);

bool is_directory(const std::string& fn);

bool is_regular_file(const std::string& fn);

Int qmkdir(const std::string& path, const mode_t mode = default_dir_mode());

Int qmkdir_p(const std::string& path_, const mode_t mode = default_dir_mode());

Int qrename(const std::string& old_path, const std::string& new_path);

Int qremove(const std::string& path);

Int qremove_all(const std::string& path);

bool check_dir(const std::string& path, const mode_t mode = default_dir_mode());

void flush();

// --------------------------

API inline Cache<std::string, bool>& get_is_directory_cache()
// Note: key should end with '/'.
{
  static Cache<std::string, bool> cache("IsDirectoryCache", 1024, 128);
  return cache;
}

inline void clear_is_directory_cache() { get_is_directory_cache().clear(); }

void remove_entry_directory_cache(const std::string& dir_);

void add_entry_directory_cache(const std::string& dir_, bool is_directory);

bool is_directory_cache(const std::string& dir_);

bool is_regular_file_cache(const std::string& fn);

bool does_file_exist_cache(const std::string& fn);

// --------------------------

Int qtruncate(const std::string& path, const Long offset=0);

inline Int ssleep(const RealD seconds)
{
  return usleep((useconds_t)(seconds * 1.0e6));
}

// --------------------------

Int qmkdir_info(const std::string& path,
                const mode_t mode = default_dir_mode());

Int qmkdir_p_info(const std::string& path,
                  const mode_t mode = default_dir_mode());

Int qrename_info(const std::string& old_path, const std::string& new_path);

Int qremove_info(const std::string& path);

Int qremove_all_info(const std::string& path);

// --------------------------

std::vector<std::string> qls_sync_node(const std::string& path,
                                       const bool is_sort = true);

std::vector<std::string> qls_all_sync_node(
    const std::string& path, const bool is_folder_before_files = false,
    const bool is_sort = true);

bool does_file_exist_sync_node(const std::string& fn);

bool is_directory_sync_node(const std::string& fn);

bool is_regular_file_sync_node(const std::string& fn);

bool does_file_exist_cache_sync_node(const std::string& fn);

bool is_directory_cache_sync_node(const std::string& fn);

bool is_regular_file_cache_sync_node(const std::string& fn);

Int qmkdir_sync_node(const std::string& path,
                     const mode_t mode = default_dir_mode());

Int qmkdir_p_sync_node(const std::string& path,
                       const mode_t mode = default_dir_mode());

Int qremove_sync_node(const std::string& path);

Int qremove_all_sync_node(const std::string& path);

// --------------------------

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

inline FILE* qopen(const std::string& path, const std::string& mode)
{
  TIMER("qopen");
  return std::fopen(path.c_str(), mode.c_str());
}

inline Int qfclose(FILE*& file)
{
  TIMER("qfclose");
  if (NULL != file) {
    FILE* tmp_file = file;
    file = NULL;
    return std::fclose(tmp_file);
  }
  return 0;
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

inline Int qfclose_info(FILE*& file)
{
  TIMER("qfclose_info");
  return qfclose(file);
}

template <class M>
Long qwrite_data(const Vector<M>& v, FILE* fp)
{
  TIMER_FLOPS("qwrite_data");
  timer.flops += v.data_size();
  return sizeof(M) * std::fwrite((void*)v.p, sizeof(M), v.n, fp);
}

template <class M>
Long qread_data(const Vector<M>& v, FILE* fp)
{
  TIMER_FLOPS("qread_data");
  timer.flops += v.data_size();
  return sizeof(M) * std::fread((void*)v.p, sizeof(M), v.n, fp);
}

inline std::string qgetline(FILE* fp)
{
  char* lineptr = NULL;
  size_t n = 0;
  const Long size = getline(&lineptr, &n, fp);
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

API inline Int& is_sigterm_received()
{
  static Int n = 0;
  return n;
}

API inline RealD& get_last_sigint_time()
{
  static RealD time = 0.0;
  return time;
}

inline void qhandler_sig(const Int signum)
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
    const RealD time = get_total_time();
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

inline Int install_qhandle_sig()
{
  TIMER_VERBOSE("install_qhandle_sig");
  struct sigaction act;
  act.sa_handler = qhandler_sig;
  return sigaction(SIGINT, &act, NULL) + sigaction(SIGTERM, &act, NULL);
}

}  // namespace qlat
