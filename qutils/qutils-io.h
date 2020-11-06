#pragma once

#include "qutils.h"
#include "qutils-vec.h"
#include "timer.h"

#include <signal.h>
#include <stdlib.h>

namespace qlat
{  //

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
  if (getline(&lineptr, &n, fp) > 0) {
    std::string ret(lineptr);
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

inline int& is_sigint_received()
{
  static int n = 0;
  return n;
}

inline void qhandler_sigint(const int signum)
{
  if (signum == SIGINT or signum == SIGTERM) {
    is_sigint_received() += 1;
    displayln(ssprintf("qhandler_sigint: triggered, current count is: %d / 10.",
                       is_sigint_received()));
    sleep(3.0);
    if (is_sigint_received() >= 10) {
      qassert(false);
    }
  } else {
    displayln(
        ssprintf("qhandler_sigint: cannot handle this signal: %d.", signum));
  }
}

inline int install_qhandle_sigint()
{
  TIMER_VERBOSE("install_qhandle_sigint");
  struct sigaction act;
  act.sa_handler = qhandler_sigint;
  return sigaction(SIGINT, &act, NULL) + sigaction(SIGTERM, &act, NULL);
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
