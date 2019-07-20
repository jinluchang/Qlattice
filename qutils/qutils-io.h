#pragma once

#include "qutils.h"
#include "timer.h"

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

}  // namespace qlat
