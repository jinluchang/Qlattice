// vim: set ts=2 sw=2 expandtab:

// Copyright (c) 2022 Luchang Jin
// All rights reserved.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <cstdarg>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <sstream>
#include <string>
#include <vector>

#define PSTRX(x) #x

#define PSTR(x) PSTRX(x)

#define API __attribute__((visibility("default")))

namespace qlat
{  //

inline std::string vssprintf(const char* fmt, va_list args)
{
  char* cstr;
  int ret = vasprintf(&cstr, fmt, args);
  assert(ret >= 0);
  const std::string str = std::string(cstr);
  std::free(cstr);
  return str;
}

inline std::string ssprintf(const char* fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  return vssprintf(fmt, args);
}

inline std::string show() { return ""; }

inline std::string show(const int& x) { return ssprintf("%d", x); }

inline std::string show(const unsigned int& x) { return ssprintf("%u", x); }

inline std::string show(const long& x) { return ssprintf("%ld", x); }

inline std::string show(const unsigned long& x) { return ssprintf("%lu", x); }

inline std::string show(const double& x) { return ssprintf("%24.17E", x); }

inline std::string show(const bool& x) { return x ? "true" : "false"; }

inline std::string show(const std::string& x)
{
  std::ostringstream out;
  out << x;
  return out.str();
}

inline std::string show_crc32(const uint32_t& x) { return ssprintf("%08X", x); }

template <class T>
std::string show_list(const std::vector<T>& vec)
{
  std::ostringstream out;
  for (long i = 0; i < (long)vec.size(); ++i) {
    out << ssprintf("%5ld: ", i) << show(vec[i]) << std::endl;
  }
  return out.str();
}

template <class T>
std::string shows(const T& x)
{
  std::ostringstream out;
  out << x;
  return out.str();
}

template <class T>
T& reads(T& x, const std::string& str)
{
  std::istringstream in(str);
  in >> x;
  return x;
}

inline long read_long(const std::string& str)
{
  long ret = 0;
  reads(ret, str);
  return ret;
}

inline double read_double(const std::string& str)
{
  double ret = 0.0;
  reads(ret, str);
  return ret;
}

inline uint32_t read_crc32(const std::string& s)
{
  uint32_t crc32;
  std::sscanf(s.c_str(), "%X", &crc32);
  return crc32;
}

inline std::string remove_trailing_newline(const std::string& s)
{
  long cur = s.size() - 1;
  while (cur >= 0 and (s[cur] == '\n' or s[cur] == '\r')) {
    cur -= 1;
  }
  return std::string(s, 0, cur + 1);
}

inline std::vector<std::string> split_into_lines(const std::string& str)
{
  const size_t len = str.length();
  std::vector<std::string> lines;
  size_t start = 0;
  size_t stop = 0;
  while (stop < len) {
    while (start < len && str[start] == '\n') {
      start += 1;
    }
    stop = start;
    while (stop < len && !(str[stop] == '\n')) {
      stop += 1;
    }
    if (stop > start) {
      lines.push_back(std::string(str, start, stop - start));
    }
    start = stop;
  }
  return lines;
}

inline bool parse_char(char& c, long& cur, const std::string& data)
{
  if ((long)data.size() <= cur) {
    return false;
  } else {
    c = data[cur];
    cur += 1;
    return true;
  }
}

inline bool parse_string(std::string& str, long& cur, const std::string& data)
{
  char c;
  if (!parse_char(c, cur, data) or c != '"') {
    return false;
  } else {
    const long start = cur;
    char c;
    while (parse_char(c, cur, data) and c != '"') {
    }
    str = std::string(data, start, cur - start - 1);
    return data[cur - 1] == '"' && cur > start;
  }
}

inline bool parse_long(long& num, long& cur, const std::string& data)
{
  const long start = cur;
  char c;
  while (parse_char(c, cur, data)) {
    if (not(('0' <= c and c <= '9') or (c == '-') or (c == '+'))) {
      cur -= 1;
      break;
    }
  }
  if (cur <= start) {
    return false;
  } else {
    const std::string str = std::string(data, start, cur - start);
    num = read_long(str);
    return true;
  }
}

inline bool parse_double(double& num, long& cur, const std::string& data)
{
  const long start = cur;
  char c;
  while (parse_char(c, cur, data)) {
    if (not(('0' <= c and c <= '9') or (c == '-') or (c == '+') or (c == '.') or
            (c == 'e') or (c == 'E'))) {
      cur -= 1;
      break;
    }
  }
  if (cur <= start) {
    return false;
  } else {
    const std::string str = std::string(data, start, cur - start);
    num = read_double(str);
    return true;
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

inline std::vector<double> read_doubles(const std::string& str)
{
  const std::vector<std::string> strs = split_line_with_spaces(str);
  std::vector<double> ret(strs.size());
  for (size_t i = 0; i < strs.size(); ++i) {
    ret[i] = read_double(strs[i]);
  }
  return ret;
}

inline std::vector<long> read_longs(const std::string& str)
{
  const std::vector<std::string> strs = split_line_with_spaces(str);
  std::vector<long> ret(strs.size());
  for (size_t i = 0; i < strs.size(); ++i) {
    ret[i] = read_long(strs[i]);
  }
  return ret;
}

inline std::string info_get_prop(const std::vector<std::string>& lines,
                                 const std::string& prop)
{
  for (size_t i = 0; i < lines.size(); ++i) {
    if (lines[i].compare(0, prop.size(), prop) == 0) {
      return std::string(lines[i], prop.size());
    }
  }
  return std::string("");
}

inline std::string info_get_prop(const std::vector<std::string>& lines,
                                 const std::string& prop,
                                 const std::string& prop1)
{
  const std::string ret = info_get_prop(lines, prop);
  if (ret != std::string("")) {
    return ret;
  } else {
    return info_get_prop(lines, prop1);
  }
}

API inline FILE*& get_output_file()
{
  static FILE* out = stdout;
  return out;
}

API inline FILE*& get_monitor_file()
{
  static FILE* out = NULL;
  return out;
}

API inline uint64_t& get_output_level()
{
  static uint64_t x = UINT64_MAX;
  return x;
}

inline void display(const std::string& str, FILE* fp = NULL)
{
  if (NULL == fp) {
    if (get_output_level() == 0) {
      return;
    }
    fp = get_monitor_file();
    if (NULL != fp) {
      fprintf(fp, "%s", str.c_str());
    }
    fp = get_output_file();
  }
  if (NULL != fp) {
    fprintf(fp, "%s", str.c_str());
  }
}

inline void displayln(const std::string& str, FILE* fp = NULL)
{
  if (NULL == fp) {
    if (get_output_level() == 0) {
      return;
    }
    fp = get_monitor_file();
    if (NULL != fp) {
      fprintf(fp, "%s\n", str.c_str());
    }
    fp = get_output_file();
  }
  if (NULL != fp) {
    fprintf(fp, "%s\n", str.c_str());
  }
}

inline const char* get_c_str(const std::string& str) { return str.c_str(); }

//////////////////////////////////////////////////////////////////

inline void fdisplay(FILE* fp, const std::string& str)
{
  fprintf(fp, "%s", str.c_str());
}

inline void fdisplayln(FILE* fp, const std::string& str)
{
  fprintf(fp, "%s\n", str.c_str());
}

}  // namespace qlat
