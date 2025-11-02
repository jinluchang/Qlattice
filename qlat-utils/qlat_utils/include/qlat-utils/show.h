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

#undef NDEBUG

#include <cassert>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

#define PSTRX(x) #x

#define PSTR(x) PSTRX(x)

#define API __attribute__((visibility("default")))

namespace qlat
{  //

using Long = int64_t;

using Int = int32_t;

using Char = int8_t;

using RealD = double;

using RealF = float;

using Real = RealD;  // default Real type should not change

inline std::string vssprintf(const char* fmt, va_list args)
{
  char* cstr;
  Int ret = vasprintf(&cstr, fmt, args);
  if (ret < 0) {
    assert(false);
  }
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

inline std::string show(const Int& x) { return ssprintf("%d", x); }

inline std::string show(const uint32_t& x) { return ssprintf("%u", x); }

inline std::string show(const Long& x) { return ssprintf("%ld", x); }

inline std::string show(const uint64_t& x) { return ssprintf("%lu", x); }

inline std::string show(const RealD& x) { return ssprintf("%24.17E", x); }

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
  for (Long i = 0; i < (Long)vec.size(); ++i) {
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

inline Long read_long(const std::string& str)
{
  Long ret = 0;
  reads(ret, str);
  return ret;
}

inline RealD read_double(const std::string& str)
{
  RealD ret = 0.0;
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
  Long cur = s.size() - 1;
  while ((cur >= 0) and (s[cur] == '\n' or s[cur] == '\r')) {
    cur -= 1;
  }
  return std::string(s, 0, cur + 1);
}

inline bool is_space(const char c)
{
  return (c == ' ') or (c == '\n') or (c == '\r') or (c == '\t');
}

inline bool parse_end(Long& cur, const std::string& data)
{
  assert(cur <= (Long)data.size());
  return cur == (Long)data.size();
}

inline bool parse_char(char& c, Long& cur, const std::string& data)
{
  assert(cur <= (Long)data.size());
  if ((Long)data.size() <= cur) {
    c = 0;
    return false;
  }
  c = data[cur];
  cur += 1;
  return true;
}

inline bool parse_char_not(char& c, Long& cur, const std::string& data,
                           const char c_match)
{
  if (not parse_char(c, cur, data)) {
    return false;
  }
  if (c != c_match) {
    return true;
  }
  cur -= 1;
  return false;
}

inline bool parse_char_space(char& c, Long& cur, const std::string& data)
{
  if (not parse_char(c, cur, data)) {
    return false;
  }
  if (is_space(c)) {
    return true;
  }
  cur -= 1;
  return false;
}

inline bool parse_char_not_space(char& c, Long& cur, const std::string& data)
{
  if (not parse_char(c, cur, data)) {
    return false;
  }
  if (not is_space(c)) {
    return true;
  }
  cur -= 1;
  return false;
}

inline bool parse_literal(Long& cur, const std::string& data,
                          const char c_match)
{
  char c;
  if (not parse_char(c, cur, data)) {
    return false;
  }
  if (c == c_match) {
    return true;
  }
  cur -= 1;
  return false;
}

inline bool parse_literal(Long& cur, const std::string& data,
                          const std::string& literal)
{
  assert(cur <= (Long)data.size());
  if (0 != data.compare(cur, literal.size(), literal)) {
    return false;
  }
  cur += literal.size();
  return true;
}

inline bool parse_len(std::string& str, Long& cur, const std::string& data,
                      const Long& len)
// not including the '"' char
{
  assert(cur <= (Long)data.size());
  if (cur + len > (Long)data.size()) {
    return false;
  }
  str = data.substr(cur, len);
  cur += len;
  return true;
}

inline bool parse_string(std::string& str, Long& cur, const std::string& data)
// not including the '"' char
{
  assert(cur <= (Long)data.size());
  const Long initial = cur;
  char c;
  if (not parse_literal(cur, data, '"')) {
    str = "";
    assert(cur == initial);
    return false;
  }
  const Long start = cur;
  while (parse_char_not(c, cur, data, '"')) {
  }
  if (not parse_literal(cur, data, '"')) {
    str = "";
    cur = initial;
    return false;
  }
  str = std::string(data, start, cur - start - 1);
  return true;
}

inline bool parse_line(std::string& str, Long& cur, const std::string& data)
// include ending '\n' (if data has it)
{
  assert(cur <= (Long)data.size());
  const Long initial = cur;
  if ((Long)data.size() <= cur) {
    str = "";
    assert(cur == initial);
    return false;
  }
  const Long start = cur;
  char c;
  while (parse_char_not(c, cur, data, '\n')) {
  }
  if (not parse_literal(cur, data, '\n')) {
    str = "";
    cur = initial;
    return false;
  }
  str = std::string(data, start, cur - start);
  return true;
}

inline bool parse_word(std::string& str, Long& cur, const std::string& data)
// not include initial spaces and do not parse space after it
{
  assert(cur <= (Long)data.size());
  const Long initial = cur;
  if ((Long)data.size() <= cur) {
    str = "";
    assert(cur == initial);
    return false;
  }
  char c;
  while (parse_char_space(c, cur, data)) {
  }
  const Long start = cur;
  while (parse_char_not_space(c, cur, data)) {
  }
  if (cur <= start) {
    str = "";
    cur = initial;
    return false;
  }
  str = std::string(data, start, cur - start);
  return true;
}

inline bool parse_long(Long& num, Long& cur, const std::string& data)
{
  assert(cur <= (Long)data.size());
  const Long initial = cur;
  const Long start = cur;
  char c;
  while (parse_char(c, cur, data)) {
    if (not((('0' <= c) and (c <= '9')) or (c == '-') or (c == '+'))) {
      cur -= 1;
      break;
    }
  }
  if (cur <= start) {
    num = 0;
    cur = initial;
    return false;
  }
  const std::string str = std::string(data, start, cur - start);
  num = read_long(str);
  return true;
}

inline bool parse_double(RealD& num, Long& cur, const std::string& data)
{
  assert(cur <= (Long)data.size());
  const Long initial = cur;
  const Long start = cur;
  char c;
  while (parse_char(c, cur, data)) {
    if (not((('0' <= c) and (c <= '9')) or (c == '-') or (c == '+') or
            (c == '.') or (c == 'e') or (c == 'E'))) {
      cur -= 1;
      break;
    }
  }
  if (cur <= start) {
    num = 0;
    cur = initial;
    return false;
  }
  const std::string str = std::string(data, start, cur - start);
  num = read_double(str);
  return true;
}

inline std::vector<std::string> split_into_lines(const std::string& str)
// results do have newline at end (if the original str has it)
{
  const size_t len = str.length();
  std::vector<std::string> lines;
  size_t start = 0;
  size_t stop = 0;
  while (stop < len) {
    while (stop < len && !(str[stop] == '\n')) {
      stop += 1;
    }
    if (stop < len && str[stop] == '\n') {
      stop += 1;
    }
    assert(stop > start);
    lines.push_back(std::string(str, start, stop - start));
    start = stop;
  }
  return lines;
}

inline std::string merge_lines(const std::vector<std::string>& lines)
{
  std::string ret;
  for (Long i = 0; i < (Long)lines.size(); ++i) {
    ret += lines[i];
  }
  return ret;
}

inline std::vector<std::string> split_line_with_spaces(const std::string& str)
// results do not have spaces
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

inline std::vector<RealD> read_doubles(const std::string& str)
{
  const std::vector<std::string> strs = split_line_with_spaces(str);
  std::vector<RealD> ret(strs.size());
  for (size_t i = 0; i < strs.size(); ++i) {
    ret[i] = read_double(strs[i]);
  }
  return ret;
}

inline std::vector<Long> read_longs(const std::string& str)
{
  const std::vector<std::string> strs = split_line_with_spaces(str);
  std::vector<Long> ret(strs.size());
  for (size_t i = 0; i < strs.size(); ++i) {
    ret[i] = read_long(strs[i]);
  }
  return ret;
}

inline std::string info_get_prop(const std::vector<std::string>& lines,
                                 const std::string& prop)
{
  for (size_t i = 0; i < lines.size(); ++i) {
    if (0 == lines[i].compare(0, prop.size(), prop)) {
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

inline bool starts_with(const std::string& str, const std::string& prefix) {
  return str.size() >= prefix.size() &&
    str.compare(0, prefix.size(), prefix) == 0;
}

inline bool ends_with(const std::string& str, const std::string& suffix) {
  if (str.size() < suffix.size()) {
    return false;
  }
  return (str.substr(str.size() - suffix.size()) == suffix);
}

inline std::string remove_prefix(const std::string& str, const std::string& prefix) {
  if (str.size() >= prefix.size() &&
      str.compare(0, prefix.size(), prefix) == 0) {
    return str.substr(prefix.size());
  }
  return str;
}

inline std::string remove_suffix(const std::string& str, const std::string& suffix) {
  if (str.size() >= suffix.size() &&
      str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0) {
    return str.substr(0, str.size() - suffix.size());
  }
  return str;
}

inline void display_c_stdout(const std::string& str)
{
  fwrite(str.c_str(), 1, str.size(), stdout);
  fflush(stdout);
}

inline void displayln_c_stdout(const std::string& str)
{
  display_c_stdout(str + "\n");
}

using DisplayPtr = void (*)(const std::string& str);

API inline DisplayPtr& get_display_ptr()
{
  static DisplayPtr ptr = display_c_stdout;
  return ptr;
}

inline void set_display_ptr() { get_display_ptr() = display_c_stdout; }

inline void set_display_ptr(DisplayPtr f) { get_display_ptr() = f; }

inline const char* get_c_str(const std::string& str) { return str.c_str(); }

//////////////////////////////////////////////////////////////////

inline void fdisplay(FILE* fp, const std::string& str)
{
  fwrite(str.c_str(), 1, str.size(), fp);
}

inline void fdisplayln(FILE* fp, const std::string& str)
{
  fwrite(str.c_str(), 1, str.size(), fp);
  fwrite("\n", 1, 1, fp);
}

}  // namespace qlat