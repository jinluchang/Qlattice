#pragma once

#include <sstream>
#include <string>
#include <cstdarg>
#include <cstring>
#include <cstdlib>
#include <cstdio>

inline std::string vssprintf(const char* fmt, va_list args)
{
  std::string str;
  char* cstr;
  vasprintf(&cstr, fmt, args);
  str += std::string(cstr);
  std::free(cstr);
  return str;
}

inline std::string ssprintf(const char* fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  return vssprintf(fmt, args);
}

inline std::string show()
{
  return "";
}

inline std::string show(const int& x)
{
  return ssprintf("%d", x);
}

inline std::string show(const unsigned int& x)
{
  return ssprintf("%u", x);
}

inline std::string show(const long& x)
{
  return ssprintf("%ld", x);
}

inline std::string show(const unsigned long& x)
{
  return ssprintf("%lu", x);
}

inline std::string show(const double& x)
{
  return ssprintf("%23.16E", x);
}

inline std::string show(const bool& x)
{
  return x ? "true" : "false";
}

inline std::string show(const std::string& x)
{
  std::ostringstream out;
  out << x;
  return out.str();
}
