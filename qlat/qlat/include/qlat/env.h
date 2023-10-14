#pragma once

#include <qlat-utils/timer.h>

namespace qlat
{  //

double get_time_limit_default();

API inline double& get_time_limit()
// qlat parameter
{
  static double limit = get_time_limit_default();
  return limit;
}

API inline double& get_default_budget()
// qlat parameter
{
  static double budget = get_env_double_default("q_budget", 15.0 * 60.0);
  return budget;
}

int get_field_init_from_env();

API inline int& get_field_init()
// qlat parameter
{
  static int t = get_field_init_from_env();
  return t;
}

API inline long& get_max_field_shift_direct_msg_size()
{
  static long size = 1024L * 1024L * 1024L;
  // static long size = 16;
  return size;
}

API inline std::string& get_lock_location()
{
  static std::string path;
  return path;
}

}  // namespace qlat
