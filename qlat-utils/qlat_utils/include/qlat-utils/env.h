#pragma once

#include <dirent.h>
#include <qlat-utils/show.h>

namespace qlat
{  //

std::string get_env(const std::string& var_name);

std::string get_env_default(const std::string& var_name, const std::string& x0);

double get_env_double_default(const std::string& var_name, const double x0);

long get_env_long_default(const std::string& var_name, const long x0);

long get_verbose_level_default();

double get_time_limit_default();

double get_time_budget_default();

API inline long& get_verbose_level()
// qlat parameter
{
  static long level = get_verbose_level_default();
  return level;
}

API inline double& get_time_limit()
// qlat parameter
{
  static double limit = get_time_limit_default();
  return limit;
}

API inline double& get_time_budget()
// qlat parameter
{
  static double budget = get_time_budget_default();
  return budget;
}

API inline long& get_qar_multi_vol_max_size()
// qlat parameter
// size in bytes
{
  static long size = get_env_long_default("q_qar_multi_vol_max_size",
                                          500L * 1000L * 1000L * 1000L);
  return size;
}

API inline mode_t& default_dir_mode()
// qlat parameter
{
  static mode_t mode = 0775;
  return mode;
}

API inline long& write_from_qfile_chunk_size()
// qlat parameter
// size in bytes
{
  static long size = 512 * 1024;
  return size;
}

}  // namespace qlat
