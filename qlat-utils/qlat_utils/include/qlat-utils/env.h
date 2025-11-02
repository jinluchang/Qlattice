#pragma once

#include <dirent.h>
#include <qlat-utils/show.h>

namespace qlat
{  //

std::string get_env(const std::string& var_name);

std::string get_env_default(const std::string& var_name, const std::string& x0);

RealD get_env_double_default(const std::string& var_name, const RealD x0);

Long get_env_long_default(const std::string& var_name, const Long x0);

Long get_verbose_level_default();

RealD get_time_limit_default();

RealD get_time_budget_default();

Long get_qar_multi_vol_max_size_default();

API inline Long& get_verbose_level()
// qlat parameter
{
  static Long level = get_verbose_level_default();
  return level;
}

API inline RealD& get_time_limit()
// qlat parameter
{
  static RealD limit = get_time_limit_default();
  return limit;
}

API inline RealD& get_time_budget()
// qlat parameter
{
  static RealD budget = get_time_budget_default();
  return budget;
}

API inline Long& get_qar_multi_vol_max_size()
// qlat parameter
// size in bytes
{
  static Long size = get_qar_multi_vol_max_size_default();
  return size;
}

API inline mode_t& default_dir_mode()
// qlat parameter
{
  static mode_t mode = 0775;
  return mode;
}

API inline Long& write_from_qfile_chunk_size()
// qlat parameter
// size in bytes
{
  static Long size =
      get_env_long_default("q_write_from_qfile_chunk_size", 4L * 1024L * 1024L);
  return size;
}

API inline std::string& get_fftw_plan_flag()
// qlat parameter
// "estimate" (default), "measure"
{
  static std::string flag = get_env_default("q_fftw_plan_flag", "estimate");
  return flag;
}

}  // namespace qlat
