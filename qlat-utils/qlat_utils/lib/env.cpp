#include <qlat-utils/env.h>
#include <qlat-utils/timer.h>

namespace qlat
{  //

std::string get_env(const std::string& var_name)
{
  const char* value = getenv(var_name.c_str());
  if (value == NULL) {
    return std::string();
  } else {
    return std::string(value);
  }
}

std::string get_env_default(const std::string& var_name, const std::string& x0)
{
  const std::string val = get_env(var_name);
  if (val == "") {
    displayln_info(0,
                   ssprintf("%s='%s' (default)", var_name.c_str(), x0.c_str()));
    return x0;
  } else {
    displayln_info(0, ssprintf("%s='%s'", var_name.c_str(), val.c_str()));
    return val;
  }
}

double get_env_double_default(const std::string& var_name, const double x0)
{
  const std::string val = get_env(var_name);
  double x;
  if (val == "") {
    x = x0;
    displayln_info(0, ssprintf("%s=%lG (default)", var_name.c_str(), x));
  } else {
    x = read_double(val);
    displayln_info(0, ssprintf("%s=%lG", var_name.c_str(), x));
  }
  return x;
}

long get_env_long_default(const std::string& var_name, const long x0)
{
  const std::string val = get_env(var_name);
  long x;
  if (val == "") {
    x = x0;
    displayln_info(0, ssprintf("%s=%ld (default)", var_name.c_str(), x));
  } else {
    x = read_long(val);
    displayln_info(0, ssprintf("%s=%ld", var_name.c_str(), x));
  }
  return x;
}

long get_verbose_level_default()
{
  const long x0 = -1;  // default verbose_level
  const std::string var_name = "q_verbose";
  const std::string val = get_env(var_name);
  long x;
  if (val == "") {
    x = x0;
  } else {
    x = read_long(val);
  }
  return x;
}

double get_time_limit_default()
{
  const double time_limit_default = 12.0 * 3600.0;
  if (get_env("q_end_time") == "") {
    return get_env_double_default("q_time_limit", time_limit_default);
  } else {
    return get_env_double_default(
               "q_end_time", get_actual_start_time() + time_limit_default) -
           get_actual_start_time();
  }
}

double get_time_budget_default()
{
  double budget = get_env_double_default("q_time_budget", 15.0 * 60.0);
  return budget;
}

}  // namespace qlat
