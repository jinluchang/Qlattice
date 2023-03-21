#include <qlat/env.h>
#include <qlat-utils/timer.h>

namespace qlat
{  //

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

int get_field_init_from_env()
{
  std::string tag = get_env_default("q_field_init", "fast");
  if (tag == "fast") {
    displayln_info("set q_field_init=fast.");
    return 0;  // do not do anything
  } else if (tag == "zero") {
    displayln_info("set q_field_init=zero.");
    return 1;  // set_zero
  } else if (tag == "random") {
    displayln_info("set q_field_init=random.");
    return 2;  // set rand
  } else {
    qassert(false);
    return -1;
  }
}

}  // namespace qlat
