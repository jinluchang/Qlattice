#include <qlat/env.h>
#include <qlat-utils/timer.h>

namespace qlat
{  //

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
