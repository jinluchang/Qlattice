#include <qlat-ext/core.h>

namespace qlat
{

int hello()
{
  TIMER_VERBOSE("hello");
  displayln_info("hello!");
  return 0;
}

}  // namespace qlat
