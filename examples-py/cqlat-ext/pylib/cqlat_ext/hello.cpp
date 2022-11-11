#include "lib.h"

EXPORT(hello_world, {
  using namespace qlat;
  {
    TIMER_VERBOSE("hello_world")
    displayln_info("Hello world!");
  }
  Py_RETURN_NONE;
})
