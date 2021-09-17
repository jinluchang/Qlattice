#include "lib.h"

EXPORT(hello_world, {
  using namespace qlat;
  displayln_info("Hello world!");
  Py_RETURN_NONE;
});
