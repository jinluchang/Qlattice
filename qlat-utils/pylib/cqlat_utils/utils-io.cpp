#include "lib.h"

EXPORT(flush, {
  using namespace qlat;
  fflush(get_output_file());
  Py_RETURN_NONE;
});
