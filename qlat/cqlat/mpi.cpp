#include "lib.h"

EXPORT(get_coor_node, {
  using namespace qlat;
  return py_convert(get_coor_node());
})

EXPORT(get_size_node, {
  using namespace qlat;
  return py_convert(get_size_node());
})

EXPORT(sync_node, {
  using namespace qlat;
  sync_node();
  Py_RETURN_NONE;
})

EXPORT(glb_sum_long, {
  using namespace qlat;
  long x = 0;
  if (!PyArg_ParseTuple(args, "l", &x)) {
    return NULL;
  }
  glb_sum(x);
  return py_convert(x);
})

EXPORT(glb_sum_double, {
  using namespace qlat;
  double x = 0.0;
  if (!PyArg_ParseTuple(args, "d", &x)) {
    return NULL;
  }
  glb_sum(x);
  return py_convert(x);
})

EXPORT(glb_sum_complex, {
  using namespace qlat;
  Complex x = 0.0;
  if (!PyArg_ParseTuple(args, "D", &x)) {
    return NULL;
  }
  glb_sum(x);
  return py_convert(x);
})
