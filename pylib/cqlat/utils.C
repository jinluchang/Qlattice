#include "lib.h"

EXPORT(timer_display, {
  using namespace qlat;
  PyObject* p_str = NULL;
  if (!PyArg_ParseTuple(args, "|O", &p_str)) {
    return NULL;
  }
  std::string str;
  if (p_str != NULL) {
    py_convert(str, p_str);
  }
  Timer::display(str);
  Py_RETURN_NONE;
});

EXPORT(timer_display_stack_always, {
  using namespace qlat;
  Timer::display_stack_always();
  Py_RETURN_NONE;
});

EXPORT(timer_display_stack, {
  using namespace qlat;
  Timer::display_stack();
  Py_RETURN_NONE;
});

EXPORT(glb_sum_long, {
  using namespace qlat;
  long x = 0;
  if (!PyArg_ParseTuple(args, "l", &x)) {
    return NULL;
  }
  glb_sum(x);
  return py_convert(x);
});

EXPORT(glb_sum_double, {
  using namespace qlat;
  double x = 0.0;
  if (!PyArg_ParseTuple(args, "d", &x)) {
    return NULL;
  }
  glb_sum(x);
  return py_convert(x);
});

EXPORT(glb_sum_complex, {
  using namespace qlat;
  Complex x = 0.0;
  if (!PyArg_ParseTuple(args, "D", &x)) {
    return NULL;
  }
  glb_sum(x);
  return py_convert(x);
});
