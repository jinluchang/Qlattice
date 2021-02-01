#include "convert.h"
#include "dispatch.h"
#include "exceptions.h"

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

