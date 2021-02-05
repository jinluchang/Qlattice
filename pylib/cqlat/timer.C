#include "lib.h"

EXPORT(mk_timer, {
  using namespace qlat;
  PyObject* p_fname = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_fname)) {
    return NULL;
  }
  std::string fname;
  py_convert(fname, p_fname);
  Timer* ptimer = new Timer(fname, false);
  return py_convert((void*)ptimer);
});

EXPORT(free_timer, {
  using namespace qlat;
  PyObject* p_timer = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_timer)) {
    return NULL;
  }
  Timer& timer = py_convert_type<Timer>(p_timer);
  delete &timer;
  Py_RETURN_NONE;
});

EXPORT(start_timer, {
  using namespace qlat;
  PyObject* p_timer = NULL;
  bool is_verbose = false;
  if (!PyArg_ParseTuple(args, "O|b", &p_timer, &is_verbose)) {
    return NULL;
  }
  Timer& timer = py_convert_type<Timer>(p_timer);
  timer.start(is_verbose);
  Py_RETURN_NONE;
});

EXPORT(stop_timer, {
  using namespace qlat;
  PyObject* p_timer = NULL;
  bool is_verbose = false;
  if (!PyArg_ParseTuple(args, "O|b", &p_timer, &is_verbose)) {
    return NULL;
  }
  Timer& timer = py_convert_type<Timer>(p_timer);
  timer.stop(is_verbose);
  Py_RETURN_NONE;
});

EXPORT(set_flops_timer, {
  using namespace qlat;
  PyObject* p_timer = NULL;
  long flops = false;
  if (!PyArg_ParseTuple(args, "Ol", &p_timer, &flops)) {
    return NULL;
  }
  Timer& timer = py_convert_type<Timer>(p_timer);
  timer.flops = flops;
  Py_RETURN_NONE;
});

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

