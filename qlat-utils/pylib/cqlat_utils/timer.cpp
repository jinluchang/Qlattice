#include "lib.h"

#include <qlat-utils/timer.h>

EXPORT(get_id_node, {
  using namespace qlat;
  return py_convert(get_id_node());
})

EXPORT(get_num_node, {
  using namespace qlat;
  return py_convert(get_num_node());
})

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
})

EXPORT(free_timer, {
  using namespace qlat;
  return free_obj<Timer>(args);
})

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
})

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
})

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
})

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
})

EXPORT(timer_autodisplay, {
  using namespace qlat;
  Timer::autodisplay();
  Py_RETURN_NONE;
})

EXPORT(timer_reset, {
  using namespace qlat;
  long max_call_times_for_always_show_info = -1;
  if (!PyArg_ParseTuple(args, "|l", &max_call_times_for_always_show_info)) {
    return NULL;
  }
  Timer::reset(max_call_times_for_always_show_info);
  Py_RETURN_NONE;
})

EXPORT(timer_fork, {
  using namespace qlat;
  long max_call_times_for_always_show_info = -1;
  if (!PyArg_ParseTuple(args, "|l", &max_call_times_for_always_show_info)) {
    return NULL;
  }
  Timer::fork(max_call_times_for_always_show_info);
  Py_RETURN_NONE;
})

EXPORT(timer_merge, {
  using namespace qlat;
  Timer::merge();
  Py_RETURN_NONE;
})

EXPORT(timer_display_stack_always, {
  using namespace qlat;
  Timer::display_stack_always();
  Py_RETURN_NONE;
})

EXPORT(timer_display_stack, {
  using namespace qlat;
  Timer::display_stack();
  Py_RETURN_NONE;
})

EXPORT(verbose_level, {
  using namespace qlat;
  PyObject* p_level = NULL;
  if (!PyArg_ParseTuple(args, "|O", &p_level)) {
    return NULL;
  }
  if (NULL != p_level) {
    const long level = py_convert_data<long>(p_level);
    verbose_level() = level;
  }
  return py_convert(verbose_level());
})

EXPORT(get_actual_start_time, {
  using namespace qlat;
  PyObject* p_time = NULL;
  if (!PyArg_ParseTuple(args, "|O", &p_time)) {
    return NULL;
  }
  if (NULL != p_time) {
    const double time = py_convert_data<double>(p_time);
    get_actual_start_time() = time;
  }
  return py_convert(get_actual_start_time());
})

EXPORT(get_start_time, {
  using namespace qlat;
  PyObject* p_time = NULL;
  if (!PyArg_ParseTuple(args, "|O", &p_time)) {
    return NULL;
  }
  if (NULL != p_time) {
    const double time = py_convert_data<double>(p_time);
    get_start_time() = time;
  }
  return py_convert(get_start_time());
})

EXPORT(get_time, {
  using namespace qlat;
  return py_convert(get_time());
})
