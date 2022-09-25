#include "lib.h"

EXPORT(mk_gauge_action, {
  using namespace qlat;
  double beta = 0.0;
  double c1 = 0.0;
  if (!PyArg_ParseTuple(args, "d|d", &beta, &c1)) {
    return NULL;
  }
  GaugeAction* pga = new GaugeAction(beta, c1);
  return py_convert((void*)pga);
})

EXPORT(free_gauge_action, {
  using namespace qlat;
  return free_obj<GaugeAction>(args);
})

EXPORT(set_gauge_action, {
  using namespace qlat;
  return set_obj<GaugeAction>(args);
})

EXPORT(get_beta_gauge_action, {
  using namespace qlat;
  PyObject* p_ga = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_ga)) {
    return NULL;
  }
  const GaugeAction& ga = py_convert_type<GaugeAction>(p_ga);
  return py_convert(ga.beta);
})

EXPORT(get_c1_gauge_action, {
  using namespace qlat;
  PyObject* p_ga = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_ga)) {
    return NULL;
  }
  const GaugeAction& ga = py_convert_type<GaugeAction>(p_ga);
  return py_convert(ga.c1);
})
