#include "convert.h"
#include "dispatch.h"
#include "exceptions.h"

EXPORT(mk_gauge_action, {
  using namespace qlat;
  double beta = 0.0;
  double c1 = 0.0;
  if (!PyArg_ParseTuple(args, "d|d", &beta, &c1)) {
    return NULL;
  }
  GaugeAction* pga = new GaugeAction(beta, c1);
  return py_convert((void*)pga);
});

EXPORT(free_gauge_action, {
  using namespace qlat;
  GaugeAction* pga = NULL;
  if (!PyArg_ParseTuple(args, "l", &pga)) {
    return NULL;
  }
  pqassert(pga != NULL);
  delete pga;
  Py_RETURN_NONE;
});

EXPORT(get_beta_gauge_action, {
  using namespace qlat;
  GaugeAction* pga = NULL;
  if (!PyArg_ParseTuple(args, "l", &pga)) {
    return NULL;
  }
  pqassert(pga != NULL);
  const GaugeAction& ga = *pga;
  return py_convert(ga.beta);
});

EXPORT(get_c1_gauge_action, {
  using namespace qlat;
  GaugeAction* pga = NULL;
  if (!PyArg_ParseTuple(args, "l", &pga)) {
    return NULL;
  }
  pqassert(pga != NULL);
  const GaugeAction& ga = *pga;
  return py_convert(ga.c1);
});
