#include "lib.h"

EXPORT(gf_wilson_flow_step, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  double epsilon = 0.0;
  double c1 = 0.0;
  if (!PyArg_ParseTuple(args, "Od|d", &p_gf, &epsilon, &c1)) {
    return NULL;
  }
  GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  gf_wilson_flow_step(gf, epsilon, c1);
  Py_RETURN_NONE;
})

EXPORT(gf_energy_density, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_gf)) {
    return NULL;
  }
  const GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  return py_convert(gf_energy_density(gf));
})
