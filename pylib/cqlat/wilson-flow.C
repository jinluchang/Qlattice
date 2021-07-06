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
});

EXPORT(gf_energy_density, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_gf)) {
    return NULL;
  }
  const GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  return py_convert(gf_energy_density(gf));
});

EXPORT(gf_ape_smear, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  PyObject* p_gf0 = NULL;
  double alpha = 0.0;
  long steps = 1;
  if (!PyArg_ParseTuple(args, "OOd|l", &p_gf, &p_gf0, &alpha, &steps)) {
    return NULL;
  }
  GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  const GaugeField& gf0 = py_convert_type<GaugeField>(p_gf0);
  gf_ape_smear(gf, gf0, alpha, steps);
  Py_RETURN_NONE;
});

EXPORT(gf_hyp_smear, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  PyObject* p_gf0 = NULL;
  // values in paper is 0.75 0.6 0.3
  // 10.1103/PhysRevD.64.034504 Eq(4)
  double alpha1 = 0.0;
  double alpha2 = 0.0;
  double alpha3 = 0.0;
  if (!PyArg_ParseTuple(args, "OOddd", &p_gf, &p_gf0, &alpha1, &alpha2,
                        &alpha3)) {
    return NULL;
  }
  GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  const GaugeField& gf0 = py_convert_type<GaugeField>(p_gf0);
  gf_hyp_smear(gf, gf0, alpha1, alpha2, alpha3);
  Py_RETURN_NONE;
});
