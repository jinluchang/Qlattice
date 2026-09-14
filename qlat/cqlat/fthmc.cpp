#include <qlat/flowed-hmc.h>

#include "lib.h"

EXPORT(free_flow_info, {
  using namespace qlat;
  PyObject* p_fi = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_fi)) {
    return NULL;
  }
  FlowInfo& fi = py_convert_type<FlowInfo>(p_fi);
  delete &fi;
  Py_RETURN_NONE;
})

EXPORT(add_flow_flow_info, {
  using namespace qlat;
  PyObject* p_fi = NULL;
  Int eo = 0;
  Int mu = 0;
  double epsilon = 0.0;
  Int flow_size = 1;
  if (!PyArg_ParseTuple(args, "Oiid|i", &p_fi, &eo, &mu, &epsilon,
                        &flow_size)) {
    return NULL;
  }
  FlowInfo& fi = py_convert_type<FlowInfo>(p_fi);
  fi.v.push_back(FlowStepInfo(eo, mu, epsilon, flow_size));
  Py_RETURN_NONE;
})

EXPORT(set_gm_force_flowed_no_det, {
  using namespace qlat;
  PyObject* p_gm_force = NULL;
  PyObject* p_gm_force_pre = NULL;
  PyObject* p_gf0 = NULL;
  PyObject* p_fi = NULL;
  if (!PyArg_ParseTuple(args, "OOOO", &p_gm_force, &p_gm_force_pre, &p_gf0,
                        &p_fi)) {
    return NULL;
  }
  GaugeMomentum& gm_force = py_convert_type<GaugeMomentum>(p_gm_force);
  const GaugeMomentum& gm_force_pre =
      py_convert_type<GaugeMomentum>(p_gm_force_pre);
  const GaugeField& gf0 = py_convert_type<GaugeField>(p_gf0);
  const FlowInfo& fi = py_convert_type<FlowInfo>(p_fi);
  set_gm_force_flowed_no_det(gm_force, gm_force_pre, gf0, fi);
  Py_RETURN_NONE;
})
