#include <qlat/flowed-hmc.h>
#include "lib.h"

EXPORT(mk_flow_info, {
  using namespace qlat;
  FlowInfo* pfi = new FlowInfo();
  return py_convert((void*)pfi);
});

EXPORT(free_flow_info, {
  using namespace qlat;
  PyObject* p_fi = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_fi)) {
    return NULL;
  }
  FlowInfo& fi = py_convert_type<FlowInfo>(p_fi);
  delete &fi;
  Py_RETURN_NONE;
});

EXPORT(show_flow_info, {
  using namespace qlat;
  PyObject* p_fi = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_fi)) {
    return NULL;
  }
  const FlowInfo& fi = py_convert_type<FlowInfo>(p_fi);
  return py_convert(show(fi));
});

EXPORT(add_flow_flow_info, {
  using namespace qlat;
  PyObject* p_fi = NULL;
  int eo = 0;
  int mu = 0;
  double epsilon = 0.0;
  int flow_size = 1;
  if (!PyArg_ParseTuple(args, "Oiid|i", &p_fi, &eo, &mu, &epsilon,
                        &flow_size)) {
    return NULL;
  }
  FlowInfo& fi = py_convert_type<FlowInfo>(p_fi);
  fi.v.push_back(FlowStepInfo(eo, mu, epsilon, flow_size));
  Py_RETURN_NONE;
});

EXPORT(add_rand_order_flow2_flow_info, {
  using namespace qlat;
  PyObject* p_fi = NULL;
  PyObject* p_rs = NULL;
  double epsilon = 0.0;
  double epsilon2 = 0.0;
  if (!PyArg_ParseTuple(args, "OOdd", &p_fi, &p_rs, &epsilon, &epsilon2)) {
    return NULL;
  }
  FlowInfo& fi = py_convert_type<FlowInfo>(p_fi);
  const RngState& rs = py_convert_type<RngState>(p_rs);
  vector_append(fi.v, mk_flow_info_step(rs, epsilon, epsilon2).v);
  Py_RETURN_NONE;
});

EXPORT(add_rand_order_flow_flow_info, {
  using namespace qlat;
  PyObject* p_fi = NULL;
  PyObject* p_rs = NULL;
  double epsilon = 0.0;
  if (!PyArg_ParseTuple(args, "OOd", &p_fi, &p_rs, &epsilon)) {
    return NULL;
  }
  FlowInfo& fi = py_convert_type<FlowInfo>(p_fi);
  const RngState& rs = py_convert_type<RngState>(p_rs);
  vector_append(fi.v, mk_flow_info_step(rs, epsilon).v);
  Py_RETURN_NONE;
});

EXPORT(gf_flow, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  PyObject* p_gf0 = NULL;
  PyObject* p_fi = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_gf, &p_gf0, &p_fi)) {
    return NULL;
  }
  GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  const GaugeField& gf0 = py_convert_type<GaugeField>(p_gf0);
  const FlowInfo& fi = py_convert_type<FlowInfo>(p_fi);
  gf_flow(gf, gf0, fi);
  Py_RETURN_NONE;
});

EXPORT(gf_flow_inv, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  PyObject* p_gf1 = NULL;
  PyObject* p_fi = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_gf, &p_gf1, &p_fi)) {
    return NULL;
  }
  GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  const GaugeField& gf1 = py_convert_type<GaugeField>(p_gf1);
  const FlowInfo& fi = py_convert_type<FlowInfo>(p_fi);
  gf_flow_inv(gf, gf1, fi);
  Py_RETURN_NONE;
});

EXPORT(gf_hamilton_flowed_node, {
  using namespace qlat;
  PyObject* p_gf0 = NULL;
  PyObject* p_ga = NULL;
  PyObject* p_fi = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_gf0, &p_ga, &p_fi)) {
    return NULL;
  }
  const GaugeField& gf0 = py_convert_type<GaugeField>(p_gf0);
  const GaugeAction& ga = py_convert_type<GaugeAction>(p_ga);
  const FlowInfo& fi = py_convert_type<FlowInfo>(p_fi);
  const double ret = gf_hamilton_flowed_node(gf0, ga, fi);
  return py_convert(ret);
});

EXPORT(set_gm_force_flowed, {
  using namespace qlat;
  PyObject* p_gm_force = NULL;
  PyObject* p_gf0 = NULL;
  PyObject* p_ga = NULL;
  PyObject* p_fi = NULL;
  if (!PyArg_ParseTuple(args, "OOOO", &p_gm_force, &p_gf0, &p_ga, &p_fi)) {
    return NULL;
  }
  GaugeMomentum& gm_force = py_convert_type<GaugeMomentum>(p_gm_force);
  const GaugeField& gf0 = py_convert_type<GaugeField>(p_gf0);
  const GaugeAction& ga = py_convert_type<GaugeAction>(p_ga);
  const FlowInfo& fi = py_convert_type<FlowInfo>(p_fi);
  set_gm_force_flowed(gm_force, gf0, ga, fi);
  Py_RETURN_NONE;
});

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
});
