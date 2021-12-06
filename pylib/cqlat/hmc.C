#include "lib.h"

EXPORT(set_rand_gauge_momentum, {
  using namespace qlat;
  PyObject* p_gm = NULL;
  double sigma = 1.0;
  PyObject* p_rng = NULL;
  if (!PyArg_ParseTuple(args, "OdO", &p_gm, &sigma, &p_rng)) {
    return NULL;
  }
  GaugeMomentum& gm = py_convert_type<GaugeMomentum>(p_gm);
  const RngState& rs = py_convert_type<RngState>(p_rng);
  set_rand_gauge_momentum(gm, sigma, rs);
  Py_RETURN_NONE;
});

EXPORT(gm_hamilton_node, {
  using namespace qlat;
  PyObject* p_gm = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_gm)) {
    return NULL;
  }
  const GaugeMomentum& gm = py_convert_type<GaugeMomentum>(p_gm);
  const double ret = gm_hamilton_node(gm);
  return py_convert(ret);
});

EXPORT(gf_hamilton_node, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  PyObject* p_ga = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_gf, &p_ga)) {
    return NULL;
  }
  const GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  const GaugeAction& ga = py_convert_type<GaugeAction>(p_ga);
  const double ret = gf_hamilton_node(gf, ga);
  return py_convert(ret);
});

EXPORT(gf_evolve, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  PyObject* p_gm = NULL;
  double step_size = 0.0;
  if (!PyArg_ParseTuple(args, "OOd", &p_gf, &p_gm, &step_size)) {
    return NULL;
  }
  GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  const GaugeMomentum& gm = py_convert_type<GaugeMomentum>(p_gm);
  gf_evolve(gf, gm, step_size);
  Py_RETURN_NONE;
});

EXPORT(set_gm_force, {
  using namespace qlat;
  PyObject* p_gm_force = NULL;
  PyObject* p_gf = NULL;
  PyObject* p_ga = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_gm_force, &p_gf, &p_ga)) {
    return NULL;
  }
  GaugeMomentum& gm_force = py_convert_type<GaugeMomentum>(p_gm_force);
  const GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  const GaugeAction& ga = py_convert_type<GaugeAction>(p_ga);
  set_gm_force(gm_force, gf, ga);
  Py_RETURN_NONE;
});

EXPORT(sf_hamilton_node, {
  using namespace qlat;
  PyObject* p_sf = NULL;
  PyObject* p_sa = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_sf, &p_sa)) {
    return NULL;
  }
  const Field<double>& sf = py_convert_type<Field<double>>(p_sf);
  const ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  const double ret = sf_hamilton_node(sf, sa);
  return py_convert(ret);
});

EXPORT(sm_hamilton_node, {
  using namespace qlat;
  PyObject* p_sm = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_sm)) {
    return NULL;
  }
  const Field<double>& sm = py_convert_type<Field<double>>(p_sm);
  const double ret = sm_hamilton_node(sm);
  return py_convert(ret);
});

EXPORT(sf_evolve, {
  using namespace qlat;
  PyObject* p_sf = NULL;
  PyObject* p_sm = NULL;
  double step_size = 0.0;
  if (!PyArg_ParseTuple(args, "OOd", &p_sf, &p_sm, &step_size)) {
    return NULL;
  }
  Field<double>& sf = py_convert_type<Field<double>>(p_sf);
  const Field<double>& sm = py_convert_type<Field<double>>(p_sm);
  sf_evolve(sf, sm, step_size);
  Py_RETURN_NONE;
});

EXPORT(set_sm_force, {
  using namespace qlat;
  PyObject* p_sm_force = NULL;
  PyObject* p_sf = NULL;
  PyObject* p_sa = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_sm_force, &p_sf, &p_sa)) {
    return NULL;
  }
  Field<double>& sm_force = py_convert_type<Field<double>>(p_sm_force);
  const Field<double>& sf = py_convert_type<Field<double>>(p_sf);
  const ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  set_sm_force(sm_force, sf, sa);
  Py_RETURN_NONE;
});
