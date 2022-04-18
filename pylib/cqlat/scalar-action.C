#include "lib.h"

EXPORT(mk_scalar_action, {
  using namespace qlat;
  double m_sq = 1.0;
  double lmbd = 1.0;
  double alpha = 0.0;
  if (!PyArg_ParseTuple(args, "d|d|d", &m_sq, &lmbd, &alpha)) {
    return NULL;
  }
  ScalarAction* psa = new ScalarAction(m_sq, lmbd, alpha);
  return py_convert((void*)psa);
});

EXPORT(free_scalar_action, {
  using namespace qlat;
  return free_obj<ScalarAction>(args);
});

EXPORT(set_scalar_action, {
  using namespace qlat;
  return set_obj<ScalarAction>(args);
});

EXPORT(get_m_sq_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_sa)) {
    return NULL;
  }
  const ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  return py_convert(sa.m_sq);
});

EXPORT(get_lmbd_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_sa)) {
    return NULL;
  }
  const ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  return py_convert(sa.lmbd);
});

EXPORT(get_alpha_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_sa)) {
    return NULL;
  }
  const ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  return py_convert(sa.alpha);
});

EXPORT(action_node_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  PyObject* p_sf = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_sa, &p_sf)) {
    return NULL;
  }
  ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  const Field<double>& sf = py_convert_type<Field<double>>(p_sf);
  const double ret = sa.action_node(sf);
  return py_convert(ret);
});

EXPORT(hmc_m_hamilton_node_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  PyObject* p_sm = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_sa, &p_sm)) {
    return NULL;
  }
  ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  const Field<Complex>& sm = py_convert_type<Field<Complex>>(p_sm);
  const double ret = sa.hmc_m_hamilton_node(sm);
  return py_convert(ret);
});

EXPORT(hmc_set_force_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  PyObject* p_sm_force = NULL;
  PyObject* p_sf = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_sa, &p_sm_force, &p_sf)) {
    return NULL;
  }
  ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  Field<double>& sm_force = py_convert_type<Field<double>>(p_sm_force);
  const Field<double>& sf = py_convert_type<Field<double>>(p_sf);
  sa.hmc_set_force(sm_force, sf);
  Py_RETURN_NONE;
});

EXPORT(hmc_field_evolve_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  PyObject* p_sf = NULL;
  PyObject* p_sm = NULL;
  double step_size = 0.0;
  if (!PyArg_ParseTuple(args, "OOOd", &p_sa, &p_sf, &p_sm, &step_size)) {
    return NULL;
  }
  ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  Field<Complex>& sf = py_convert_type<Field<Complex>>(p_sf);
  const Field<Complex>& sm = py_convert_type<Field<Complex>>(p_sm);
  sa.hmc_field_evolve(sf, sm, step_size);
  Py_RETURN_NONE;
});

EXPORT(axial_current_node_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  PyObject* p_axial_cur = NULL;
  PyObject* p_sf = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_sa, &p_axial_cur, &p_sf)) {
    return NULL;
  }
  ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  Field<double>& axial_cur = py_convert_type<Field<double>>(p_axial_cur);
  const Field<double>& sf = py_convert_type<Field<double>>(p_sf);
  sa.axial_current_node(axial_cur, sf);
  Py_RETURN_NONE;
});

EXPORT(sum_sq_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  PyObject* p_f = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_sa, &p_f)) {
    return NULL;
  }
  ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  const Field<double>& f = py_convert_type<Field<double>>(p_f);
  const double ret = sa.sum_sq(f);
  return py_convert(ret);
});

EXPORT(hmc_set_rand_momentum_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  PyObject* p_sm = NULL;
  PyObject* p_rs = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_sa, &p_sm, &p_rs)) {
    return NULL;
  }
  ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  Field<Complex>& sm = py_convert_type<Field<Complex>>(p_sm);
  const RngState& rs = py_convert_type<RngState>(p_rs);
  sa.hmc_set_rand_momentum(sm, rs);
  Py_RETURN_NONE;
});
