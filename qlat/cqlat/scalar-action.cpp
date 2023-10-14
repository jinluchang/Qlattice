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
})

EXPORT(free_scalar_action, {
  using namespace qlat;
  return free_obj<ScalarAction>(args);
})

EXPORT(set_scalar_action, {
  using namespace qlat;
  return set_obj<ScalarAction>(args);
})

EXPORT(get_m_sq_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_sa)) {
    return NULL;
  }
  const ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  return py_convert(sa.m_sq);
})

EXPORT(get_lmbd_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_sa)) {
    return NULL;
  }
  const ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  return py_convert(sa.lmbd);
})

EXPORT(get_alpha_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_sa)) {
    return NULL;
  }
  const ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  return py_convert(sa.alpha);
})

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
})

EXPORT(hmc_estimate_mass_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  PyObject* p_masses = NULL;
  PyObject* p_field_ft = NULL;
  PyObject* p_force_ft = NULL;
  double phi0 = 0.0;
  if (!PyArg_ParseTuple(args, "OOOOd", &p_sa, &p_masses, &p_field_ft, &p_force_ft, &phi0)) {
    return NULL;
  }
  ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  Field<double>& masses = py_convert_type<Field<double>>(p_masses);
  const Field<Complex>& field_ft = py_convert_type<Field<Complex>>(p_field_ft);
  const Field<Complex>& force_ft = py_convert_type<Field<Complex>>(p_force_ft);
  sa.hmc_estimate_mass(masses, field_ft, force_ft, phi0);
  Py_RETURN_NONE;
})

EXPORT(to_mass_factor_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  PyObject* p_sin_domega = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_sa, &p_sin_domega)) {
    return NULL;
  }
  ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  Field<double>& sin_domega = py_convert_type<Field<double>>(p_sin_domega);
  sa.to_mass_factor(sin_domega);
  Py_RETURN_NONE;
})

EXPORT(hmc_m_hamilton_node_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  PyObject* p_sm = NULL;
  PyObject* p_masses = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_sa, &p_sm, &p_masses)) {
    return NULL;
  }
  ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  const Field<Complex>& sm = py_convert_type<Field<Complex>>(p_sm);
  const Field<double>& masses = py_convert_type<Field<double>>(p_masses);
  const double ret = sa.hmc_m_hamilton_node(sm, masses);
  return py_convert(ret);
})

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
})

EXPORT(hmc_field_evolve_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  PyObject* p_sf = NULL;
  PyObject* p_sm = NULL;
  PyObject* p_masses = NULL;
  double step_size = 0.0;
  if (!PyArg_ParseTuple(args, "OOOOd", &p_sa, &p_sf, &p_sm, &p_masses, &step_size)) {
    return NULL;
  }
  ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  Field<Complex>& sf = py_convert_type<Field<Complex>>(p_sf);
  const Field<Complex>& sm = py_convert_type<Field<Complex>>(p_sm);
  const Field<double>& masses = py_convert_type<Field<double>>(p_masses);
  sa.hmc_field_evolve(sf, sm, masses, step_size);
  Py_RETURN_NONE;
})

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
})

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
})

EXPORT(hmc_set_rand_momentum_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  PyObject* p_sm = NULL;
  PyObject* p_masses = NULL;
  PyObject* p_rs = NULL;
  if (!PyArg_ParseTuple(args, "OOOO", &p_sa, &p_sm, &p_masses, &p_rs)) {
    return NULL;
  }
  ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  Field<Complex>& sm = py_convert_type<Field<Complex>>(p_sm);
  Field<double>& masses = py_convert_type<Field<double>>(p_masses);
  const RngState& rs = py_convert_type<RngState>(p_rs);
  sa.hmc_set_rand_momentum(sm, masses, rs);
  Py_RETURN_NONE;
})

EXPORT(hmc_predict_field_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  PyObject* p_sf_ft = NULL;
  PyObject* p_sm_ft = NULL;
  PyObject* p_masses = NULL;
  double vev_sigma = 0.0;
  if (!PyArg_ParseTuple(args, "OOOOd", &p_sa, &p_sf_ft, &p_sm_ft, &p_masses, &vev_sigma)) {
    return NULL;
  }
  ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  Field<Complex>& sf_ft = py_convert_type<Field<Complex>>(p_sf_ft);
  const Field<Complex>& sm_ft = py_convert_type<Field<Complex>>(p_sm_ft);
  const Field<double>& masses = py_convert_type<Field<double>>(p_masses);
  sa.hmc_predict_field(sf_ft, sm_ft, masses, vev_sigma);
  Py_RETURN_NONE;
})

EXPORT(get_polar_field_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  PyObject* p_pf = NULL;
  PyObject* p_sf = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_sa, &p_pf, &p_sf)) {
    return NULL;
  }
  ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  Field<double>& pf = py_convert_type<Field<double>>(p_pf);
  const Field<double>& sf = py_convert_type<Field<double>>(p_sf);
  sa.get_polar_field(pf, sf);
  Py_RETURN_NONE;
})
