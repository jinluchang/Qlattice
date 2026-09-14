#include "lib.h"

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

EXPORT(hmc_estimate_mass_scalar_action, {
  using namespace qlat;
  PyObject* p_sa = NULL;
  PyObject* p_masses = NULL;
  PyObject* p_field_ft = NULL;
  PyObject* p_force_ft = NULL;
  RealD phi0 = 0.0;
  if (!PyArg_ParseTuple(args, "OOOOd", &p_sa, &p_masses, &p_field_ft,
                        &p_force_ft, &phi0)) {
    return NULL;
  }
  ScalarAction& sa = py_convert_type<ScalarAction>(p_sa);
  Field<RealD>& masses = py_convert_type<Field<RealD>>(p_masses);
  const Field<ComplexD>& field_ft =
      py_convert_type<Field<ComplexD>>(p_field_ft);
  const Field<ComplexD>& force_ft =
      py_convert_type<Field<ComplexD>>(p_force_ft);
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
  Field<RealD>& sin_domega = py_convert_type<Field<RealD>>(p_sin_domega);
  sa.to_mass_factor(sin_domega);
  Py_RETURN_NONE;
})

