#include "lib.h"

EXPORT(mk_scalar_action, {
  using namespace qlat;
  double m_sq = 0.0;
  double lmbd = 0.0;
  if (!PyArg_ParseTuple(args, "d|d", &m_sq, &lmbd)) {
    return NULL;
  }
  ScalarAction* psa = new ScalarAction(m_sq, lmbd);
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
