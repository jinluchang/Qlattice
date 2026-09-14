#include "lib.h"

EXPORT(free_fermion_action, {
  using namespace qlat;
  return free_obj<FermionAction>(args);
})

EXPORT(set_fermion_action, {
  using namespace qlat;
  return set_obj<FermionAction>(args);
})

EXPORT(get_ls_fermion_action, {
  using namespace qlat;
  PyObject* p_fa = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_fa)) {
    return NULL;
  }
  const FermionAction& fa = py_convert_type<FermionAction>(p_fa);
  return py_convert(fa.ls);
})

EXPORT(get_omega_fermion_action, {
  using namespace qlat;
  PyObject* p_fa = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_fa)) {
    return NULL;
  }
  const FermionAction& fa = py_convert_type<FermionAction>(p_fa);
  if (fa.is_using_zmobius) {
    std::vector<ComplexD> omega(fa.bs.size());
    for (Int i = 0; i < (int)omega.size(); ++i) {
      omega[i] = 1.0 / (fa.bs[i] + fa.cs[i]);
    }
    return py_convert(omega);
  } else {
    Py_RETURN_NONE;
  }
})

EXPORT(get_mobius_scale_fermion_action, {
  using namespace qlat;
  PyObject* p_fa = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_fa)) {
    return NULL;
  }
  const FermionAction& fa = py_convert_type<FermionAction>(p_fa);
  if (fa.is_using_zmobius) {
    qassert(fa.mobius_scale == 0.0);
  } else {
    qassert(fa.mobius_scale != 0.0);
  }
  return py_convert(fa.mobius_scale);
})
