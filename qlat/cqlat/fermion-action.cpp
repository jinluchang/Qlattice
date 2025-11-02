#include "lib.h"

EXPORT(mk_fermion_action_mobius, {
  using namespace qlat;
  RealD mass = 0.0;
  Int ls = 0;
  RealD m5 = 0.0;
  RealD mobius_scale = 0.0;
  if (!PyArg_ParseTuple(args, "didd", &mass, &ls, &m5, &mobius_scale)) {
    return NULL;
  }
  FermionAction* pfa = new FermionAction(mass, ls, m5, mobius_scale, true, true);
  return py_convert((void*)pfa);
})

EXPORT(mk_fermion_action_zmobius, {
  using namespace qlat;
  RealD mass = 0.0;
  RealD m5 = 0.0;
  PyObject* p_omega = NULL;
  if (!PyArg_ParseTuple(args, "ddO", &mass, &m5, &p_omega)) {
    return NULL;
  }
  std::vector<ComplexD> omega;
  py_convert(omega, p_omega);
  const Int ls = omega.size();
  FermionAction* pfa = new FermionAction(mass, ls, m5, 0.0, true, true);
  FermionAction& fa = *pfa;
  qassert(fa.bs.size() == omega.size());
  qassert(fa.cs.size() == omega.size());
  for (Int i = 0; i < (int)omega.size(); ++i) {
    fa.bs[i] = 0.5 * (1.0 / omega[i] + 1.0);
    fa.cs[i] = fa.bs[i] - 1.0;
  }
  return py_convert((void*)pfa);
})

EXPORT(free_fermion_action, {
  using namespace qlat;
  return free_obj<FermionAction>(args);
})

EXPORT(set_fermion_action, {
  using namespace qlat;
  return set_obj<FermionAction>(args);
})

EXPORT(get_mass_fermion_action, {
  using namespace qlat;
  PyObject* p_fa = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_fa)) {
    return NULL;
  }
  const FermionAction& fa = py_convert_type<FermionAction>(p_fa);
  return py_convert(fa.mass);
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

EXPORT(get_m5_fermion_action, {
  using namespace qlat;
  PyObject* p_fa = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_fa)) {
    return NULL;
  }
  const FermionAction& fa = py_convert_type<FermionAction>(p_fa);
  return py_convert(fa.m5);
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
