#include "lib.h"

EXPORT(contract_chvp_16_field, {
  using namespace qlat;
  PyObject* p_chvp_16 = NULL;
  PyObject* p_prop1 = NULL;
  PyObject* p_prop2 = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_chvp_16, &p_prop1, &p_prop2)) {
    return NULL;
  }
  FieldM<Complex, 16>& chvp = py_convert_type<FieldM<Complex, 16> >(p_chvp_16);
  const Propagator4d& prop1 = py_convert_type<Propagator4d>(p_prop1);
  const Propagator4d& prop2 = py_convert_type<Propagator4d>(p_prop2);
  contract_chvp_16(chvp, prop1, prop2);
  Py_RETURN_NONE;
})
