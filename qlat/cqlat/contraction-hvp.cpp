#include "lib.h"

EXPORT(contract_chvp3_sfield, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  PyObject* p_prop1 = NULL;
  PyObject* p_prop2 = NULL;
  Long tslice_src = -1;
  if (!PyArg_ParseTuple(args, "OOOi", &p_ld, &p_prop1, &p_prop2, &tslice_src)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  const SelProp& prop1 = py_convert_type<SelProp>(p_prop1);
  const SelProp& prop2 = py_convert_type<SelProp>(p_prop2);
  QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_prop1, "fsel");
  QLAT_DIAGNOSTIC_POP;
  qassert(&fsel == &(py_convert_type<FieldSelection>(p_prop2, "fsel")));
  ld = contract_chvp3(prop1, prop2, tslice_src, fsel);
  Py_RETURN_NONE;
})
