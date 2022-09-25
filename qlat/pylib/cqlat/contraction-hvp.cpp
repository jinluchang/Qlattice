#include "lib.h"

EXPORT(contract_chvp3_sfield, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  PyObject* p_prop1 = NULL;
  PyObject* p_prop2 = NULL;
  long tslice_src = -1;
  if (!PyArg_ParseTuple(args, "OOOi", &p_ld, &p_prop1, &p_prop2, &tslice_src)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  const SelProp& prop1 = py_convert_type<SelProp>(p_prop1);
  const SelProp& prop2 = py_convert_type<SelProp>(p_prop2);
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_prop1, "fsel");
  pqassert(&fsel == &(py_convert_type<FieldSelection>(p_prop2, "fsel")));
  ld = contract_chvp3(prop1, prop2, tslice_src, fsel);
  Py_RETURN_NONE;
})
