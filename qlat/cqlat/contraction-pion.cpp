#include "lib.h"

EXPORT(contract_pion_field, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  PyObject* p_prop = NULL;
  Int tslice_src = -1;
  if (!PyArg_ParseTuple(args, "OOi", &p_ld, &p_prop, &tslice_src)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  const Propagator4d& prop = py_convert_type<Propagator4d>(p_prop);
  ld = contract_pion(prop, tslice_src);
  Py_RETURN_NONE;
})

EXPORT(contract_pion_sfield, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  PyObject* p_prop = NULL;
  Int tslice_src = -1;
  if (!PyArg_ParseTuple(args, "OOi", &p_ld, &p_prop, &tslice_src)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  const SelProp& prop = py_convert_type<SelProp>(p_prop);
  QLAT_PUSH_DIAGNOSTIC_DISABLE_DANGLING_REF;
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_prop, "fsel");
  QLAT_DIAGNOSTIC_POP;
  ld = contract_pion(prop, tslice_src, fsel);
  Py_RETURN_NONE;
})
