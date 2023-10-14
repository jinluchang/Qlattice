#include "lib.h"

EXPORT(contract_pion_field, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  PyObject* p_prop = NULL;
  int tslice_src = -1;
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
  int tslice_src = -1;
  if (!PyArg_ParseTuple(args, "OOi", &p_ld, &p_prop, &tslice_src)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  const SelProp& prop = py_convert_type<SelProp>(p_prop);
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_prop, "fsel");
  ld = contract_pion(prop, tslice_src, fsel);
  Py_RETURN_NONE;
})
