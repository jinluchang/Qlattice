#include "lib.h"

EXPORT(bcast_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_ld)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  bcast(ld);
  Py_RETURN_NONE;
})

EXPORT(glb_sum_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_ld)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  qlat::glb_sum_lat_data(ld);
  Py_RETURN_NONE;
})
