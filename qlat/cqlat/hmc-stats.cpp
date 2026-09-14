#include "lib.h"

EXPORT(get_gm_force_magnitudes, {
  using namespace qlat;
  PyObject* p_gm_force = NULL;
  Int n_elems = 0;
  if (!PyArg_ParseTuple(args, "Oi", &p_gm_force, &n_elems)) {
    return NULL;
  }
  const GaugeMomentum& gm_force = py_convert_type<GaugeMomentum>(p_gm_force);
  std::vector<RealD> ret = get_gm_force_magnitudes(gm_force, n_elems);
  return py_convert(ret);
})

