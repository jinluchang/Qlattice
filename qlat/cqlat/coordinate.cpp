#include "lib.h"

EXPORT(mod_coordinate, {
  using namespace qlat;
  PyObject* p_xg = NULL;
  PyObject* p_size = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_xg, &p_size)) {
    return NULL;
  }
  const Coordinate xg = py_convert_data<Coordinate>(p_xg);
  const Coordinate size = py_convert_data<Coordinate>(p_size);
  return py_convert(mod(xg, size));
})

EXPORT(smod_coordinate, {
  using namespace qlat;
  PyObject* p_xg = NULL;
  PyObject* p_size = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_xg, &p_size)) {
    return NULL;
  }
  const Coordinate xg = py_convert_data<Coordinate>(p_xg);
  const Coordinate size = py_convert_data<Coordinate>(p_size);
  return py_convert(smod(xg, size));
})

EXPORT(c_rand_gen, {
  using namespace qlat;
  PyObject* p_rng = NULL;
  PyObject* p_size = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_rng, &p_size)) {
    return NULL;
  }
  RngState& rng = py_convert_type<RngState>(p_rng);
  Coordinate size;
  py_convert(size, p_size);
  const Coordinate x = c_rand_gen(rng, size);
  return py_convert(x);
})
