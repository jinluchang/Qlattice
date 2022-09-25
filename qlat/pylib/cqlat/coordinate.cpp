#include "lib.h"

EXPORT(index_from_coordinate, {
  using namespace qlat;
  PyObject* p_coor = NULL;
  PyObject* p_size = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_coor, &p_size)) {
    return NULL;
  }
  Coordinate coor, size;
  py_convert(coor, p_coor);
  py_convert(size, p_size);
  return py_convert(index_from_coordinate(coor, size));
})

EXPORT(coordinate_from_index, {
  using namespace qlat;
  long index = -1;
  PyObject* p_size = NULL;
  if (!PyArg_ParseTuple(args, "lO", &index, &p_size)) {
    return NULL;
  }
  Coordinate size;
  py_convert(size, p_size);
  return py_convert(coordinate_from_index(index, size));
})

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
