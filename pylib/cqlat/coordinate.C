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
});

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
});
