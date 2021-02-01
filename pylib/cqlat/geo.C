#include "exceptions.h"
#include "convert.h"

EXPORT(mk_geo, {
  using namespace qlat;
  PyObject* p_total_site = NULL;
  int multiplicity = 1;
  if (!PyArg_ParseTuple(args, "O|i", &p_total_site, &multiplicity)) {
    return NULL;
  }
  Coordinate total_site;
  py_convert(total_site, p_total_site);
  Geometry* p_geo = new Geometry(total_site, multiplicity);
  return PyLong_FromVoidPtr(p_geo);
});

EXPORT(free_geo, {
  using namespace qlat;
  Geometry* p_geo = NULL;
  if (!PyArg_ParseTuple(args, "l", &p_geo)) {
    return NULL;
  }
  pqassert(p_geo != NULL);
  delete p_geo;
  Py_RETURN_NONE;
});
