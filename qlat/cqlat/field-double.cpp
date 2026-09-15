#include <qlat/field-double.h>

#include "lib.h"

namespace qlat
{  //

template <class M>
PyObject* multiply_double_field_ctype(PyObject* p_sf, PyObject* p_factor)
{
  Field<M>& sf = py_convert_type_field<M>(p_sf);
  Field<RealD>& factor = py_convert_type_field<RealD>(p_factor);
  multiply_double(sf, factor);
  Py_RETURN_NONE;
}

template <class M>
PyObject* invert_double_field_ctype(PyObject* p_sf)
{
  Field<M>& sf = py_convert_type_field<M>(p_sf);
  invert_double(sf);
  Py_RETURN_NONE;
}

}  // namespace qlat

EXPORT(multiply_double_field, {  // tested: cqlat-fields
  using namespace qlat;
  PyObject* p_sf = NULL;
  PyObject* p_factor = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_sf, &p_factor)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_sf);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, multiply_double_field_ctype, ctype, p_sf, p_factor);
  return p_ret;
})

EXPORT(invert_double_field, {  // tested: cqlat-fields
  using namespace qlat;
  PyObject* p_sf = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_sf)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_sf);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, invert_double_field_ctype, ctype, p_sf);
  return p_ret;
})
