#include "lib.h"

namespace qlat
{  //

template <class M>
PyObject* set_checkers_double_field_ctype(PyObject* p_field)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  set_checkers_double(f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_complex_from_double_field_ctype(PyObject* p_field, PyObject* p_sf)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  Field<double>& sf = py_convert_type_field<double>(p_sf);
  set_complex_from_double(f, sf);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_double_from_complex_field_ctype(PyObject* p_field, PyObject* p_cf)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  Field<Complex>& cf = py_convert_type_field<Complex>(p_cf);
  set_double_from_complex(f, cf);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_abs_from_complex_field_ctype(PyObject* p_field, PyObject* p_cf)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  Field<Complex>& cf = py_convert_type_field<Complex>(p_cf);
  set_abs_from_complex(f, cf);
  Py_RETURN_NONE;
}

template <class M>
PyObject* set_ratio_double_field_ctype(PyObject* p_field, PyObject* p_sf1, PyObject* p_sf2)
{
  Field<M>& f = py_convert_type_field<M>(p_field);
  Field<double>& sf1 = py_convert_type_field<double>(p_sf1);
  Field<double>& sf2 = py_convert_type_field<double>(p_sf2);
  set_ratio_double(f, sf1, sf2);
  Py_RETURN_NONE;
}

template <class M>
PyObject* less_than_double_field_ctype(PyObject* p_sf1, PyObject* p_sf2, PyObject* p_mask)
{
  Field<M>& sf1 = py_convert_type_field<M>(p_sf1);
  Field<double>& sf2 = py_convert_type_field<double>(p_sf2);
  Field<double>& mask = py_convert_type_field<double>(p_mask);
  less_than_double(sf1, sf2, mask);
  Py_RETURN_NONE;
}

template <class M>
PyObject* invert_double_field_ctype(PyObject* p_sf)
{
  Field<M>& sf = py_convert_type_field<M>(p_sf);
  invert_double(sf);
  Py_RETURN_NONE;
}

} // namespace qlat

EXPORT(set_checkers_double_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_checkers_double_field_ctype, ctype, p_field);
  return p_ret;
})

EXPORT(set_complex_from_double_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_sf = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_sf)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_complex_from_double_field_ctype, ctype, p_field, p_sf);
  return p_ret;
})

EXPORT(set_double_from_complex_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_cf = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_cf)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_double_from_complex_field_ctype, ctype, p_field, p_cf);
  return p_ret;
})

EXPORT(set_abs_from_complex_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_cf = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_cf)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_abs_from_complex_field_ctype, ctype, p_field, p_cf);
  return p_ret;
})

EXPORT(set_ratio_double_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_sf1 = NULL;
  PyObject* p_sf2 = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_field, &p_sf1, &p_sf2)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, set_ratio_double_field_ctype, ctype, p_field, p_sf1, p_sf2);
  return p_ret;
})

EXPORT(less_than_double_field, {
  using namespace qlat;
  PyObject* p_sf1 = NULL;
  PyObject* p_sf2 = NULL;
  PyObject* p_mask = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_sf1, &p_sf2, &p_mask)) {
    return NULL;
  }
  const std::string ctype = py_get_ctype(p_sf1);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, less_than_double_field_ctype, ctype, p_sf1, p_sf2, p_mask);
  return p_ret;
})

EXPORT(invert_double_field, {
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
