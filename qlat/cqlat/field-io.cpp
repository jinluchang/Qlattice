#include "lib.h"

namespace qlat
{  //

template <class M>
PyObject* save_field_ctype(PyField& pf, const std::string& path,
                           const Coordinate& new_size_node)
{
  const Field<M>& f = *(Field<M>*)pf.cdata;
  const long ret = write_field(f, path, new_size_node);
  return py_convert(ret);
}

template <class M>
PyObject* load_field_ctype(PyField& pf, const std::string& path)
{
  Field<M>& f = *(Field<M>*)pf.cdata;
  const long ret = read_field(f, path);
  return py_convert(ret);
}

template <class M>
PyObject* convert_float_from_double_field_ctype(PyField& pf_new, PyField& pf)
{
  qassert(pf_new.ctype == "Float");
  Field<float>& f_new = *(Field<float>*)pf_new.cdata;
  const Field<M>& f = *(Field<M>*)pf.cdata;
  convert_field_float_from_double(f_new, f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* convert_double_from_float_field_ctype(PyField& pf_new, PyField& pf)
{
  qassert(pf.ctype == "Float");
  const Field<float>& f = *(Field<float>*)pf.cdata;
  Field<M>& f_new = *(Field<M>*)pf_new.cdata;
  convert_field_double_from_float(f_new, f);
  Py_RETURN_NONE;
}

template <class M>
PyObject* to_from_endianness_field_ctype(PyField& pf,
                                         const std::string& endianness_tag)
{
  Field<M>& f = *(Field<M>*)pf.cdata;
  if ("big_32" == endianness_tag) {
    to_from_big_endian_32(get_data(f));
  } else if ("big_64" == endianness_tag) {
    to_from_big_endian_64(get_data(f));
  } else if ("little_32" == endianness_tag) {
    to_from_little_endian_32(get_data(f));
  } else if ("little_64" == endianness_tag) {
    to_from_little_endian_64(get_data(f));
  } else {
    qassert(false);
  }
  Py_RETURN_NONE;
}

}  // namespace qlat

EXPORT(save_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_path = NULL;
  PyObject* p_new_size_node = NULL;
  if (!PyArg_ParseTuple(args, "OO|O", &p_field, &p_path, &p_new_size_node)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  std::string path;
  py_convert(path, p_path);
  Coordinate new_size_node;
  if (NULL != p_new_size_node) {
    py_convert(new_size_node, p_new_size_node);
  }
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, save_field_ctype, pf.ctype, pf, path, new_size_node);
  return p_ret;
})

EXPORT(load_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_path)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  std::string path;
  py_convert(path, p_path);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, load_field_ctype, pf.ctype, pf, path);
  return p_ret;
})

EXPORT(convert_float_from_double_field, {
  using namespace qlat;
  PyObject* p_field_new = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field_new, &p_field)) {
    return NULL;
  }
  PyField pf_new = py_convert_field(p_field_new);
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, convert_float_from_double_field_ctype, pf.ctype, pf_new,
                 pf);
  return p_ret;
})

EXPORT(convert_double_from_float_field, {
  using namespace qlat;
  PyObject* p_field_new = NULL;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field_new, &p_field)) {
    return NULL;
  }
  PyField pf_new = py_convert_field(p_field_new);
  PyField pf = py_convert_field(p_field);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, convert_double_from_float_field_ctype, pf_new.ctype,
                 pf_new, pf);
  return p_ret;
})

EXPORT(to_from_endianness_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_endianness_tag = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_field, &p_endianness_tag)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  std::string endianness_tag;
  py_convert(endianness_tag, p_endianness_tag);
  PyObject* p_ret = NULL;
  FIELD_DISPATCH(p_ret, to_from_endianness_field_ctype, pf.ctype, pf,
                 endianness_tag);
  return p_ret;
})
