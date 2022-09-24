#pragma once

#include <Python.h>
#include <qlat-utils/py_convert.h>
#include <qlat/qlat.h>

#define FIELD_DISPATCH(p_ret, fname, ctype, ...)                            \
  {                                                                         \
    if ("ColorMatrix" == (ctype)) {                                         \
      (p_ret) = fname<ColorMatrix>(__VA_ARGS__);                            \
    } else if ("WilsonMatrix" == (ctype)) {                                 \
      (p_ret) = fname<WilsonMatrix>(__VA_ARGS__);                           \
    } else if ("NonRelWilsonMatrix" == (ctype)) {                           \
      (p_ret) = fname<NonRelWilsonMatrix>(__VA_ARGS__);                     \
    } else if ("SpinMatrix" == (ctype)) {                                   \
      (p_ret) = fname<SpinMatrix>(__VA_ARGS__);                             \
    } else if ("WilsonVector" == (ctype)) {                                 \
      (p_ret) = fname<WilsonVector>(__VA_ARGS__);                           \
    } else if ("Complex" == (ctype)) {                                      \
      (p_ret) = fname<Complex>(__VA_ARGS__);                                \
    } else if ("ComplexF" == (ctype)) {                                     \
      (p_ret) = fname<ComplexF>(__VA_ARGS__);                               \
    } else if ("double" == (ctype)) {                                       \
      (p_ret) = fname<double>(__VA_ARGS__);                                 \
    } else if ("float" == (ctype)) {                                        \
      (p_ret) = fname<float>(__VA_ARGS__);                                  \
    } else if ("long" == (ctype)) {                                         \
      (p_ret) = fname<long>(__VA_ARGS__);                                   \
    } else if ("int64_t" == (ctype)) {                                      \
      (p_ret) = fname<int64_t>(__VA_ARGS__);                                \
    } else if ("char" == (ctype)) {                                         \
      (p_ret) = fname<char>(__VA_ARGS__);                                   \
    } else if ("int8_t" == (ctype)) {                                       \
      (p_ret) = fname<int8_t>(__VA_ARGS__);                                 \
    } else {                                                                \
      pqerr("%s %s='%s' does not exist.", #fname, #ctype, (ctype).c_str()); \
      (p_ret) = NULL;                                                       \
    }                                                                       \
  }

namespace qlat
{  //

template <>
inline bool check_ctype_name<ColorMatrix>(const std::string& ctype)
{
  return "ColorMatrix" == ctype;
}

template <>
inline bool check_ctype_name<WilsonMatrix>(const std::string& ctype)
{
  return "WilsonMatrix" == ctype;
}

template <>
inline bool check_ctype_name<NonRelWilsonMatrix>(const std::string& ctype)
{
  return "NonRelWilsonMatrix" == ctype;
}

template <>
inline bool check_ctype_name<SpinMatrix>(const std::string& ctype)
{
  return "SpinMatrix" == ctype;
}

template <>
inline bool check_ctype_name<WilsonVector>(const std::string& ctype)
{
  return "WilsonVector" == ctype;
}

inline void py_convert(Coordinate& out, PyObject* in)
{
  if (PyList_Check(in)) {
    pqassert(DIMN == PyList_Size(in));
    for (size_t i = 0; i < out.size(); i++) {
      py_convert(out[i], PyList_GetItem(in, i));
    }
  } else if (PyTuple_Check(in)) {
    pqassert(DIMN == PyTuple_Size(in));
    for (size_t i = 0; i < out.size(); i++) {
      py_convert(out[i], PyTuple_GetItem(in, i));
    }
  } else {
    pqassert(false);
  }
}

inline void py_convert(CoordinateD& out, PyObject* in)
{
  if (PyList_Check(in)) {
    pqassert(DIMN == PyList_Size(in));
    for (size_t i = 0; i < out.size(); i++) {
      py_convert(out[i], PyList_GetItem(in, i));
    }
  } else if (PyTuple_Check(in)) {
    pqassert(DIMN == PyTuple_Size(in));
    for (size_t i = 0; i < out.size(); i++) {
      py_convert(out[i], PyTuple_GetItem(in, i));
    }
  } else {
    pqassert(false);
  }
}

struct PyField {
  // deprecated
  std::string ctype;
  void* cdata;
};

inline void py_convert(PyField& out, PyObject* in)
// deprecated
{
  PyObject* p_ctype = PyObject_GetAttrString(in, "ctype");
  PyObject* p_cdata = PyObject_GetAttrString(in, "cdata");
  py_convert(out.ctype, p_ctype);
  py_convert((long&)out.cdata, p_cdata);
  Py_DECREF(p_ctype);
  Py_DECREF(p_cdata);
}

inline PyField py_convert_field(PyObject* in)
// deprecated
{
  return py_convert_data<PyField>(in);
}

template <class M>
Field<M>& py_convert_type_field(PyObject* in)
// py_convert_type<Field<M> >(in);
{
  const std::string ctype = py_get_ctype(in);
  pqassert(check_ctype_name<M>(ctype));
  Field<M>& f = py_convert_type<Field<M> >(in);
  return f;
}

template <class M, int multiplicity>
FieldM<M, multiplicity>& py_convert_type_field(PyObject* in)
// py_convert_type<FieldM<M, multiplicity> >(in);
{
  const std::string ctype = py_get_ctype(in);
  pqassert(check_ctype_name<M>(ctype));
  FieldM<M, multiplicity>& f = py_convert_type<FieldM<M, multiplicity> >(in);
  if (is_initialized(f)) {
    pqassert(multiplicity == f.geo().multiplicity);
  }
  return f;
}

template <class M>
SelectedField<M>& py_convert_type_sfield(PyObject* in)
// interface
// py_convert_type<SelectedField<M> >(in);
{
  const std::string ctype = py_get_ctype(in);
  pqassert(check_ctype_name<M>(ctype));
  SelectedField<M>& f = py_convert_type<SelectedField<M> >(in);
  return f;
}

template <class M>
SelectedPoints<M>& py_convert_type_spoints(PyObject* in)
// interface
// py_convert_type<SelectedPoints<M> >(in);
{
  const std::string ctype = py_get_ctype(in);
  pqassert(check_ctype_name<M>(ctype));
  SelectedPoints<M>& f = py_convert_type<SelectedPoints<M> >(in);
  return f;
}

template <>
inline CommMarks& py_convert_type(PyObject* in)
{
  Field<int8_t>& f = py_convert_type_field<int8_t>(in);
  CommMarks& ret = static_cast<CommMarks&>(f);
  return ret;
}

template <>
inline Propagator4d& py_convert_type(PyObject* in)
{
  FieldM<WilsonMatrix, 1>& f = py_convert_type_field<WilsonMatrix, 1>(in);
  Propagator4d& ret = static_cast<Propagator4d&>(f);
  return ret;
}

template <>
inline GaugeField& py_convert_type(PyObject* in)
{
  FieldM<ColorMatrix, 4>& f = py_convert_type_field<ColorMatrix, 4>(in);
  GaugeField& ret = static_cast<GaugeField&>(f);
  return ret;
}

// -------------------------------------------------------------

inline PyObject* py_convert(const Coordinate& coor)
{
  PyObject* ret = PyList_New(coor.size());
  for (long i = 0; i < (long)coor.size(); i++) {
    PyList_SetItem(ret, i, PyLong_FromLong(coor[i]));
  }
  return ret;
}

inline PyObject* py_convert(const ColorMatrix& x) {
  return py_convert(get_data_complex(get_data_one_elem(x)));
}

inline PyObject* py_convert(const WilsonMatrix& x) {
  return py_convert(get_data_complex(get_data_one_elem(x)));
}

inline PyObject* py_convert(const NonRelWilsonMatrix& x) {
  return py_convert(get_data_complex(get_data_one_elem(x)));
}

inline PyObject* py_convert(const SpinMatrix& x) {
  return py_convert(get_data_complex(get_data_one_elem(x)));
}

inline PyObject* py_convert(const WilsonVector& x) {
  return py_convert(get_data_complex(get_data_one_elem(x)));
}

}  // namespace qlat
