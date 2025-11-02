#pragma once

#include <Python.h>
#include <qlat-utils/py_convert.h>
#include <qlat/qlat.h>

#define FIELD_DISPATCH(p_ret, fname, ctype, ...)                        \
  {                                                                     \
    if ("ColorMatrix" == (ctype)) {                                     \
      (p_ret) = fname<ColorMatrix>(__VA_ARGS__);                        \
    } else if ("WilsonMatrix" == (ctype)) {                             \
      (p_ret) = fname<WilsonMatrix>(__VA_ARGS__);                       \
    } else if ("NonRelWilsonMatrix" == (ctype)) {                       \
      (p_ret) = fname<NonRelWilsonMatrix>(__VA_ARGS__);                 \
    } else if ("IsospinMatrix" == (ctype)) {                            \
      (p_ret) = fname<IsospinMatrix>(__VA_ARGS__);                      \
    } else if ("SpinMatrix" == (ctype)) {                               \
      (p_ret) = fname<SpinMatrix>(__VA_ARGS__);                         \
    } else if ("WilsonVector" == (ctype)) {                             \
      (p_ret) = fname<WilsonVector>(__VA_ARGS__);                       \
    } else if ("ComplexD" == (ctype)) {                                 \
      (p_ret) = fname<ComplexD>(__VA_ARGS__);                           \
    } else if ("ComplexF" == (ctype)) {                                 \
      (p_ret) = fname<ComplexF>(__VA_ARGS__);                           \
    } else if ("RealD" == (ctype)) {                                    \
      (p_ret) = fname<RealD>(__VA_ARGS__);                              \
    } else if ("RealF" == (ctype)) {                                    \
      (p_ret) = fname<RealF>(__VA_ARGS__);                              \
    } else if ("Long" == (ctype)) {                                     \
      (p_ret) = fname<Long>(__VA_ARGS__);                               \
    } else if ("Int" == (ctype)) {                                      \
      (p_ret) = fname<Int>(__VA_ARGS__);                                \
    } else if ("Char" == (ctype)) {                                     \
      (p_ret) = fname<Char>(__VA_ARGS__);                               \
    } else if ("Int64t" == (ctype)) {                                   \
      (p_ret) = fname<int64_t>(__VA_ARGS__);                            \
    } else if ("Int32t" == (ctype)) {                                   \
      (p_ret) = fname<int32_t>(__VA_ARGS__);                            \
    } else if ("Int8t" == (ctype)) {                                    \
      (p_ret) = fname<int8_t>(__VA_ARGS__);                             \
    } else {                                                            \
      qerr(qlat::ssprintf("%s %s='%s' does not exist.", #fname, #ctype, \
                          (ctype).c_str()));                            \
      (p_ret) = NULL;                                                   \
    }                                                                   \
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
inline bool check_ctype_name<IsospinMatrix>(const std::string& ctype)
{
  return "IsospinMatrix" == ctype;
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
    qassert(DIMN == PyList_Size(in));
    for (size_t i = 0; i < out.size(); i++) {
      py_convert(out[i], PyList_GetItem(in, i));
    }
  } else if (PyTuple_Check(in)) {
    qassert(DIMN == PyTuple_Size(in));
    for (size_t i = 0; i < out.size(); i++) {
      py_convert(out[i], PyTuple_GetItem(in, i));
    }
  } else {
    qassert(false);
  }
}

inline void py_convert(CoordinateD& out, PyObject* in)
{
  if (PyList_Check(in)) {
    qassert(DIMN == PyList_Size(in));
    for (size_t i = 0; i < out.size(); i++) {
      py_convert(out[i], PyList_GetItem(in, i));
    }
  } else if (PyTuple_Check(in)) {
    qassert(DIMN == PyTuple_Size(in));
    for (size_t i = 0; i < out.size(); i++) {
      py_convert(out[i], PyTuple_GetItem(in, i));
    }
  } else {
    qassert(false);
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
  out.ctype = py_get_ctype(in);
  out.cdata = (void*)py_convert_data<Long>(in, "cdata");
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
  qassert(check_ctype_name<M>(ctype));
  Field<M>& f = py_convert_type<Field<M>>(in);
  return f;
}

template <class M, Int multiplicity>
FieldM<M, multiplicity>& py_convert_type_field(PyObject* in)
// py_convert_type<FieldM<M, multiplicity> >(in);
{
  const std::string ctype = py_get_ctype(in);
  qassert(check_ctype_name<M>(ctype));
  FieldM<M, multiplicity>& f = py_convert_type<FieldM<M, multiplicity>>(in);
  if (is_initialized(f)) {
    qassert(multiplicity == f.multiplicity);
  }
  return f;
}

template <class M>
SelectedField<M>& py_convert_type_sfield(PyObject* in)
// interface
// py_convert_type<SelectedField<M> >(in);
{
  const std::string ctype = py_get_ctype(in);
  qassert(check_ctype_name<M>(ctype));
  SelectedField<M>& f = py_convert_type<SelectedField<M>>(in);
  return f;
}

template <class M>
SelectedPoints<M>& py_convert_type_spoints(PyObject* in)
// interface
// py_convert_type<SelectedPoints<M> >(in);
{
  const std::string ctype = py_get_ctype(in);
  qassert(check_ctype_name<M>(ctype));
  SelectedPoints<M>& f = py_convert_type<SelectedPoints<M>>(in);
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
  for (Long i = 0; i < (Long)coor.size(); i++) {
    PyList_SetItem(ret, i, PyLong_FromLong(coor[i]));
  }
  return ret;
}

inline PyObject* py_convert(const ColorMatrix& x)
{
  return py_convert(get_data_complex_d(x));
}

inline PyObject* py_convert(const WilsonMatrix& x)
{
  return py_convert(get_data_complex_d(x));
}

inline PyObject* py_convert(const NonRelWilsonMatrix& x)
{
  return py_convert(get_data_complex_d(x));
}

inline PyObject* py_convert(const IsospinMatrix& x)
{
  return py_convert(get_data_complex_d(x));
}

inline PyObject* py_convert(const SpinMatrix& x)
{
  return py_convert(get_data_complex_d(x));
}

inline PyObject* py_convert(const WilsonVector& x)
{
  return py_convert(get_data_complex_d(x));
}

}  // namespace qlat
