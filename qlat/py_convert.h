#pragma once

#include <Python.h>
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

template <class M>
bool check_ctype_name(const std::string& ctype)
{
  return false;
}

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

template <>
inline bool check_ctype_name<Complex>(const std::string& ctype)
{
  return "Complex" == ctype;
}

template <>
inline bool check_ctype_name<double>(const std::string& ctype)
{
  return "double" == ctype;
}

template <>
inline bool check_ctype_name<float>(const std::string& ctype)
{
  return "float" == ctype;
}

template <>
inline bool check_ctype_name<int64_t>(const std::string& ctype)
{
  return "int64_t" == ctype or "long" == ctype;
}

template <>
inline bool check_ctype_name<char>(const std::string& ctype)
{
  return "char" == ctype;
}

template <>
inline bool check_ctype_name<int8_t>(const std::string& ctype)
{
  return "int8_t" == ctype;
}

inline void py_convert(PyObject*& out, PyObject* in) { out = in; }

inline void py_convert(int& out, PyObject* in)
{
  pqassert(PyLong_Check(in));
  out = PyLong_AsLong(in);
}

inline void py_convert(long& out, PyObject* in)
{
  pqassert(PyLong_Check(in));
  out = PyLong_AsLong(in);
}

inline void py_convert(double& out, PyObject* in)
{
  if (PyFloat_Check(in)) {
    out = PyFloat_AsDouble(in);
  } else if (PyLong_Check(in)) {
    out = PyLong_AsLong(in);
  } else {
    pqassert(false);
  }
}

inline void py_convert(Complex& out, PyObject* in)
{
  if (PyLong_Check(in)) {
    out = PyLong_AsLong(in);
  } else if (PyFloat_Check(in)) {
    out = PyFloat_AsDouble(in);
  } else {
    pqassert(PyComplex_Check(in));
    Py_complex& py_out = (Py_complex&)out;
    py_out = PyComplex_AsCComplex(in);
  }
}

inline void py_convert(bool& out, PyObject* in)
{
  pqassert(PyBool_Check(in));
  out = in == Py_True;
}

inline void py_convert(std::string& s, PyObject* in)
{
  if (PyType_Check(in)) {
    s = ((PyTypeObject*)in)->tp_name;
  } else if (PyBytes_Check(in)) {
    s = PyBytes_AsString(in);
  } else if (PyUnicode_Check(in)) {
    PyObject* temp = PyUnicode_AsEncodedString(in, "UTF-8", "strict");
    pqassert(temp);
    s = PyBytes_AS_STRING(temp);
    Py_DECREF(temp);
  } else {
    pqassert(false);
  }
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

template <class M>
void py_convert(std::vector<M>& out, PyObject* in)
{
  if (PyList_Check(in)) {
    out.resize(PyList_Size(in));
    for (size_t i = 0; i < out.size(); i++) {
      py_convert(out[i], PyList_GetItem(in, i));
    }
  } else if (PyTuple_Check(in)) {
    out.resize(PyTuple_Size(in));
    for (size_t i = 0; i < out.size(); i++) {
      py_convert(out[i], PyTuple_GetItem(in, i));
    }
  } else {
    pqassert(false);
  }
}

template <>
inline void py_convert<bool>(std::vector<bool>& out, PyObject* in)
{
  if (PyList_Check(in)) {
    out.resize(PyList_Size(in));
    for (size_t i = 0; i < out.size(); i++) {
      bool v;
      py_convert(v, PyList_GetItem(in, i));
      out[i] = v;
    }
  } else if (PyTuple_Check(in)) {
    out.resize(PyTuple_Size(in));
    for (size_t i = 0; i < out.size(); i++) {
      bool v;
      py_convert(v, PyTuple_GetItem(in, i));
      out[i] = v;
    }
  } else {
    pqassert(false);
  }
}

template <class M>
void py_convert(Vector<M> out, PyObject* in)
{
  if (PyList_Check(in)) {
    pqassert(out.size() == PyList_Size(in));
    for (long i = 0; i < out.size(); i++) {
      py_convert(out[i], PyList_GetItem(in, i));
    }
  } else if (PyTuple_Check(in)) {
    pqassert(out.size() == PyTuple_Size(in));
    for (long i = 0; i < out.size(); i++) {
      py_convert(out[i], PyTuple_GetItem(in, i));
    }
  } else {
    pqassert(false);
  }
}

template <class T>
T py_convert_data(PyObject* in)
// examples:
// py_convert_data<int>(in)
// py_convert_data<long>(in)
// py_convert_data<double>(in)
// py_convert_data<Complex>(in)
// py_convert_data<bool>(in)
// py_convert_data<std::string>(in)
// py_convert_data<Coordinate>(in)
// py_convert_data<CoordinateD>(in)
// py_convert_data<std::vector<M> >(in)
{
  T x;
  py_convert(x, in);
  return x;
}

template <class T>
T py_convert_data(PyObject* in, const std::string& attr)
// py_convert_data<std::string>(in, "ctype")
{
  PyObject* p_obj = PyObject_GetAttrString(in, attr.c_str());
  T x = py_convert_data<T>(p_obj);
  Py_DECREF(p_obj);
  return x;
}

inline std::string py_get_ctype(PyObject* in)
{
  return py_convert_data<std::string>(in, "ctype");
}

struct PyField {
  std::string ctype;
  void* cdata;
};

inline void py_convert(PyField& out, PyObject* in)
{
  PyObject* p_ctype = PyObject_GetAttrString(in, "ctype");
  PyObject* p_cdata = PyObject_GetAttrString(in, "cdata");
  py_convert(out.ctype, p_ctype);
  py_convert((long&)out.cdata, p_cdata);
  Py_DECREF(p_ctype);
  Py_DECREF(p_cdata);
}

inline PyField py_convert_field(PyObject* in)
{
  return py_convert_data<PyField>(in);
}

template <class T>
T& py_convert_type(PyObject* in)
// use cdata property of PyObject* in as pointer
// examples:
// py_convert_type<Geometry>(in);
// py_convert_type<LatData>(in);
// py_convert_type<RngState>(in);
// py_convert_type<PointSelection>(in);
// py_convert_type<CommPlan>(in);
// specifications:
// py_convert_type<Propagator4d>(in);
// py_convert_type<GaugeField>(in);
// py_convert_type<CommMarks>(in);
{
  PyObject* p_cdata = PyObject_GetAttrString(in, "cdata");
  T* out;
  py_convert((long&)out, p_cdata);
  Py_DECREF(p_cdata);
  return *out;
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
// py_convert_type<SelectedField<M> >(in);
{
  const std::string ctype = py_get_ctype(in);
  pqassert(check_ctype_name<M>(ctype));
  SelectedField<M>& f = py_convert_type<SelectedField<M> >(in);
  return f;
}

template <class M>
SelectedPoints<M>& py_convert_type_spoints(PyObject* in)
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

inline PyObject* py_convert(const char& x) { return PyLong_FromLong((long)x); }

inline PyObject* py_convert(const bool& x) { return PyBool_FromLong((long)x); }

inline PyObject* py_convert(const int& x) { return PyLong_FromLong((long)x); }

inline PyObject* py_convert(const long& x) { return PyLong_FromLong(x); }

inline PyObject* py_convert(const long long& x)
{
  return PyLong_FromLongLong(x);
}

inline PyObject* py_convert(const unsigned int& x)
{
  return PyLong_FromUnsignedLong((unsigned long)x);
}

inline PyObject* py_convert(const unsigned long& x)
{
  return PyLong_FromUnsignedLong(x);
}

inline PyObject* py_convert(const unsigned long long& x)
{
  return PyLong_FromUnsignedLongLong(x);
}

inline PyObject* py_convert(const float& x)
{
  return PyFloat_FromDouble((double)x);
}

inline PyObject* py_convert(const double& x) { return PyFloat_FromDouble(x); }

inline PyObject* py_convert(void* x) { return PyLong_FromVoidPtr(x); }

inline PyObject* py_convert(const Complex& x)
{
  return PyComplex_FromCComplex((Py_complex&)x);
}

inline PyObject* py_convert(const Coordinate& coor)
{
  PyObject* ret = PyList_New(coor.size());
  for (long i = 0; i < (long)coor.size(); i++) {
    PyList_SetItem(ret, i, PyLong_FromLong(coor[i]));
  }
  return ret;
}

inline PyObject* py_convert(const std::string& x)
{
  return PyUnicode_FromStringAndSize(x.c_str(), x.size());
}

template <class M>
PyObject* py_convert(const std::vector<M>& vec)
{
  PyObject* ret = PyList_New(vec.size());
  for (long i = 0; i < (long)vec.size(); i++) {
    PyList_SetItem(ret, i, py_convert(vec[i]));
  }
  return ret;
}

template <class M>
PyObject* py_convert(const Vector<M>& vec)
{
  PyObject* ret = PyList_New(vec.size());
  for (long i = 0; i < (long)vec.size(); i++) {
    PyList_SetItem(ret, i, py_convert(vec[i]));
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

template <class M>
PyObject* py_convert_mview(const Vector<M>& x)
{
  return PyMemoryView_FromMemory((char*)x.data(), x.data_size(), PyBUF_WRITE);
}

}  // namespace qlat
