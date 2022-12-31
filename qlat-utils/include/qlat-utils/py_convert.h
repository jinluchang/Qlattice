#pragma once

#include <qlat-utils/core.h>
#include <qlat-utils/qutils.h>
#include <qlat-utils/qutils-vec.h>

#include <Python.h>

namespace qlat
{  //

template <class M>
bool check_ctype_name(const std::string& ctype)
{
  (void)ctype;
  return false;
}

template <>
inline bool check_ctype_name<Complex>(const std::string& ctype)
{
  return "Complex" == ctype;
}

template <>
inline bool check_ctype_name<ComplexF>(const std::string& ctype)
{
  return "ComplexF" == ctype;
}

template <>
inline bool check_ctype_name<double>(const std::string& ctype)
{
  return "Double" == ctype;
}

template <>
inline bool check_ctype_name<float>(const std::string& ctype)
{
  return "Float" == ctype;
}

template <>
inline bool check_ctype_name<int64_t>(const std::string& ctype)
{
  return "Int64t" == ctype or "Long" == ctype;
}

template <>
inline bool check_ctype_name<char>(const std::string& ctype)
{
  return "Char" == ctype;
}

template <>
inline bool check_ctype_name<int8_t>(const std::string& ctype)
{
  return "Int8t" == ctype;
}

inline void py_convert(PyObject*& out, PyObject* in) { out = in; }

inline void py_convert(int& out, PyObject* in)
{
  qassert(PyLong_Check(in));
  out = PyLong_AsLong(in);
}

inline void py_convert(long& out, PyObject* in)
{
  qassert(PyLong_Check(in));
  out = PyLong_AsLong(in);
}

inline void py_convert(double& out, PyObject* in)
{
  if (PyFloat_Check(in)) {
    out = PyFloat_AsDouble(in);
  } else if (PyLong_Check(in)) {
    out = PyLong_AsLong(in);
  } else {
    qassert(false);
  }
}

inline void py_convert(Complex& out, PyObject* in)
{
  if (PyLong_Check(in)) {
    out = PyLong_AsLong(in);
  } else if (PyFloat_Check(in)) {
    out = PyFloat_AsDouble(in);
  } else {
    qassert(PyComplex_Check(in));
    Py_complex& py_out = (Py_complex&)out;
    py_out = PyComplex_AsCComplex(in);
  }
}

inline void py_convert(bool& out, PyObject* in)
{
  qassert(PyBool_Check(in));
  out = in == Py_True;
}

inline void py_convert(std::string& s, PyObject* in)
{
  if (PyType_Check(in)) {
    s = ((PyTypeObject*)in)->tp_name;
  } else if (PyBytes_Check(in)) {
    const long size = PyBytes_Size(in);
    s = std::string(PyBytes_AsString(in), size);
  } else if (PyUnicode_Check(in)) {
    PyObject* temp = PyUnicode_AsEncodedString(in, "UTF-8", "strict");
    qassert(temp);
    s = PyBytes_AS_STRING(temp);
    Py_DECREF(temp);
  } else {
    qassert(false);
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
    qassert(false);
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
    qassert(false);
  }
}

template <class M>
void py_convert(Vector<M> out, PyObject* in)
{
  if (PyList_Check(in)) {
    qassert(out.size() == PyList_Size(in));
    for (long i = 0; i < out.size(); i++) {
      py_convert(out[i], PyList_GetItem(in, i));
    }
  } else if (PyTuple_Check(in)) {
    qassert(out.size() == PyTuple_Size(in));
    for (long i = 0; i < out.size(); i++) {
      py_convert(out[i], PyTuple_GetItem(in, i));
    }
  } else {
    qassert(false);
  }
}

template <class M, unsigned long N>
void py_convert(array<M, N>& out, PyObject* in)
{
  if (PyList_Check(in)) {
    qassert(out.size() == PyList_Size(in));
    for (long i = 0; i < out.size(); i++) {
      py_convert(out[i], PyList_GetItem(in, i));
    }
  } else if (PyTuple_Check(in)) {
    qassert(out.size() == PyTuple_Size(in));
    for (long i = 0; i < out.size(); i++) {
      py_convert(out[i], PyTuple_GetItem(in, i));
    }
  } else {
    qassert(false);
  }
}

template <class T>
T py_convert_data(PyObject* in)
// interface
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
// interface
// py_convert_data<std::string>(in, "ctype")
// py_convert_data<long>(in, "cdata")
{
  PyObject* p_obj = PyObject_GetAttrString(in, attr.c_str());
  T x = py_convert_data<T>(p_obj);
  Py_DECREF(p_obj);
  return x;
}

template <class T>
T py_convert_data(PyObject* in, const std::string& attr, const std::string& attr1)
// interface
{
  PyObject* p_obj = PyObject_GetAttrString(in, attr.c_str());
  T x = py_convert_data<T>(p_obj, attr1);
  Py_DECREF(p_obj);
  return x;
}

inline std::string py_get_ctype(PyObject* in)
// interface
// py_convert_data<std::string>(in, "ctype", "name")
{
  return py_convert_data<std::string>(in, "ctype", "name");
}

template <class T>
T& py_convert_type(PyObject* in)
// interface
// use cdata property of PyObject* in as pointer
// examples:
// py_convert_type<Geometry>(in);
// py_convert_type<LatData>(in);
// py_convert_type<RngState>(in);
// py_convert_type<PointSelection>(in);
// py_convert_type<CommPlan>(in);
// py_convert_type<FieldSelection>(in);
// specifications:
// py_convert_type<Propagator4d>(in);
// py_convert_type<GaugeField>(in);
// py_convert_type<CommMarks>(in);
{
  T* out = (T*)py_convert_data<long>(in, "cdata");
  return *out;
}

template <class T>
T& py_convert_type(PyObject* in, const std::string& attr)
// interface
// py_convert_type<PointSelection>(in, "psel")
// py_convert_type<FieldSelection>(in, "fsel")
// py_convert_type<Geometry>(in, "geo")
{
  PyObject* p_obj = PyObject_GetAttrString(in, attr.c_str());
  T& x = py_convert_type<T>(p_obj);
  Py_DECREF(p_obj);
  return x;
}

template <class T>
T& py_convert_type(PyObject* in, const std::string& attr, const std::string& attr1)
// interface
// py_convert_type<Geometry>(in, "psel", "geo")
{
  PyObject* p_obj = PyObject_GetAttrString(in, attr.c_str());
  T& x = py_convert_type<T>(p_obj, attr1);
  Py_DECREF(p_obj);
  return x;
}

// -------------------------------------------------------------

inline PyObject* py_convert(PyObject* x) { return x; }

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

inline PyObject* py_convert(const Complex& x)
{
  return PyComplex_FromCComplex((Py_complex&)x);
}

inline PyObject* py_convert(void* x) { return PyLong_FromVoidPtr(x); }

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

template <>
inline PyObject* py_convert(const Vector<char>& vec)
{
  return PyBytes_FromStringAndSize(vec.data(), vec.size());
}

template <class M>
PyObject* py_convert_mview(const Vector<M>& x)
{
  return PyMemoryView_FromMemory((char*)x.data(), x.data_size(), PyBUF_WRITE);
}

}  // namespace qlat
