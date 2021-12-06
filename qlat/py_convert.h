#pragma once

#include <Python.h>
#include <qlat/qlat.h>

namespace qlat
{  //

template <class M>
std::string get_ctype_name()
{
  return "M";
}

template <>
inline std::string get_ctype_name<ColorMatrix>()
{
  return "ColorMatrix";
}

template <>
inline std::string get_ctype_name<WilsonMatrix>()
{
  return "WilsonMatrix";
}

template <>
inline std::string get_ctype_name<NonRelWilsonMatrix>()
{
  return "NonRelWilsonMatrix";
}

template <>
inline std::string get_ctype_name<SpinMatrix>()
{
  return "SpinMatrix";
}

template <>
inline std::string get_ctype_name<WilsonVector>()
{
  return "WilsonVector";
}

template <>
inline std::string get_ctype_name<Complex>()
{
  return "Complex";
}

template <>
inline std::string get_ctype_name<double>()
{
  return "double";
}

template <>
inline std::string get_ctype_name<float>()
{
  return "float";
}

template <>
inline std::string get_ctype_name<int64_t>()
{
  return "int64_t";
}

template <>
inline std::string get_ctype_name<char>()
{
  return "char";
}

template <>
inline std::string get_ctype_name<int8_t>()
{
  return "int8_t";
}

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
// py_convert_data<PyField>(in)
{
  T x;
  py_convert(x, in);
  return x;
}

template <class T>
T py_convert_data(PyObject* in, const std::string& attr)
// py_convert_data<std::string>(in, "ctype")
{
  return py_convert_data<T>(PyObject_GetAttrString(in, attr.c_str()));
}

inline std::string py_get_ctype(PyObject* in)
{
  return py_convert_data<std::string>(in, "ctype");
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
// py_convert_type<CommPlan>(in);
// specifications:
// py_convert_type<CommMarks>(in);
{
  PyObject* p_cdata = PyObject_GetAttrString(in, "cdata");
  T* out;
  py_convert((long&)out, p_cdata);
  return *out;
}

template <class M>
Field<M>& py_convert_type_field(PyObject* in)
// py_convert_ftype<Field<M> >(in);
{
  PyField pf = py_convert_field(in);
  pqassert(pf.ctype == get_ctype_name<M>());
  Field<M>& f = *(Field<M>*)pf.cdata;
  return f;
}

template <class M, int multiplicity>
FieldM<M, multiplicity>& py_convert_type_field(PyObject* in)
// py_convert_ftype<FieldM<M, multiplicity> >(in);
{
  PyField pf = py_convert_field(in);
  pqassert(pf.ctype == get_ctype_name<M>());
  FieldM<M, multiplicity>& f = *(FieldM<M, multiplicity>*)pf.cdata;
  pqassert(multiplicity == f.geo().multiplicity);
  return f;
}

template <>
inline CommMarks& py_convert_type<CommMarks>(PyObject* in)
{
  Field<int8_t>& f = py_convert_type_field<int8_t>(in);
  CommMarks& marks = static_cast<CommMarks&>(f);
  return marks;
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
