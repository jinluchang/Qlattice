#pragma once

#include <Python.h>
#include <qlat/py_exceptions.h>
#include <qlat/qlat.h>

namespace qlat
{  //

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
  pqassert(PyFloat_Check(in));
  out = PyFloat_AsDouble(in);
}

inline void py_convert(Complex& out, PyObject* in)
{
  if (PyFloat_Check(in)) {
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

inline PyField py_convert_field(PyObject* in)
{
  PyField out;
  py_convert(out, in);
  return out;
}

template <class T>
T& py_convert_type(PyObject* in)
{
  PyObject* p_cdata = PyObject_GetAttrString(in, "cdata");
  T* out;
  py_convert((long&)out, p_cdata);
  return *out;
}

inline PyObject* py_convert(const Coordinate& coor)
{
  PyObject* ret = PyList_New(coor.size());
  for (long i = 0; i < (long)coor.size(); i++) {
    PyList_SetItem(ret, i, PyLong_FromLong(coor[i]));
  }
  return ret;
}

inline PyObject* py_convert(const int& x) { return PyLong_FromLong((long)x); }

inline PyObject* py_convert(const long& x) { return PyLong_FromLong(x); }

inline PyObject* py_convert(const long long& x) { return PyLong_FromLongLong(x); }

inline PyObject* py_convert(const unsigned int& x) { return PyLong_FromUnsignedLong((unsigned long)x); }

inline PyObject* py_convert(const unsigned long& x) { return PyLong_FromUnsignedLong(x); }

inline PyObject* py_convert(const unsigned long long& x) { return PyLong_FromUnsignedLongLong(x); }

inline PyObject* py_convert(const double& x) { return PyFloat_FromDouble(x); }

inline PyObject* py_convert(void* x) { return PyLong_FromVoidPtr(x); }

inline PyObject* py_convert(const Complex& x)
{
  return PyComplex_FromCComplex((Py_complex&)x);
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

template <class M>
PyObject* py_convert_mview(const Vector<M>& x)
{
  return PyMemoryView_FromMemory((char*)x.data(), x.data_size(), PyBUF_WRITE);
}

}  // namespace qlat
