#pragma once

#include <Python.h>
#include <qlat/qlat.h>
#include "exceptions.h"

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

template <typename t>
void py_convert(std::vector<t>& out, PyObject* in)
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

inline void py_convert(qlat::Coordinate& out, PyObject* in)
{
  if (PyList_Check(in)) {
    pqassert(qlat::DIMN == PyList_Size(in));
    for (size_t i = 0; i < out.size(); i++) {
      py_convert(out[i], PyList_GetItem(in, i));
    }
  } else if (PyTuple_Check(in)) {
    pqassert(qlat::DIMN == PyTuple_Size(in));
    for (size_t i = 0; i < out.size(); i++) {
      py_convert(out[i], PyTuple_GetItem(in, i));
    }
  } else {
    pqassert(false);
  }
}

inline PyObject* py_convert(const qlat::Coordinate& coor)
{
  PyObject* ret = PyList_New(coor.size());
  for (long i = 0; i < (long)coor.size(); i++) {
    PyList_SetItem(ret, i, PyLong_FromLong(coor[i]));
  }
  return ret;
}

}  // namespace qlat
