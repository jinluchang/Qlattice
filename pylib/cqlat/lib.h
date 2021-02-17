#pragma once

// From https://github.com/lehner/gpt/blob/master/lib/cgpt/lib/exception.h
// Original author Christoph Lehner
// With modifications from Luchang Jin

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <qlat/py_convert.h>
#include "dispatch.h"

#define EXPORT(name, ...)                                \
  PyObject* cqlat_##name(PyObject* self, PyObject* args) \
  {                                                      \
    try {                                                \
      __VA_ARGS__;                                       \
      return NULL;                                       \
    } catch (std::string err) {                          \
      fprintf(stderr, "ERR: %s\n", err.c_str());         \
      PyErr_SetString(PyExc_RuntimeError, err.c_str());  \
      return NULL;                                       \
    }                                                    \
  }

namespace qlat
{  //

template <class T>
PyObject* free_obj(PyObject* args)
{
  PyObject* p_obj = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_obj)) {
    return NULL;
  }
  T& obj = py_convert_type<T>(p_obj);
  delete &obj;
  Py_RETURN_NONE;
}

}  // namespace qlat
