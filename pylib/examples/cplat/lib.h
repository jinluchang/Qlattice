#pragma once

// From https://github.com/lehner/gpt/blob/master/lib/cgpt/lib/exception.h
// Original author Christoph Lehner

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <qlat/py_convert.h>

#define EXPORT(name, ...)                                  \
  PyObject* cplat##_##name(PyObject* self, PyObject* args) \
  {                                                        \
    try {                                                  \
      __VA_ARGS__;                                         \
      return NULL;                                         \
    } catch (std::string err) {                            \
      fprintf(stderr, "ERR: %s\n", err.c_str());           \
      PyErr_SetString(PyExc_RuntimeError, err.c_str());    \
      return NULL;                                         \
    }                                                      \
  }
