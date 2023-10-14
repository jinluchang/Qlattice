#pragma once

// From https://github.com/lehner/gpt/blob/master/lib/cgpt/lib/exception.h
// Original author Christoph Lehner
// With modifications from Luchang Jin

#define PY_PKG_NAME cqlat

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <qlat-utils/show.h>
#include <qlat/py_convert.h>

#define PKG_PREFIX(name) PKG_PREFIX_X(PY_PKG_NAME, name)
#define PKG_PREFIX_X(pname, name) PKG_PREFIX_XX(pname, name)
#define PKG_PREFIX_XX(pname, name) pname##_##name

#define EXPORT(name, ...)                                      \
                                                               \
  extern "C" {                                                 \
  PyObject* PKG_PREFIX(name)(PyObject* self, PyObject* args);  \
  }                                                            \
                                                               \
  PyObject* PKG_PREFIX(name)(PyObject* self, PyObject* args)   \
  {                                                            \
    (void)self;                                                \
    (void)args;                                                \
    try {                                                      \
      {__VA_ARGS__};                                           \
      return NULL;                                             \
    } catch (std::string err) {                                \
      fprintf(stderr, "ERR: %s\n", err.c_str());               \
      PyErr_SetString(PyExc_RuntimeError, err.c_str());        \
      return NULL;                                             \
    }                                                          \
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

template <class T>
PyObject* set_obj(PyObject* args)
{
  PyObject* p_obj_new = NULL;
  PyObject* p_obj = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_obj_new, &p_obj)) {
    return NULL;
  }
  T& obj_new = py_convert_type<T>(p_obj_new);
  const T& obj = py_convert_type<T>(p_obj);
  obj_new = obj;
  Py_RETURN_NONE;
}

}  // namespace qlat
