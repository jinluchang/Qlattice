#pragma once

// From https://github.com/lehner/gpt/blob/master/lib/cgpt/lib/exception.h
// Original author Christoph Lehner
// With modifications from Luchang Jin

#define PY_PKG_NAME cqlat_grid_io

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <qlat-utils/show.h>
#include <qlat/py_convert.h>

#define EXPORT(name, ...) EXPORT_X(PY_PKG_NAME, name, {__VA_ARGS__})
#define EXPORT_X(pname, name, ...) EXPORT_XX(pname, name, {__VA_ARGS__})
#define EXPORT_XX(pname, name, ...)                         \
                                                            \
  extern "C" {                                              \
  PyObject* pname##_##name(PyObject* self, PyObject* args); \
  }                                                         \
                                                            \
  PyObject* pname##_##name(PyObject* self, PyObject* args)  \
  {                                                         \
    (void)self;                                             \
    (void)args;                                             \
    try {                                                   \
      {__VA_ARGS__};                                        \
      return NULL;                                          \
    } catch (std::string err) {                             \
      fprintf(stderr, "ERR: %s\n", err.c_str());            \
      PyErr_SetString(PyExc_RuntimeError, err.c_str());     \
      return NULL;                                          \
    }                                                       \
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