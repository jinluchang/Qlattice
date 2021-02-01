#pragma once

// From https://github.com/lehner/gpt/blob/master/lib/cgpt/lib/exception.h
// Original author Christoph Lehner

#include <Python.h>
#include <qlat/qlat.h>

#define STRX(x) #x

#define STR(x) STRX(x)

#define pqassert(x)                                              \
  {                                                              \
    if (!(x))                                                    \
      throw std::string("Assert " #x " failed in file " __FILE__ \
                        ":" STR(__LINE__));                      \
  };

#define pqerr(...)                                                      \
  {                                                                   \
    const std::string msg =                                           \
        qlat::ssprintf(msg, __VA_ARGS__) +                            \
        qlat::ssprintf(" in from '%s' line %d ", __FILE__, __LINE__); \
    throw std::string(msg);                                           \
  };

#define EXPORT(name, ...)                               \
  PyObject* cqlat_##name(PyObject* self, PyObject* args) \
  {                                                     \
    try {                                               \
      __VA_ARGS__;                                      \
      return NULL;                                      \
    } catch (std::string err) {                         \
      fprintf(stderr, "ERR: %s\n", err.c_str());        \
      PyErr_SetString(PyExc_RuntimeError, err.c_str()); \
      return NULL;                                      \
    }                                                   \
  }
