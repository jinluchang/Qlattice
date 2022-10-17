
// See https://github.com/lehner/gpt/blob/master/lib/cgpt/lib/lib.cc
// Original author Christoph Lehner

#include "lib.h"

// declare
#define EXPORT_FUNCTION(name) \
  extern PyObject* PKG_PREFIX(name)(PyObject* self, PyObject * args);
extern "C" {
#include "exports.h"
}
#undef EXPORT_FUNCTION

// add to module functions
#define EXPORT_FUNCTION(name) {#name, PKG_PREFIX(name), METH_VARARGS, #name},
static PyMethodDef module_functions[] = {
#include "exports.h"
    {NULL, NULL, 0, NULL}};
#undef EXPORT_FUNCTION

// on exit
void free_module(void* self) { (void)self; }

// module definition
static struct PyModuleDef module_def = {
    PyModuleDef_HEAD_INIT,
    PSTR(PY_PKG_NAME),                          /* m_name */
    "The C++ interface for " PSTR(PY_PKG_NAME), /* m_doc */
    -1,                                         /* m_size */
    module_functions,                           /* m_methods */
    NULL,                                       /* m_reload */
    NULL,                                       /* m_traverse */
    NULL,                                       /* m_clear */
    free_module,                                /* m_free */
};

// export module creation
#define PY_INIT_PY_PKG_NAME(name) PY_INIT_PY_PKG_NAME_X(name)
#define PY_INIT_PY_PKG_NAME_X(name) PyInit_##name
PyMODINIT_FUNC PY_INIT_PY_PKG_NAME(PY_PKG_NAME)(void)
{
  return PyModule_Create(&module_def);
}
#undef PY_INIT_PY_PKG_NAME_X
#undef PY_INIT_PY_PKG_NAME
