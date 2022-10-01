
// See https://github.com/lehner/gpt/blob/master/lib/cgpt/lib/lib.cc
// Original author Christoph Lehner

#include "lib.h"

#define EXPORT_FUNCTION(name) EXPORT_FUNCTION_X(PY_PKG_NAME, name)
#define EXPORT_FUNCTION_X(pname, name) EXPORT_FUNCTION_XX(pname, name)

// declare
#define EXPORT_FUNCTION_XX(pname, name) \
  extern PyObject* pname##_##name(PyObject* self, PyObject* args);
#include "exports.h"
#undef EXPORT_FUNCTION_XX

// add to module functions
#define EXPORT_FUNCTION_XX(pname, name) {#name, pname##_##name, METH_VARARGS, #name},
static PyMethodDef module_functions[] = {
#include "exports.h"
    {NULL, NULL, 0, NULL}};
#undef EXPORT_FUNCTION_XX

#undef EXPORT_FUNCTION_X
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
#define PY_INIT_PY_PKG_NAME_X(name) PyInit_##name
#define PY_INIT_PY_PKG_NAME(name) PY_INIT_PY_PKG_NAME_X(name)
PyMODINIT_FUNC PY_INIT_PY_PKG_NAME(PY_PKG_NAME)(void)
{
  return PyModule_Create(&module_def);
}
#undef PY_INIT_PY_PKG_NAME_X
#undef PY_INIT_PY_PKG_NAME
