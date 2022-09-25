
// See https://github.com/lehner/gpt/blob/master/lib/cgpt/lib/lib.cc
// Original author Christoph Lehner

#include "lib.h"

// declare
#define EXPORT_FUNCTION(name) \
  PyObject* cqlat_utils_##name(PyObject* self, PyObject* args);
extern "C" {
#include "exports.h"
}
#undef EXPORT_FUNCTION

// add to module functions
#define EXPORT_FUNCTION(name) {#name, cqlat_utils_##name, METH_VARARGS, #name},
static PyMethodDef module_functions[] = {
#include "exports.h"
    {NULL, NULL, 0, NULL}};
#undef EXPORT_FUNCTION

// on exit
void free_module(void* self) { (void)self; }

// module definition
static struct PyModuleDef module_def = {
    PyModuleDef_HEAD_INIT,
    "cqlat_utils",                      /* m_name */
    "The C++ interface for qlat_utils", /* m_doc */
    -1,                                 /* m_size */
    module_functions,                   /* m_methods */
    NULL,                               /* m_reload */
    NULL,                               /* m_traverse */
    NULL,                               /* m_clear */
    free_module,                        /* m_free */
};

// export module creation
PyMODINIT_FUNC PyInit_cqlat_utils(void) { return PyModule_Create(&module_def); }
