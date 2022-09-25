#include "lib.h"

#include <qlat-utils/cache.h>

EXPORT(random_permute, {
  using namespace qlat;
  PyObject* p_list = NULL;
  PyObject* p_rng = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_list, &p_rng)) {
    return NULL;
  }
  const RngState& rs = py_convert_type<RngState>(p_rng);
  pqassert(PyList_Check(p_list));
  const long len = PyList_Size(p_list);
  pqassert(len >= 0);
  std::vector<PyObject*> vec(len);
  for (long i = 0; i < len; ++i) {
    vec[i] = PyList_GetItem(p_list, i);
  }
  random_permute(vec, rs);
  PyObject* ret = PyList_New(len);
  for (long i = 0; i < len; ++i) {
    Py_INCREF(vec[i]);
    PyList_SetItem(ret, i, vec[i]);
  }
  return ret;
})

EXPORT(get_all_caches_info, {
  using namespace qlat;
  const std::vector<std::string> ret = get_all_caches_info();
  return py_convert(ret);
})

EXPORT(clear_all_caches, {
  using namespace qlat;
  clear_all_caches();
  Py_RETURN_NONE;
})

EXPORT(malloc_stats, {
  using namespace qlat;
  malloc_stats();
  Py_RETURN_NONE;
})
