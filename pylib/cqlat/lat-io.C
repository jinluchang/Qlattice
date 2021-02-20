#include "lib.h"

EXPORT(mk_lat_data, {
  using namespace qlat;
  LatData* pld = new LatData();
  return py_convert((void*)pld);
});

EXPORT(free_lat_data, {
  using namespace qlat;
  return free_obj<LatData>(args);
});

EXPORT(set_lat_data, {
  using namespace qlat;
  return set_obj<LatData>(args);
});

EXPORT(set_zero_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_ld)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  set_zero(ld);
  Py_RETURN_NONE;
});

EXPORT(show_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_ld)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  return py_convert(show(ld));
});

EXPORT(load_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_ld, &p_path)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  std::string path;
  py_convert(path, p_path);
  ld.load(path);
  Py_RETURN_NONE;
});

EXPORT(save_lat_data, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_ld, &p_path)) {
    return NULL;
  }
  const LatData& ld = py_convert_type<LatData>(p_ld);
  std::string path;
  py_convert(path, p_path);
  ld.save(path);
  Py_RETURN_NONE;
});

