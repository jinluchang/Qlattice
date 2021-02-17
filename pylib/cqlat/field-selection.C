#include "lib.h"

EXPORT(mk_psel, {
  using namespace qlat;
  PyObject* p_coordinate_list = NULL;
  if (!PyArg_ParseTuple(args, "|O", &p_coordinate_list)) {
    return NULL;
  }
  PointSelection* ppsel = new PointSelection();
  if (NULL != p_coordinate_list) {
    PointSelection& psel = *ppsel;
    py_convert(psel, p_coordinate_list);
  }
  return py_convert((void*)ppsel);
});

EXPORT(free_psel, {
  using namespace qlat;
  return free_obj<PointSelection>(args);
});

EXPORT(load_psel, {
  using namespace qlat;
  PyObject* p_psel = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_psel, &p_path)) {
    return NULL;
  }
  PointSelection& psel = py_convert_type<PointSelection>(p_psel);
  std::string path;
  py_convert(path, p_path);
  psel = load_point_selection_info(path);
  Py_RETURN_NONE;
});

EXPORT(save_psel, {
  using namespace qlat;
  PyObject* p_psel = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_psel, &p_path)) {
    return NULL;
  }
  const PointSelection& psel = py_convert_type<PointSelection>(p_psel);
  std::string path;
  py_convert(path, p_path);
  save_point_selection_info(psel, path);
  Py_RETURN_NONE;
});

EXPORT(mk_fsel, {
  using namespace qlat;
  PyObject* p_total_site = NULL;
  long n_per_tslice = 0;
  PyObject* p_rng = NULL;
  PyObject* p_psel = NULL;
  if (!PyArg_ParseTuple(args, "|OlOO", &p_total_site, &n_per_tslice, &p_rng,
                        &p_psel)) {
    return NULL;
  }
  FieldSelection* pfsel = new FieldSelection();
  if (NULL != p_total_site) {
    FieldSelection& fsel = *pfsel;
    Coordinate total_site;
    py_convert(total_site, p_total_site);
    pqassert(NULL != p_rng);
    RngState& rs = py_convert_type<RngState>(p_rng);
    mk_field_selection(fsel.f_rank, total_site, n_per_tslice, rs);
    if (NULL != p_psel) {
      PointSelection& psel = py_convert_type<PointSelection>(p_psel);
      add_field_selection(fsel.f_rank, psel);
    }
  }
  return py_convert((void*)pfsel);
});

EXPORT(free_fsel, {
  using namespace qlat;
  return free_obj<FieldSelection>(args);
});

EXPORT(load_fsel, {
  using namespace qlat;
  PyObject* p_fsel = NULL;
  PyObject* p_path = NULL;
  long n_per_tslice = 0;
  if (!PyArg_ParseTuple(args, "OOl", &p_fsel, &p_path, &n_per_tslice)) {
    return NULL;
  }
  FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  std::string path;
  py_convert(path, p_path);
  read_field_selection(fsel, path, n_per_tslice);
  Py_RETURN_NONE;
});

EXPORT(save_fsel, {
  using namespace qlat;
  PyObject* p_fsel = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_fsel, &p_path)) {
    return NULL;
  }
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  std::string path;
  py_convert(path, p_path);
  write_field_selection(fsel, path);
  Py_RETURN_NONE;
});

EXPORT(add_psel_fsel, {
  using namespace qlat;
  PyObject* p_fsel = NULL;
  PyObject* p_psel = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_fsel, &p_psel)) {
    return NULL;
  }
  FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  PointSelection& psel = py_convert_type<PointSelection>(p_psel);
  add_field_selection(fsel.f_rank, psel);
  Py_RETURN_NONE;
});

EXPORT(update_fsel, {
  using namespace qlat;
  PyObject* p_fsel = NULL;
  long n_per_tslice = -1;
  if (!PyArg_ParseTuple(args, "O|l", &p_fsel, &n_per_tslice)) {
    return NULL;
  }
  FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  if (n_per_tslice < 0) {
    // only update various indices
    update_field_selection(fsel);
  } else {
    // only update parameters
    update_field_selection(fsel, n_per_tslice);
  }
  Py_RETURN_NONE;
});
