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

EXPORT(set_psel, {
  using namespace qlat;
  return set_obj<PointSelection>(args);
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

EXPORT(mk_list_psel, {
  using namespace qlat;
  PyObject* p_psel = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_psel)) {
    return NULL;
  }
  const PointSelection& psel = py_convert_type<PointSelection>(p_psel);
  return py_convert(psel);
});

EXPORT(set_list_psel, {
  using namespace qlat;
  PyObject* p_psel = NULL;
  PyObject* p_coordinate_list = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_psel, &p_coordinate_list)) {
    return NULL;
  }
  PointSelection& psel = py_convert_type<PointSelection>(p_psel);
  py_convert(psel, p_coordinate_list);
  Py_RETURN_NONE;
});

EXPORT(set_rand_psel, {
  using namespace qlat;
  PyObject* p_psel = NULL;
  PyObject* p_rng = NULL;
  PyObject* p_total_site = NULL;
  long n_points = 0;
  if (!PyArg_ParseTuple(args, "OOOl", &p_psel, &p_rng, &p_total_site, &n_points)) {
    return NULL;
  }
  PointSelection& psel = py_convert_type<PointSelection>(p_psel);
  const RngState& rs = py_convert_type<RngState>(p_rng);
  Coordinate total_site;
  py_convert(total_site, p_total_site);
  psel = mk_random_point_selection(total_site, n_points, rs);
  Py_RETURN_NONE;
});

EXPORT(set_tslice_psel, {
  using namespace qlat;
  PyObject* p_psel = NULL;
  PyObject* p_total_site = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_psel, &p_total_site)) {
    return NULL;
  }
  PointSelection& psel = py_convert_type<PointSelection>(p_psel);
  Coordinate total_site;
  py_convert(total_site, p_total_site);
  psel = mk_tslice_point_selection(total_site);
  Py_RETURN_NONE;
});

EXPORT(mk_fsel, {
  using namespace qlat;
  FieldSelection* pfsel = new FieldSelection();
  return py_convert((void*)pfsel);
});

EXPORT(free_fsel, {
  using namespace qlat;
  return free_obj<FieldSelection>(args);
});

EXPORT(set_fsel, {
  using namespace qlat;
  return set_obj<FieldSelection>(args);
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
  const long ret = read_field_selection(fsel, path, n_per_tslice);
  return py_convert(ret);
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
  const long ret = write_field_selection(fsel, path);
  return py_convert(ret);
});

EXPORT(set_uniform_fsel, {
  using namespace qlat;
  PyObject* p_fsel = NULL;
  PyObject* p_total_site = NULL;
  int64_t val = 0;
  if (!PyArg_ParseTuple(args, "OO|l", &p_fsel, &p_total_site, &val)) {
    return NULL;
  }
  FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  Coordinate total_site;
  py_convert(total_site, p_total_site);
  mk_field_selection(fsel.f_rank, total_site, val);
  Py_RETURN_NONE;
});

EXPORT(set_rand_fsel, {
  using namespace qlat;
  PyObject* p_fsel = NULL;
  PyObject* p_rng = NULL;
  PyObject* p_total_site = NULL;
  long n_per_tslice = 0;
  if (!PyArg_ParseTuple(args, "OOOl", &p_fsel, &p_rng, &p_total_site,
                        &n_per_tslice)) {
    return NULL;
  }
  FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  const RngState& rs = py_convert_type<RngState>(p_rng);
  Coordinate total_site;
  py_convert(total_site, p_total_site);
  mk_field_selection(fsel.f_rank, total_site, n_per_tslice, rs);
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

EXPORT(set_geo_fsel, {
  using namespace qlat;
  PyObject* p_geo = NULL;
  PyObject* p_fsel = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_geo, &p_fsel)) {
    return NULL;
  }
  Geometry& geo = py_convert_type<Geometry>(p_geo);
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  geo = fsel.f_rank.geo();
  Py_RETURN_NONE;
});

EXPORT(get_total_site_fsel, {
  using namespace qlat;
  PyObject* p_fsel = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_fsel)) {
    return NULL;
  }
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  return py_convert(fsel.f_rank.geo().total_site());
});

EXPORT(get_n_elems_fsel, {
  using namespace qlat;
  PyObject* p_fsel = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_fsel)) {
    return NULL;
  }
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  const long n_elems = fsel.n_elems;
  return py_convert(n_elems);
});

EXPORT(get_n_per_tslice_fsel, {
  using namespace qlat;
  PyObject* p_fsel = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_fsel)) {
    return NULL;
  }
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  const long n_per_tslice= fsel.n_per_tslice;
  return py_convert(n_per_tslice);
});

EXPORT(get_prob_fsel, {
  using namespace qlat;
  PyObject* p_fsel = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_fsel)) {
    return NULL;
  }
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  const double prob = fsel.prob;
  return py_convert(prob);
});
