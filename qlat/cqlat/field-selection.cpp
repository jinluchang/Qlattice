#include "lib.h"

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
})

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
})

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
})

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
})

EXPORT(add_psel_fsel, {
  using namespace qlat;
  PyObject* p_fsel = NULL;
  PyObject* p_psel = NULL;
  long rank_psel = 1024L * 1024L * 1024L * 1024L * 1024L;
  if (!PyArg_ParseTuple(args, "OO|l", &p_fsel, &p_psel, &rank_psel)) {
    return NULL;
  }
  FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  PointsSelection& psel = py_convert_type<PointsSelection>(p_psel);
  add_field_selection(fsel.f_rank, psel, rank_psel);
  Py_RETURN_NONE;
})

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
})

EXPORT(select_rank_range_fsel, {
  using namespace qlat;
  PyObject* p_fsel = NULL;
  PyObject* p_fsel0 = NULL;
  long rank_start = 0;
  long rank_stop = -1;
  if (!PyArg_ParseTuple(args, "OO|ll", &p_fsel, &p_fsel0, &rank_start,
                        &rank_stop)) {
    return NULL;
  }
  FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  const FieldSelection& fsel0 = py_convert_type<FieldSelection>(p_fsel0);
  fsel.f_rank = fsel0.f_rank;
  select_rank_range(fsel.f_rank, rank_start, rank_stop);
  Py_RETURN_NONE;
})

EXPORT(select_t_range_fsel, {
  using namespace qlat;
  PyObject* p_fsel = NULL;
  PyObject* p_fsel0 = NULL;
  long t_start = 0;
  long t_stop = -1;
  if (!PyArg_ParseTuple(args, "OO|ll", &p_fsel, &p_fsel0, &t_start, &t_stop)) {
    return NULL;
  }
  FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  const FieldSelection& fsel0 = py_convert_type<FieldSelection>(p_fsel0);
  fsel.f_rank = fsel0.f_rank;
  select_t_range(fsel.f_rank, t_start, t_stop);
  Py_RETURN_NONE;
})

EXPORT(set_psel_fsel, {
  using namespace qlat;
  PyObject* p_psel = NULL;
  PyObject* p_fsel = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_psel, &p_fsel)) {
    return NULL;
  }
  PointsSelection& psel = py_convert_type<PointsSelection>(p_psel);
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  psel = psel_from_fsel(fsel);
  Py_RETURN_NONE;
})

EXPORT(set_psel_fsel_local, {
  using namespace qlat;
  PyObject* p_psel = NULL;
  PyObject* p_fsel = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_psel, &p_fsel)) {
    return NULL;
  }
  PointsSelection& psel = py_convert_type<PointsSelection>(p_psel);
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  psel = psel_from_fsel_local(fsel);
  Py_RETURN_NONE;
})

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
})

EXPORT(get_total_site_fsel, {
  using namespace qlat;
  PyObject* p_fsel = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_fsel)) {
    return NULL;
  }
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  return py_convert(fsel.f_rank.geo().total_site());
})

EXPORT(get_n_elems_fsel, {
  using namespace qlat;
  PyObject* p_fsel = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_fsel)) {
    return NULL;
  }
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  const long n_elems = fsel.n_elems;
  return py_convert(n_elems);
})

EXPORT(get_n_per_tslice_fsel, {
  using namespace qlat;
  PyObject* p_fsel = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_fsel)) {
    return NULL;
  }
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  const long n_per_tslice= fsel.n_per_tslice;
  return py_convert(n_per_tslice);
})

EXPORT(get_prob_fsel, {
  using namespace qlat;
  PyObject* p_fsel = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_fsel)) {
    return NULL;
  }
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  const double prob = fsel.prob;
  return py_convert(prob);
})

EXPORT(get_coordinate_from_idx_fsel, {
  // return global coordinate
  using namespace qlat;
  PyObject* p_fsel = NULL;
  long idx = -1;
  if (!PyArg_ParseTuple(args, "Ol", &p_fsel, &idx)) {
    return NULL;
  }
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  const Geometry& geo = fsel.f_rank.geo();
  const long index = fsel.indices[idx];
  const Coordinate xl = geo.coordinate_from_index(index);
  const Coordinate xg = geo.coordinate_g_from_l(xl);
  return py_convert(xg);
})

EXPORT(get_idx_from_coordinate_fsel, {
  // need global coordinate
  using namespace qlat;
  PyObject* p_fsel = NULL;
  PyObject* p_xg = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_fsel, &p_xg)) {
    return NULL;
  }
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  const Geometry& geo = fsel.f_rank.geo();
  const Coordinate xg = py_convert_data<Coordinate>(p_xg);
  const Coordinate xl = geo.coordinate_l_from_g(xg);
  const long idx = fsel.f_local_idx.get_elem(xl);
  return py_convert(idx);
})
