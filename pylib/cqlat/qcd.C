#include "lib.h"

EXPORT(gf_show_info, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_gf)) {
    return NULL;
  }
  const GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  gf_show_info(gf);
  Py_RETURN_NONE;
});

EXPORT(gf_avg_plaq, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_gf)) {
    return NULL;
  }
  const GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  const double ret = gf_avg_plaq(gf);
  return py_convert(ret);
});

EXPORT(gf_avg_link_trace, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_gf)) {
    return NULL;
  }
  const GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  const double ret = gf_avg_link_trace(gf);
  return py_convert(ret);
});

EXPORT(set_g_rand_color_matrix_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  PyObject* p_rng = NULL;
  double sigma = 1.0;
  int n_step = 1;
  if (!PyArg_ParseTuple(args, "OOd|i", &p_field, &p_rng, &sigma, &n_step)) {
    return NULL;
  }
  Field<ColorMatrix>& field = py_convert_type<Field<ColorMatrix> >(p_field);
  const RngState& rng = py_convert_type<RngState>(p_rng);
  set_g_rand_color_matrix_field(field, rng, sigma, n_step);
  Py_RETURN_NONE;
});

EXPORT(unitarize_color_matrix_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  Field<ColorMatrix>& field = py_convert_type<Field<ColorMatrix> >(p_field);
  unitarize(field);
  Py_RETURN_NONE;
});

EXPORT(save_gauge_fiel, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_gf, &p_path)) {
    return NULL;
  }
  const GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  std::string path;
  py_convert(path, p_path);
  const long ret = save_gauge_field(gf, path);
  return py_convert(ret);
});

EXPORT(load_gauge_field, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_gf, &p_path)) {
    return NULL;
  }
  GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  std::string path;
  py_convert(path, p_path);
  const long ret = load_gauge_field_par(gf, path);
  return py_convert(ret);
});

