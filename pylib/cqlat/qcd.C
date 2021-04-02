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

EXPORT(save_gauge_field, {
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
  const long ret = load_gauge_field(gf, path);
  return py_convert(ret);
});

EXPORT(gf_twist_boundary_at_boundary, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  double lmom = 0.0;
  int mu = -1;
  if (!PyArg_ParseTuple(args, "Odi", &p_gf, &lmom, &mu)) {
    return NULL;
  }
  GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  pqassert(0 <= mu and mu < 4);
  twist_boundary_at_boundary(gf, lmom, mu);
  Py_RETURN_NONE;
});

EXPORT(apply_gt_gt, {
  using namespace qlat;
  PyObject* p_gt = NULL;
  PyObject* p_gt1 = NULL;
  PyObject* p_gt0 = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_gt, &p_gt1, &p_gt0)) {
    return NULL;
  }
  // p_gt <- p_gt1 * p_gt0
  GaugeTransform& gt = py_convert_type<GaugeTransform>(p_gt);
  const GaugeTransform& gt1 = py_convert_type<GaugeTransform>(p_gt1);
  const GaugeTransform& gt0 = py_convert_type<GaugeTransform>(p_gt0);
  gt_apply_gauge_transformation(gt, gt0, gt1);
  Py_RETURN_NONE;
});

EXPORT(apply_gt_gf, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  PyObject* p_gt = NULL;
  PyObject* p_gf0 = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_gf, &p_gt, &p_gf0)) {
    return NULL;
  }
  // p_gf <- p_gt * p_gf0
  GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  const GaugeTransform& gt = py_convert_type<GaugeTransform>(p_gt);
  const GaugeField& gf0 = py_convert_type<GaugeField>(p_gf0);
  gf_apply_gauge_transformation(gf, gf0, gt);
  Py_RETURN_NONE;
});

EXPORT(apply_gt_prop, {
  using namespace qlat;
  PyObject* p_prop = NULL;
  PyObject* p_gt = NULL;
  PyObject* p_prop0 = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_prop, &p_gt, &p_prop0)) {
    return NULL;
  }
  // p_prop <- p_gt * p_prop0
  Propagator4d& prop = py_convert_type<Propagator4d>(p_prop);
  const GaugeTransform& gt = py_convert_type<GaugeTransform>(p_gt);
  const Propagator4d& prop0 = py_convert_type<Propagator4d>(p_prop0);
  prop_apply_gauge_transformation(prop, prop0, gt);
  Py_RETURN_NONE;
});

EXPORT(apply_gt_sprop, {
  using namespace qlat;
  PyObject* p_prop = NULL;
  PyObject* p_gt = NULL;
  PyObject* p_prop0 = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_prop, &p_gt, &p_prop0)) {
    return NULL;
  }
  // p_prop <- p_gt * p_prop0
  SelectedField<WilsonMatrix>& prop =
      py_convert_type<SelectedField<WilsonMatrix> >(p_prop);
  const GaugeTransform& gt = py_convert_type<GaugeTransform>(p_gt);
  const SelectedField<WilsonMatrix>& prop0 =
      py_convert_type<SelectedField<WilsonMatrix> >(p_prop0);
  PyObject* p_fsel = PyObject_GetAttrString(p_prop0, "fsel");
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_fsel);
  prop_apply_gauge_transformation(prop, prop0, gt, fsel);
  Py_RETURN_NONE;
});

EXPORT(apply_gt_psprop, {
  using namespace qlat;
  PyObject* p_prop = NULL;
  PyObject* p_gt = NULL;
  PyObject* p_prop0 = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_prop, &p_gt, &p_prop0)) {
    return NULL;
  }
  // p_prop <- p_gt * p_prop0
  SelectedPoints<WilsonMatrix>& prop =
      py_convert_type<SelectedPoints<WilsonMatrix> >(p_prop);
  const GaugeTransform& gt = py_convert_type<GaugeTransform>(p_gt);
  const SelectedPoints<WilsonMatrix>& prop0 =
      py_convert_type<SelectedPoints<WilsonMatrix> >(p_prop0);
  PyObject* p_psel = PyObject_GetAttrString(p_prop0, "psel");
  const PointSelection& psel = py_convert_type<PointSelection>(p_psel);
  prop_apply_gauge_transformation(prop, prop0, gt, psel);
  Py_RETURN_NONE;
});

EXPORT(gt_invert, {
  using namespace qlat;
  PyObject* p_gt = NULL;
  PyObject* p_gt0 = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_gt, &p_gt0)) {
    return NULL;
  }
  // p_gt <- p_gt0^{-1}
  GaugeTransform& gt = py_convert_type<GaugeTransform>(p_gt);
  const GaugeTransform& gt0 = py_convert_type<GaugeTransform>(p_gt0);
  gt_invert(gt, gt0);
  Py_RETURN_NONE;
});
