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
})

EXPORT(gf_avg_plaq, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_gf)) {
    return NULL;
  }
  const GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  const double ret = gf_avg_plaq(gf);
  return py_convert(ret);
})

EXPORT(gf_avg_link_trace, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_gf)) {
    return NULL;
  }
  const GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  const double ret = gf_avg_link_trace(gf);
  return py_convert(ret);
})

EXPORT(gf_avg_wilson_loop_normalized_tr, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  int l = 0;
  int t = 0;
  if (!PyArg_ParseTuple(args, "Oii", &p_gf, &l, &t)) {
    return NULL;
  }
  const GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  const double ret = matrix_trace(gf_avg_wilson_loop(gf, l, t)).real() / 3.0;
  return py_convert(ret);
})

EXPORT(gf_wilson_line_no_comm, {
  using namespace qlat;
  PyObject* p_wilson_line_field = NULL;
  int wilson_line_field_m = 0;
  PyObject* p_gf_ext = NULL;
  PyObject* p_path = NULL;
  PyObject* p_path_n = NULL;
  if (!PyArg_ParseTuple(args, "OiOO|O", &p_wilson_line_field,
                        &wilson_line_field_m, &p_gf_ext, &p_path, &p_path_n)) {
    return NULL;
  }
  Field<ColorMatrix>& wilson_line_field =
      py_convert_type_field<ColorMatrix>(p_wilson_line_field);
  const GaugeField& gf_ext = py_convert_type<GaugeField>(p_gf_ext);
  const std::vector<int> path = py_convert_data<std::vector<int> >(p_path);
  if (p_path_n != NULL) {
    const std::vector<int> path_n =
        py_convert_data<std::vector<int> >(p_path_n);
    gf_wilson_line_no_comm(wilson_line_field, wilson_line_field_m, gf_ext, path,
                           path_n);
  } else {
    gf_wilson_line_no_comm(wilson_line_field, wilson_line_field_m, gf_ext,
                           path);
  }
  Py_RETURN_NONE;
})

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
})

EXPORT(unitarize_color_matrix_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_field)) {
    return NULL;
  }
  Field<ColorMatrix>& field = py_convert_type<Field<ColorMatrix> >(p_field);
  unitarize(field);
  Py_RETURN_NONE;
})

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
})

EXPORT(load_gauge_field, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_gf, &p_path)) {
    return NULL;
  }
  GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  const std::string path = py_convert_data<std::string>(p_path);
  const long ret = load_gauge_field(gf, path);
  return py_convert(ret);
})

EXPORT(gf_twist_boundary_at_boundary, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  double lmom = 0.0;
  int mu = -1;
  if (!PyArg_ParseTuple(args, "Odi", &p_gf, &lmom, &mu)) {
    return NULL;
  }
  GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  qassert(0 <= mu and mu < 4);
  twist_boundary_at_boundary(gf, lmom, mu);
  Py_RETURN_NONE;
})

EXPORT(save_gauge_transform_cps, {
  using namespace qlat;
  PyObject* p_gt = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_gt, &p_path)) {
    return NULL;
  }
  const GaugeTransform& gt = py_convert_type<GaugeTransform>(p_gt);
  const std::string path = py_convert_data<std::string>(p_path);
  const long ret = save_gauge_transform_cps(gt, path);
  return py_convert(ret);
})

EXPORT(load_gauge_transform_cps, {
  using namespace qlat;
  PyObject* p_gt = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_gt, &p_path)) {
    return NULL;
  }
  GaugeTransform& gt = py_convert_type<GaugeTransform>(p_gt);
  const std::string path = py_convert_data<std::string>(p_path);
  const long ret = load_gauge_transform_cps(gt, path);
  return py_convert(ret);
})

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
})

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
})

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
})

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
  const FieldSelection& fsel = py_convert_type<FieldSelection>(p_prop0, "fsel");
  prop_apply_gauge_transformation(prop, prop0, gt, fsel);
  Py_RETURN_NONE;
})

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
  const PointSelection& psel = py_convert_type<PointSelection>(p_prop0, "psel");
  prop_apply_gauge_transformation(prop, prop0, gt, psel);
  Py_RETURN_NONE;
})

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
})
