#include "lib.h"

EXPORT(gf_wilson_line_no_comm, {  // tested: cqlat-qcd
  using namespace qlat;
  PyObject* p_wilson_line_field = NULL;
  Int wilson_line_field_m = 0;
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
  const std::vector<int> path = py_convert_data<std::vector<int>>(p_path);
  if (p_path_n != NULL) {
    const std::vector<int> path_n = py_convert_data<std::vector<int>>(p_path_n);
    gf_wilson_line_no_comm(wilson_line_field, wilson_line_field_m, gf_ext, path,
                           path_n);
  } else {
    gf_wilson_line_no_comm(wilson_line_field, wilson_line_field_m, gf_ext,
                           path);
  }
  Py_RETURN_NONE;
})

EXPORT(gf_twist_boundary_at_boundary, {  // tested: cqlat-qcd
  using namespace qlat;
  PyObject* p_gf = NULL;
  RealD lmom = 0.0;
  Int mu = -1;
  if (!PyArg_ParseTuple(args, "Odi", &p_gf, &lmom, &mu)) {
    return NULL;
  }
  GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  qassert(0 <= mu and mu < 4);
  twist_boundary_at_boundary(gf, lmom, mu);
  Py_RETURN_NONE;
})

EXPORT(save_gauge_transform_cps, {  // tested: cqlat-qcd
  using namespace qlat;
  PyObject* p_gt = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_gt, &p_path)) {
    return NULL;
  }
  const GaugeTransform& gt = py_convert_type<GaugeTransform>(p_gt);
  const std::string path = py_convert_data<std::string>(p_path);
  const Long ret = save_gauge_transform_cps(gt, path);
  return py_convert(ret);
})

EXPORT(load_gauge_transform_cps, {  // tested: cqlat-qcd
  using namespace qlat;
  PyObject* p_gt = NULL;
  PyObject* p_path = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_gt, &p_path)) {
    return NULL;
  }
  GaugeTransform& gt = py_convert_type<GaugeTransform>(p_gt);
  const std::string path = py_convert_data<std::string>(p_path);
  const Long ret = load_gauge_transform_cps(gt, path);
  return py_convert(ret);
})
