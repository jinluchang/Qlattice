#include "lib.h"

EXPORT(get_gm_force_magnitudes, {
  using namespace qlat;
  PyObject* p_gm_force = NULL;
  int n_elems = 0;
  if (!PyArg_ParseTuple(args, "Oi", &p_gm_force, &n_elems)) {
    return NULL;
  }
  const GaugeMomentum& gm_force = py_convert_type<GaugeMomentum>(p_gm_force);
  std::vector<double> ret = get_gm_force_magnitudes(gm_force, n_elems);
  return py_convert(ret);
})

EXPORT(display_gm_force_magnitudes, {
  using namespace qlat;
  PyObject* p_gm_force = NULL;
  int n_elems = 0;
  if (!PyArg_ParseTuple(args, "Oi", &p_gm_force, &n_elems)) {
    return NULL;
  }
  const GaugeMomentum& gm_force = py_convert_type<GaugeMomentum>(p_gm_force);
  display_gm_force_magnitudes(gm_force, n_elems);
  Py_RETURN_NONE;
})

EXPORT(save_gm_force_magnitudes_list, {
  using namespace qlat;
  PyObject* p_fn = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_fn)) {
    return NULL;
  }
  std::string fn;
  py_convert(fn, p_fn);
  save_gm_force_magnitudes_list(fn);
  Py_RETURN_NONE;
})

EXPORT(display_gauge_field_info_table_with_wilson_flow, {
  using namespace qlat;
  PyObject* p_fn_gf_info = NULL;
  PyObject* p_fn_wilson_flow_energy = NULL;
  PyObject* p_gf = NULL;
  double flow_time = 0.0;
  int flow_steps = 0.0;
  int steps = 0.0;
  double c1 = 0.0;
  if (!PyArg_ParseTuple(args, "OOOdii|d", &p_fn_gf_info,
                        &p_fn_wilson_flow_energy, &p_gf, &flow_time,
                        &flow_steps, &steps, &c1)) {
    return NULL;
  }
  std::string fn_gf_info, fn_wilson_flow_energy;
  py_convert(fn_gf_info, p_fn_gf_info);
  py_convert(fn_wilson_flow_energy, p_fn_wilson_flow_energy);
  const GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  display_gauge_field_info_table_with_wilson_flow(
      fn_gf_info, fn_wilson_flow_energy, gf, flow_time, flow_steps, steps, c1);
  Py_RETURN_NONE;
})
