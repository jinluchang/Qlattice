#include "lib.h"

EXPORT(mk_inverter_domain_wall, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  PyObject* p_fa = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_gf, &p_fa)) {
    return NULL;
  }
  InverterDomainWall* pinv = new InverterDomainWall();
  InverterDomainWall& inv = *pinv;
  const GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  const FermionAction& fa = py_convert_type<FermionAction>(p_fa);
  setup_inverter(inv, gf, fa);
  return py_convert((void*)pinv);
})

EXPORT(free_inverter_domain_wall, {
  using namespace qlat;
  return free_obj<InverterDomainWall>(args);
})

EXPORT(invert_inverter_domain_wall, {
  using namespace qlat;
  PyObject* p_prop_sol = NULL;
  PyObject* p_prop_src = NULL;
  PyObject* p_inv = NULL;
  if (!PyArg_ParseTuple(args, "OOO", &p_prop_sol, &p_prop_src, &p_inv)) {
    return NULL;
  }
  Propagator4d& prop_sol = py_convert_type<Propagator4d>(p_prop_sol);
  const Propagator4d& prop_src = py_convert_type<Propagator4d>(p_prop_src);
  const InverterDomainWall& inv = py_convert_type<InverterDomainWall>(p_inv);
  invert(prop_sol, prop_src, inv);
  Py_RETURN_NONE;
})

EXPORT(get_stop_rsd_inverter_domain_wall, {
  using namespace qlat;
  PyObject* p_inv = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_inv)) {
    return NULL;
  }
  const InverterDomainWall& inv = py_convert_type<InverterDomainWall>(p_inv);
  return py_convert(inv.stop_rsd());
})

EXPORT(set_stop_rsd_inverter_domain_wall, {
  using namespace qlat;
  PyObject* p_inv = NULL;
  double stop_rsd = 1e-8;
  if (!PyArg_ParseTuple(args, "Od", &p_inv, &stop_rsd)) {
    return NULL;
  }
  InverterDomainWall& inv = py_convert_type<InverterDomainWall>(p_inv);
  inv.stop_rsd() = stop_rsd;
  Py_RETURN_NONE;
})

EXPORT(get_max_num_iter_inverter_domain_wall, {
  using namespace qlat;
  PyObject* p_inv = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_inv)) {
    return NULL;
  }
  const InverterDomainWall& inv = py_convert_type<InverterDomainWall>(p_inv);
  return py_convert(inv.max_num_iter());
})

EXPORT(set_max_num_iter_inverter_domain_wall, {
  using namespace qlat;
  PyObject* p_inv = NULL;
  long max_num_iter = 200;
  if (!PyArg_ParseTuple(args, "Ol", &p_inv, &max_num_iter)) {
    return NULL;
  }
  InverterDomainWall& inv = py_convert_type<InverterDomainWall>(p_inv);
  inv.max_num_iter() = max_num_iter;
  Py_RETURN_NONE;
})

EXPORT(get_max_mixed_precision_cycle_inverter_domain_wall, {
  using namespace qlat;
  PyObject* p_inv = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_inv)) {
    return NULL;
  }
  const InverterDomainWall& inv = py_convert_type<InverterDomainWall>(p_inv);
  return py_convert(inv.max_mixed_precision_cycle());
})

EXPORT(set_max_mixed_precision_cycle_inverter_domain_wall, {
  using namespace qlat;
  PyObject* p_inv = NULL;
  long max_mixed_precision_cycle = 300;
  if (!PyArg_ParseTuple(args, "Ol", &p_inv, &max_mixed_precision_cycle)) {
    return NULL;
  }
  InverterDomainWall& inv = py_convert_type<InverterDomainWall>(p_inv);
  inv.max_mixed_precision_cycle() = max_mixed_precision_cycle;
  Py_RETURN_NONE;
})
