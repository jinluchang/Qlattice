#include "lib.h"

EXPORT(free_inverter_domain_wall, {
  using namespace qlat;
  return free_obj<InverterDomainWall>(args);
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
  RealD stop_rsd = 1e-8;
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
  Long max_num_iter = 200;
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
  Long max_mixed_precision_cycle = 300;
  if (!PyArg_ParseTuple(args, "Ol", &p_inv, &max_mixed_precision_cycle)) {
    return NULL;
  }
  InverterDomainWall& inv = py_convert_type<InverterDomainWall>(p_inv);
  inv.max_mixed_precision_cycle() = max_mixed_precision_cycle;
  Py_RETURN_NONE;
})
