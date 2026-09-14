#include "lib.h"

EXPORT(free_inverter_domain_wall, {  // tested: cqlat-dwf-inverter
  using namespace qlat;
  return free_obj<InverterDomainWall>(args);
})

EXPORT(get_stop_rsd_inverter_domain_wall, {  // tested: cqlat-dwf-inverter
  using namespace qlat;
  PyObject* p_inv = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_inv)) {
    return NULL;
  }
  const InverterDomainWall& inv = py_convert_type<InverterDomainWall>(p_inv);
  return py_convert(inv.stop_rsd());
})

EXPORT(set_stop_rsd_inverter_domain_wall, {  // tested: cqlat-dwf-inverter
  using namespace qlat;
  PyObject* p_inv = NULL;
  // the default value (from InverterParams::init()) is only in effect until
  // this setter is called; the value argument is required
  RealD stop_rsd = 0.0;
  if (!PyArg_ParseTuple(args, "Od", &p_inv, &stop_rsd)) {
    return NULL;
  }
  InverterDomainWall& inv = py_convert_type<InverterDomainWall>(p_inv);
  inv.stop_rsd() = stop_rsd;
  Py_RETURN_NONE;
})

EXPORT(get_max_num_iter_inverter_domain_wall, {  // tested: cqlat-dwf-inverter
  using namespace qlat;
  PyObject* p_inv = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_inv)) {
    return NULL;
  }
  const InverterDomainWall& inv = py_convert_type<InverterDomainWall>(p_inv);
  return py_convert(inv.max_num_iter());
})

EXPORT(set_max_num_iter_inverter_domain_wall, {  // tested: cqlat-dwf-inverter
  using namespace qlat;
  PyObject* p_inv = NULL;
  // required argument (no default is declared here; InverterParams::init()
  // provides the initial value); "L" matches Long == int64_t on every platform
  Long max_num_iter = 0;
  if (!PyArg_ParseTuple(args, "OL", &p_inv, &max_num_iter)) {
    return NULL;
  }
  InverterDomainWall& inv = py_convert_type<InverterDomainWall>(p_inv);
  inv.max_num_iter() = max_num_iter;
  Py_RETURN_NONE;
})

EXPORT(get_max_mixed_precision_cycle_inverter_domain_wall, {  // tested: cqlat-dwf-inverter
  using namespace qlat;
  PyObject* p_inv = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_inv)) {
    return NULL;
  }
  const InverterDomainWall& inv = py_convert_type<InverterDomainWall>(p_inv);
  return py_convert(inv.max_mixed_precision_cycle());
})

EXPORT(set_max_mixed_precision_cycle_inverter_domain_wall, {  // tested: cqlat-dwf-inverter
  using namespace qlat;
  PyObject* p_inv = NULL;
  // required argument (no default is declared here; InverterParams::init()
  // provides the initial value); "L" matches Long == int64_t on every platform
  Long max_mixed_precision_cycle = 0;
  if (!PyArg_ParseTuple(args, "OL", &p_inv, &max_mixed_precision_cycle)) {
    return NULL;
  }
  InverterDomainWall& inv = py_convert_type<InverterDomainWall>(p_inv);
  inv.max_mixed_precision_cycle() = max_mixed_precision_cycle;
  Py_RETURN_NONE;
})
