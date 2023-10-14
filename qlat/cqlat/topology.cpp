#include "lib.h"

EXPORT(gf_topology_field_clf, {
  // Use the basic gf_clover_leaf_field
  // NOT using 5 loop improved definition
  using namespace qlat;
  PyObject* p_topf = NULL;
  PyObject* p_gf = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_topf, &p_gf)) {
    return NULL;
  }
  const GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  qassert(py_get_ctype(p_topf) == "Double");
  FieldM<double, 1>& topf = py_convert_type<FieldM<double, 1> >(p_topf);
  qassert(topf.geo().multiplicity == 1);
  clf_topology_field(topf, gf);
  Py_RETURN_NONE;
})

EXPORT(gf_topology_field, {
  // using the 5 loop improved definition
  // https://arxiv.org/pdf/hep-lat/9701012v2.pdf
  using namespace qlat;
  PyObject* p_topf = NULL;
  PyObject* p_gf = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_topf, &p_gf)) {
    return NULL;
  }
  const GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  qassert(py_get_ctype(p_topf) == "Double");
  FieldM<double, 1>& topf = py_convert_type<FieldM<double, 1> >(p_topf);
  qassert(topf.geo().multiplicity == 1);
  clf_topology_field_5(topf, gf);
  Py_RETURN_NONE;
})

EXPORT(gf_topology_terms_field, {
  // using the 5 loop improved definition (each term's contribution individually)
  // https://arxiv.org/pdf/hep-lat/9701012v2.pdf
  using namespace qlat;
  PyObject* p_topf = NULL;
  PyObject* p_gf = NULL;
  if (!PyArg_ParseTuple(args, "OO", &p_topf, &p_gf)) {
    return NULL;
  }
  const GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  qassert(py_get_ctype(p_topf) == "Double");
  FieldM<double, 5>& topf = py_convert_type<FieldM<double, 5> >(p_topf);
  qassert(topf.geo().multiplicity == 5);
  clf_topology_field_5_terms(topf, gf);
  Py_RETURN_NONE;
})
