#include "convert.h"
#include "dispatch.h"
#include "exceptions.h"

EXPORT(set_rand_gauge_momentum, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  GaugeMomentum* pfield = NULL;
  double sigma = 1.0;
  RngState* prng = NULL;
  if (!PyArg_ParseTuple(args, "Oldl", &p_ctype, &pfield, &sigma, &prng)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  pqassert(ctype == "ColorMatrix");
  GaugeMomentum& gm = *pfield;
  const RngState& rs = *prng;
  set_rand_gauge_momentum(gm, sigma, rs);
  Py_RETURN_NONE;
});

EXPORT(gm_hamilton_node, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  GaugeMomentum* pfield = NULL;
  if (!PyArg_ParseTuple(args, "Ol", &p_ctype, &pfield)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  pqassert(ctype == "ColorMatrix");
  const GaugeMomentum& gm = *pfield;
  const double ret = gm_hamilton_node(gm);
  return py_convert(ret);
});

EXPORT(gf_hamilton_node, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  GaugeField* pfield = NULL;
  GaugeAction* pga = NULL;
  if (!PyArg_ParseTuple(args, "Oll", &p_ctype, &pfield, &pga)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  pqassert(ctype == "ColorMatrix");
  const GaugeField& gf = *pfield;
  const GaugeAction& ga = *pga;
  const double ret = gf_hamilton_node(gf, ga);
  return py_convert(ret);
});

EXPORT(gf_evolve, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  GaugeField* pfield = NULL;
  PyObject* p_ctype_1 = NULL;
  GaugeMomentum* pfield_1 = NULL;
  double step_size = 0.0;
  if (!PyArg_ParseTuple(args, "OlOld", &p_ctype, &pfield, &p_ctype_1, &pfield_1,
                        &step_size)) {
    return NULL;
  }
  std::string ctype, ctype_1;
  py_convert(ctype, p_ctype);
  py_convert(ctype_1, p_ctype_1);
  pqassert(ctype == "ColorMatrix");
  pqassert(ctype_1 == "ColorMatrix");
  GaugeField& gf = *pfield;
  const GaugeMomentum& gm = *pfield_1;
  gf_evolve(gf, gm, step_size);
  Py_RETURN_NONE;
});

EXPORT(set_gm_force, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  GaugeMomentum* pfield = NULL;
  PyObject* p_ctype_1 = NULL;
  GaugeField* pfield_1 = NULL;
  GaugeAction* pga = NULL;
  if (!PyArg_ParseTuple(args, "OlOll", &p_ctype, &pfield, &p_ctype_1, &pfield_1,
                        &pga)) {
    return NULL;
  }
  std::string ctype, ctype_1;
  py_convert(ctype, p_ctype);
  py_convert(ctype_1, p_ctype_1);
  pqassert(ctype == "ColorMatrix");
  pqassert(ctype_1 == "ColorMatrix");
  GaugeMomentum& gm_force = *pfield;
  const GaugeField& gf = *pfield_1;
  const GaugeAction& ga = *pga;
  set_gm_force(gm_force, gf, ga);
  Py_RETURN_NONE;
});
