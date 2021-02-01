#include "convert.h"
#include "dispatch.h"
#include "exceptions.h"

EXPORT(gf_show_info, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  void* pfield = NULL;
  if (!PyArg_ParseTuple(args, "Ol", &p_ctype, &pfield)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  pqassert(ctype == "ColorMatrix");
  const GaugeField& gf = *(GaugeField*)pfield;
  gf_show_info(gf);
  Py_RETURN_NONE;
});

EXPORT(set_g_rand_color_matrix_field, {
  using namespace qlat;
  PyObject* p_ctype = NULL;
  void* pfield = NULL;
  void* prng = NULL;
  double sigma = 1.0;
  int n_step = 1;
  if (!PyArg_ParseTuple(args, "Olld|i", &p_ctype, &pfield, &prng, &sigma,
                        &n_step)) {
    return NULL;
  }
  std::string ctype;
  py_convert(ctype, p_ctype);
  pqassert(ctype == "ColorMatrix");
  Field<ColorMatrix>& field = *(Field<ColorMatrix>*)pfield;
  const RngState& rng = *(RngState*)prng;
  set_g_rand_color_matrix_field(field, rng, sigma, n_step);
  Py_RETURN_NONE;
});
