#include "lib.h"
#include <qlat/vector_utils/utils_smear_vecs.h>

EXPORT(gf_ape_smear, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  PyObject* p_gf0 = NULL;
  double alpha = 0.0;
  long steps = 1;
  if (!PyArg_ParseTuple(args, "OOd|l", &p_gf, &p_gf0, &alpha, &steps)) {
    return NULL;
  }
  GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  const GaugeField& gf0 = py_convert_type<GaugeField>(p_gf0);
  gf_ape_smear(gf, gf0, alpha, steps);
  Py_RETURN_NONE;
})

EXPORT(gf_spatial_ape_smear, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  PyObject* p_gf0 = NULL;
  // used value: alpha = 0.5, steps = 30
  double alpha = 0.0;
  long steps = 1;
  if (!PyArg_ParseTuple(args, "OOd|l", &p_gf, &p_gf0, &alpha, &steps)) {
    return NULL;
  }
  GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  const GaugeField& gf0 = py_convert_type<GaugeField>(p_gf0);
  gf_spatial_ape_smear(gf, gf0, alpha, steps);
  Py_RETURN_NONE;
})

EXPORT(gf_hyp_smear, {
  using namespace qlat;
  PyObject* p_gf = NULL;
  PyObject* p_gf0 = NULL;
  // values in paper is 0.75 0.6 0.3
  // 10.1103/PhysRevD.64.034504 Eq(4)
  double alpha1 = 0.0;
  double alpha2 = 0.0;
  double alpha3 = 0.0;
  if (!PyArg_ParseTuple(args, "OOddd", &p_gf, &p_gf0, &alpha1, &alpha2,
                        &alpha3)) {
    return NULL;
  }
  GaugeField& gf = py_convert_type<GaugeField>(p_gf);
  const GaugeField& gf0 = py_convert_type<GaugeField>(p_gf0);
  gf_hyp_smear(gf, gf0, alpha1, alpha2, alpha3);
  Py_RETURN_NONE;
})

EXPORT(prop_smear, {
  using namespace qlat;
  // prop is of normal size
  PyObject* p_prop = NULL;
  // gf1 is left_expanded and refreshed
  // set_left_expanded_gauge_field(gf1, gf)
  PyObject* p_gf1 = NULL;
  // smear params:
  // 24D: coef = 0.9375, step = 10
  // 48I: coef = 0.9375, step = 29
  // mom ratio = 0.5
  double coef = 0;
  int step = 0;
  int mode_smear = 0;
  // momentum in lattice unit 1/a (not unit of lattice momentum 2 pi / L / a)
  PyObject* p_mom = NULL;
  bool smear_in_time_dir = false;
  if (!PyArg_ParseTuple(args, "OOdi|Obi", &p_prop, &p_gf1, &coef, &step, &p_mom,
                        &smear_in_time_dir, &mode_smear)) {
    return NULL;
  }
  Propagator4d& prop = py_convert_type<Propagator4d>(p_prop);
  const GaugeField& gf1 = py_convert_type<GaugeField>(p_gf1);
  const CoordinateD mom =
      NULL == p_mom ? CoordinateD() : py_convert_data<CoordinateD>(p_mom);
  if(mode_smear == 0){
    smear_propagator(prop, gf1, coef, step, mom, smear_in_time_dir);}
  if(mode_smear >= 1){
    smear_propagator_qlat_convension(prop, gf1, coef, step, mom, smear_in_time_dir, mode_smear);
  }
  Py_RETURN_NONE;
})
