#include "lib.h"

#include <qlat/qlat.h>

EXPORT(set_pion_four_point_mom_field, {
  using namespace qlat;
  PyObject* p_field = NULL;
  double pion_mass = 0.0;
  PyObject* p_tag = NULL;
  double r_pi = 0.0;
  if (!PyArg_ParseTuple(args, "OdOd", &p_field, &pion_mass, &p_tag, &r_pi)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  pqassert(pf.ctype == "Complex");
  Field<Complex>& f = *(Field<Complex>*)pf.cdata;
  const Geometry& geo = f.geo();
  pqassert(f.geo().multiplicity == 16);
  std::string tag;
  py_convert(tag, p_tag);
  const Coordinate total_site = geo.total_site();
  const CoordinateD momtwist;
  qacc_for(index, geo.local_volume(), {
    const Coordinate kl = geo.coordinate_from_index(index);
    Coordinate kg = geo.coordinate_g_from_l(kl);
    CoordinateD kk, ks;
    double s2 = 0.0;
    for (int i = 0; i < DIMN; i++) {
      kg[i] = smod(kg[i], total_site[i]);
      kk[i] = 2.0 * PI * (kg[i] + momtwist[i]) / (double)total_site[i];
      ks[i] = 2.0 * std::sin(kk[i] / 2.0);
      s2 += sqr(ks[i]);
    }
    Vector<Complex> v = f.get_elems(kl);
    set_zero(v);
    // TODO
  });
  Py_RETURN_NONE;
});


