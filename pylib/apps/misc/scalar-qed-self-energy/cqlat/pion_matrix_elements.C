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
  int tag_index = 0;
  if ("" == tag) {
    tag_index = 0;
  } else {
    qassert(false);
  }
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
    array<Complex, 4> tp_p_q, tp_m_q;
    set_zero(tp_p_q);
    set_zero(tp_m_q);
    tp_p_q[0] = 2 * ii * pion_mass;
    tp_m_q[0] = 2 * ii * pion_mass;
    for (int mu = 0; mu < 4; ++mu) {
      tp_p_q[mu] += ks[mu];
      tp_m_q[mu] -= ks[mu];
    }
    for (int mu = 0; mu < 4; ++mu) {
      for (int nu = 0; nu < 4; ++nu) {
        int mu_nu = mu * 4 + nu;
        if (mu == nu) {
          v[mu_nu] += 2.0;
        }

      }
    }
    // TODO
  });
  Py_RETURN_NONE;
});
