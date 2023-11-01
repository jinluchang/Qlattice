#include "lib.h"

#include <qlat/qlat.h>

namespace qlat
{  //

inline void acc_four_point_func_em(LatData& ld,
                                   const Field<ComplexD>& hf,
                                   const int type,
                                   const double r_scaling_factor)
{
  TIMER("acc_four_point_func_em");
  const Geometry geo = hf.geo();
  const int multiplicity = geo.multiplicity;
  int ndir;
  if (16 == multiplicity) {
    ndir = 4;
  } else if (64 == multiplicity) {
    ndir = 8;
  } else {
    qassert(false);
  }
  const Coordinate total_site = geo.total_site();
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    const Coordinate xg = geo.coordinate_g_from_l(xl);
    const Coordinate xrel = smod(xg, total_site);
    const CoordinateD xreld = smod_sym(xg, total_site);
    const int t = std::abs(xrel[3]);
    const int sqr_x = sqr(xrel[0]) + sqr(xrel[1]) + sqr(xrel[2]);
    const double rd = r_scaling_factor * std::sqrt((double)sqr_x);
    const int r1 = std::floor(rd);
    const int r2 = r1 + 1;
    const double coef1 = r2 - rd;
    const double coef2 = rd - r1;
    qassert(r2 < ld.info[2].size);
    if (t < ld.info[1].size) {
      const Vector<ComplexD> hv = hf.get_elems_const(xl);
      std::array<ComplexD, 4> vt;
      vt[0] = 0.0;
      for (int mu = 0; mu < 4; ++mu) {
        vt[0] += hv[ndir * mu + mu];
      }
      vt[1] = hv[ndir * 3 + 3];
      vt[2] = 0.0;
      for (int i = 0; i < 3; ++i) {
        vt[2] += hv[ndir * i + i];
      }
      vt[3] = 0.0;
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          vt[3] += xreld[i] * xreld[j] * hv[ndir * i + j];
        }
      }
      Vector<ComplexD> v1 =
          lat_data_complex_get(ld, make_array<int>(type, t, r1));
      Vector<ComplexD> v2 =
          lat_data_complex_get(ld, make_array<int>(type, t, r2));
      for (int m = 0; m < 4; ++m) {
        v1[m] += coef1 * vt[m];
        v2[m] += coef2 * vt[m];
      }
    }
  }
}

}  // namespace qlat

EXPORT(set_pion_four_point_mom_field, {
  TIMER("set_pion_four_point_mom_field");
  using namespace qlat;
  PyObject* p_field = NULL;
  double pion_mass = 0.0;
  PyObject* p_tag = NULL;
  double r_pi = 0.0;
  if (!PyArg_ParseTuple(args, "OdOd", &p_field, &pion_mass, &p_tag, &r_pi)) {
    return NULL;
  }
  PyField pf = py_convert_field(p_field);
  pqassert(pf.ctype == "ComplexD");
  Field<ComplexD>& f = *(Field<Complex>*)pf.cdata;
  const Geometry& geo = f.geo();
  pqassert(f.geo().multiplicity == 16);
  std::string tag;
  py_convert(tag, p_tag);
  const Coordinate total_site = geo.total_site();
  const CoordinateD momtwist;
  int tag_index = 0;
  if ("" == tag) {
    tag_index = 0;
  } else if ("pole" == tag) {
    tag_index = 1;
  } else if ("pole_p" == tag) {
    tag_index = 2;
  } else if ("linear" == tag) {
    tag_index = 3;
  } else if ("0-pole" == tag) {
    tag_index = -1;
  } else if ("0-linear" == tag) {
    tag_index = -2;
  } else {
    qassert(false);
  }
  qacc_for(index, geo.local_volume(), {
    const Coordinate kl = geo.coordinate_from_index(index);
    const Coordinate kg = geo.coordinate_g_from_l(kl);
    const Coordinate kg_rel = smod(kg, total_site);
    const CoordinateD kg_rel_sym = smod_sym(kg, total_site);
    CoordinateD kk, ks;
    double s2 = 0.0;
    CoordinateD kk_sym, ks_sym;
    for (int i = 0; i < DIMN; i++) {
      kk[i] = 2.0 * PI * (kg_rel[i] + momtwist[i]) / (double)total_site[i];
      ks[i] = 2.0 * std::sin(kk[i] / 2.0);
      s2 += sqr(ks[i]);
      kk_sym[i] = 2.0 * PI * (kg_rel_sym[i] + momtwist[i]) / (double)total_site[i];
      ks_sym[i] = 2.0 * std::sin(kk_sym[i] / 2.0);
    }
    double f_pi_sq = 0.0;
    if (0 == tag_index) {
      f_pi_sq = 1.0;
    } else if (1 == tag_index) {
      f_pi_sq = sqr(1.0 / (1.0 + sqr(r_pi) / 6.0 * s2));
    } else if (2 == tag_index) {
      f_pi_sq = sqr(1.0 / (1.0 + sqr(r_pi) / 6.0 * s2 + 0.01 * sqr(sqr(r_pi) / 6.0 * s2)));
    } else if (3 == tag_index) {
      f_pi_sq = sqr(1.0 - sqr(r_pi) / 6.0 * s2);
    } else if (-1 == tag_index) {
      f_pi_sq = sqr(sqr(r_pi) / 6.0 * s2 / (1.0 + sqr(r_pi) / 6.0 * s2));
    } else if (-2 == tag_index) {
      f_pi_sq = sqr(sqr(r_pi) / 6.0 * s2);
    }
    Vector<ComplexD> v = f.get_elems(kl);
    set_zero(v);
    for (int mu = 0; mu < 4; ++mu) {
      for (int nu = 0; nu < 4; ++nu) {
        const int mu_nu = mu * 4 + nu;
        if (mu == nu) {
          v[mu_nu] += 2.0;
        }
      }
    }
    const bool is_only_contact = false;
    if (is_only_contact) {
      // no need to do anything if we only include the contact term
    } else if (kg == Coordinate()) {
      const int mu = 3;
      const int nu = 3;
      const int mu_nu = mu * 4 + nu;
      v[mu_nu] = total_site[3] * 2 * pion_mass;
    } else {
      array<ComplexD, 4> p_p_q, p_m_q, tp_p_q, tp_m_q, tp_p_q_sym, tp_m_q_sym;
      set_zero(p_p_q);
      set_zero(p_m_q);
      set_zero(tp_p_q);
      set_zero(tp_m_q);
      set_zero(tp_p_q_sym);
      set_zero(tp_m_q_sym);
      const double mh = 2 * std::sinh(pion_mass / 2);
      p_p_q[3] = ii * mh;
      p_m_q[3] = ii * mh;
      tp_p_q[3] = 2.0 * ii * mh;
      tp_m_q[3] = 2.0 * ii * mh;
      tp_p_q_sym[3] = 2.0 * ii * mh;
      tp_m_q_sym[3] = 2.0 * ii * mh;
      ComplexD s_p = sqr(mh);
      ComplexD s_m = sqr(mh);
      for (int mu = 0; mu < 4; ++mu) {
        p_p_q[mu] += ks[mu];
        p_m_q[mu] -= ks[mu];
        tp_p_q[mu] += ks[mu];
        tp_m_q[mu] -= ks[mu];
        tp_p_q_sym[mu] += ks_sym[mu];
        tp_m_q_sym[mu] -= ks_sym[mu];
        s_p += sqr(p_p_q[mu]);
        s_m += sqr(p_m_q[mu]);
      }
      for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
          const int mu_nu = mu * 4 + nu;
          if (mu == nu) {
            v[mu_nu] -= tp_p_q[mu] * tp_p_q[nu] / s_p;
            v[mu_nu] -= tp_m_q[mu] * tp_m_q[nu] / s_m;
          } else {
            v[mu_nu] -= tp_p_q_sym[mu] * tp_p_q_sym[nu] / s_p;
            v[mu_nu] -= tp_m_q_sym[mu] * tp_m_q_sym[nu] / s_m;
          }
        }
      }
    }
    for (int mu = 0; mu < 4; ++mu) {
      for (int nu = 0; nu < 4; ++nu) {
        const int mu_nu = mu * 4 + nu;
        v[mu_nu] *= f_pi_sq / (2 * pion_mass);
      }
    }
  });
  Py_RETURN_NONE;
});

EXPORT(acc_four_point_func_em, {
  // There is no glb_sum in this function
  // need to perform glb_sum explicitly after calling this function
  using namespace qlat;
  PyObject* p_ld = NULL;
  PyObject* p_field = NULL;
  int type = -1;
  double r_scaling_factor = 0.0;
  if (!PyArg_ParseTuple(args, "OOid", &p_ld, &p_field, &type, &r_scaling_factor)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  PyField pf = py_convert_field(p_field);
  pqassert(pf.ctype == "ComplexD");
  Field<ComplexD>& f = *(Field<Complex>*)pf.cdata;
  qlat::acc_four_point_func_em(ld, f, type, r_scaling_factor);
  Py_RETURN_NONE;
});
