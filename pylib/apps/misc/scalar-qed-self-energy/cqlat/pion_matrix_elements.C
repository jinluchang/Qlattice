#include "lib.h"

#include <qlat/qlat.h>

namespace qlat
{  //

inline double r_scaling_factor() { return 5.0; }

inline LatData mk_four_point_em_table(const Coordinate& total_site, const int n_types)
{
  LatData ld;
  ld.info.push_back(lat_dim_number("type", 0, n_types - 1));
  ld.info.push_back(lat_dim_number("t", 0, total_site[3] / 2));
  ld.info.push_back(
      lat_dim_number("r", 0,
                     (long)std::ceil(1.0 + r_scaling_factor() * sqrt(3.0) *
                                               (double)total_site[0] / 2.0)));
  ld.info.push_back(
      lat_dim_string("em", make_array<std::string>("mm", "tt", "ii", "xx")));
  ld.info.push_back(lat_dim_re_im());
  lat_data_alloc(ld);
  set_zero(ld);
  return ld;
}

inline void acc_four_point_func_em(LatData& ld,
                                   const Field<Complex>& hf,
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
      const Vector<Complex> hv = hf.get_elems_const(xl);
      std::array<Complex, 4> vt;
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
      Vector<Complex> v1 =
          lat_data_complex_get(ld, make_array<int>(type, t, r1));
      Vector<Complex> v2 =
          lat_data_complex_get(ld, make_array<int>(type, t, r2));
      for (int m = 0; m < 4; ++m) {
        v1[m] += coef1 * vt[m];
        v2[m] += coef2 * vt[m];
      }
    }
  }
}

inline void partial_sum_r_four_point_func_em(LatData& ld)
{
  TIMER("partial_sum_r_four_point_func_em");
  for (int type = 0; type < ld.info[0].size; ++type) {
    for (int t = 0; t < ld.info[1].size; ++t) {
      for (int r = 1; r < ld.info[2].size; ++r) {
        const Vector<Complex> v1 =
            lat_data_complex_get(ld, make_array<int>(type, t, r - 1));
        Vector<Complex> v2 =
            lat_data_complex_get(ld, make_array<int>(type, t, r));
        qassert(v1.size() == ld.info[3].size);
        qassert(v2.size() == ld.info[3].size);
        for (int em = 0; em < ld.info[3].size; ++em) {
          v2[em] += v1[em];
        }
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
  } else if ("pole" == tag) {
    tag_index = 1;
  } else if ("0-pole" == tag) {
    tag_index = -1;
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
    double f_pi_sq = 0.0;
    if (0 == tag_index) {
      f_pi_sq = 1.0;
    } else if (1 == tag_index) {
      f_pi_sq = sqr(1.0 / (1.0 + sqr(r_pi) / 6.0 * s2));
    } else if (-1 == tag_index) {
      f_pi_sq = sqr(sqr(r_pi) / 6.0 * s2 / (1.0 + sqr(r_pi) / 6.0 * s2));
    }
    Vector<Complex> v = f.get_elems(kl);
    set_zero(v);
    if (kg == Coordinate()) {
      int mu_nu = 3 * 4 + 3;
      v[mu_nu] = total_site[3];
      v[mu_nu] *= f_pi_sq;
      continue;
    }
    array<Complex, 4> p_p_q, p_m_q, tp_p_q, tp_m_q;
    set_zero(p_p_q);
    set_zero(p_m_q);
    set_zero(tp_p_q);
    set_zero(tp_m_q);
    const double mh = 2 * std::sinh(pion_mass / 2);
    p_p_q[0] = ii * mh;
    p_m_q[0] = ii * mh;
    tp_p_q[0] = 2.0 * ii * mh;
    tp_m_q[0] = 2.0 * ii * mh;
    Complex s_p = sqr(mh);
    Complex s_m = sqr(mh);
    for (int mu = 0; mu < 4; ++mu) {
      p_p_q[mu] += ks[mu];
      p_m_q[mu] -= ks[mu];
      tp_p_q[mu] += ks[mu];
      tp_m_q[mu] -= ks[mu];
      s_p += sqr(p_p_q[mu]);
      s_m += sqr(p_m_q[mu]);
    }
    for (int mu = 0; mu < 4; ++mu) {
      for (int nu = 0; nu < 4; ++nu) {
        int mu_nu = mu * 4 + nu;
        if (mu == nu) {
          v[mu_nu] += 2.0;
        }
        v[mu_nu] -= tp_p_q[mu] * tp_p_q[nu] / s_p;
        v[mu_nu] -= tp_m_q[mu] * tp_m_q[nu] / s_m;
        v[mu_nu] *= f_pi_sq;
      }
    }
  });
  Py_RETURN_NONE;
});

EXPORT(acc_four_point_func_em, {
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
  pqassert(pf.ctype == "Complex");
  Field<Complex>& f = *(Field<Complex>*)pf.cdata;
  qlat::acc_four_point_func_em(ld, f, type, r_scaling_factor);
  Py_RETURN_NONE;
});

EXPORT(partial_sum_r_four_point_func_em, {
  using namespace qlat;
  PyObject* p_ld = NULL;
  if (!PyArg_ParseTuple(args, "O", &p_ld)) {
    return NULL;
  }
  LatData& ld = py_convert_type<LatData>(p_ld);
  qlat::partial_sum_r_four_point_func_em(ld);
  Py_RETURN_NONE;
});
