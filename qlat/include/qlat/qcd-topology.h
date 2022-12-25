#pragma once

#include <qlat-utils/matrix.h>
#include <qlat/qcd-utils.h>
#include <qlat/qcd.h>

namespace qlat
{  //

struct CloverLeafField : FieldM<ColorMatrix, 6> {
};

qacc ColorMatrix gf_clover_leaf_no_comm(const GaugeField& gf1,
                                        const Coordinate& xl, const int mu,
                                        const int nu)
{
  ColorMatrix m;
  set_zero(m);
  m += gf_wilson_line_no_comm(gf1, xl,
                              make_array<int>(mu, nu, -mu - 1, -nu - 1));
  m += gf_wilson_line_no_comm(gf1, xl,
                              make_array<int>(-mu - 1, -nu - 1, mu, nu));
  m += gf_wilson_line_no_comm(gf1, xl,
                              make_array<int>(nu, -mu - 1, -nu - 1, mu));
  m += gf_wilson_line_no_comm(gf1, xl,
                              make_array<int>(-nu - 1, mu, nu, -mu - 1));
  return (Complex)0.25 * m;
}

qacc ColorMatrix gf_clover_leaf_m_n_no_comm(const GaugeField& gf1,
                                            const Coordinate& xl, const int mu,
                                            const int m, const int nu,
                                            const int n)
{
  ColorMatrix cm;
  set_zero(cm);
  cm +=
      gf_wilson_line_no_comm(gf1, xl, make_array<int>(mu, nu, -mu - 1, -nu - 1),
                             make_array<int>(m, n, m, n));
  cm +=
      gf_wilson_line_no_comm(gf1, xl, make_array<int>(-mu - 1, -nu - 1, mu, nu),
                             make_array<int>(m, n, m, n));
  cm +=
      gf_wilson_line_no_comm(gf1, xl, make_array<int>(nu, -mu - 1, -nu - 1, mu),
                             make_array<int>(n, m, n, m));
  cm +=
      gf_wilson_line_no_comm(gf1, xl, make_array<int>(-nu - 1, mu, nu, -mu - 1),
                             make_array<int>(n, m, n, m));
  return (Complex)0.25 * cm;
}

inline void gf_clover_leaf_field_no_comm(CloverLeafField& clf,
                                         const GaugeField& gf1)
// F_01, F_02, F_03, F_12, F_13, F_23
{
  TIMER("gf_clover_leaf_field_no_comm");
  const Geometry geo = geo_reform(gf1.geo(), 6, 0);
  clf.init(geo);
  qassert(is_matching_geo_mult(clf.geo(), geo));
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = clf.get_elems(xl);
    v[0] = gf_clover_leaf_no_comm(gf1, xl, 0, 1);
    v[1] = gf_clover_leaf_no_comm(gf1, xl, 0, 2);
    v[2] = gf_clover_leaf_no_comm(gf1, xl, 0, 3);
    v[3] = gf_clover_leaf_no_comm(gf1, xl, 1, 2);
    v[4] = gf_clover_leaf_no_comm(gf1, xl, 1, 3);
    v[5] = gf_clover_leaf_no_comm(gf1, xl, 2, 3);
  });
}

inline void gf_clover_leaf_field_m_n_no_comm(CloverLeafField& clf,
                                             const GaugeField& gf1, const int m,
                                             const int n)
// F_01, F_02, F_03, F_12, F_13, F_23
{
  TIMER("gf_clover_leaf_field_m_n_no_comm");
  const Geometry geo = geo_reform(gf1.geo(), 6, 0);
  clf.init(geo);
  qassert(is_matching_geo_mult(clf.geo(), geo));
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = clf.get_elems(xl);
    v[0] = gf_clover_leaf_m_n_no_comm(gf1, xl, 0, m, 1, n);
    v[1] = gf_clover_leaf_m_n_no_comm(gf1, xl, 0, m, 2, n);
    v[2] = gf_clover_leaf_m_n_no_comm(gf1, xl, 0, m, 3, n);
    v[3] = gf_clover_leaf_m_n_no_comm(gf1, xl, 1, m, 2, n);
    v[4] = gf_clover_leaf_m_n_no_comm(gf1, xl, 1, m, 3, n);
    v[5] = gf_clover_leaf_m_n_no_comm(gf1, xl, 2, m, 3, n);
    if (m != n) {
      v[0] += gf_clover_leaf_m_n_no_comm(gf1, xl, 0, n, 1, m);
      v[1] += gf_clover_leaf_m_n_no_comm(gf1, xl, 0, n, 2, m);
      v[2] += gf_clover_leaf_m_n_no_comm(gf1, xl, 0, n, 3, m);
      v[3] += gf_clover_leaf_m_n_no_comm(gf1, xl, 1, n, 2, m);
      v[4] += gf_clover_leaf_m_n_no_comm(gf1, xl, 1, n, 3, m);
      v[5] += gf_clover_leaf_m_n_no_comm(gf1, xl, 2, n, 3, m);
      for (int i = 0; i < 6; ++i) {
        v[i] *= 0.5;
      }
    }
  });
}

inline void gf_clover_leaf_field(CloverLeafField& clf, const GaugeField& gf)
{
  TIMER("gf_clover_leaf_field");
  GaugeField gf1;
  gf1.init(geo_resize(gf.geo(), 1));
  gf1 = gf;
  refresh_expanded(gf1);
  gf_clover_leaf_field_no_comm(clf, gf1);
}

inline void gf_clover_leaf_field_5(CloverLeafField& clf1, CloverLeafField& clf2,
                                   CloverLeafField& clf3, CloverLeafField& clf4,
                                   CloverLeafField& clf5, const GaugeField& gf)
{
  TIMER("gf_clover_leaf_field_5");
  GaugeField gf1;
  gf1.init(geo_resize(gf.geo(), 3));
  gf1 = gf;
  refresh_expanded(gf1);
  gf_clover_leaf_field_m_n_no_comm(clf1, gf1, 1, 1);
  gf_clover_leaf_field_m_n_no_comm(clf2, gf1, 2, 2);
  gf_clover_leaf_field_m_n_no_comm(clf3, gf1, 1, 2);
  gf_clover_leaf_field_m_n_no_comm(clf4, gf1, 1, 3);
  gf_clover_leaf_field_m_n_no_comm(clf5, gf1, 3, 3);
}

qacc double clf_plaq_action_density(const CloverLeafField& clf,
                                    const Coordinate& xl)
// \sum_P (1 - 1/3 * Re Tr U_P)
//
// Action = beta * total_volume() * action_density
// Single instanton action = 8 * sqr(PI) / g^2
// beta = 6/g^2
{
  const Vector<ColorMatrix> v = clf.get_elems_const(xl);
  double sum = 0.0;
  for (int i = 0; i < 6; ++i) {
    sum += 1.0 - 1.0 / 3.0 * matrix_trace(v[i]).real();
  }
  return sum;
}

qacc double clf_spatial_plaq_action_density(const CloverLeafField& clf,
                                            const Coordinate& xl)
// \sum_P(spatial only) (1 - 1/3 * Re Tr U_P)
{
  const Vector<ColorMatrix> v = clf.get_elems_const(xl);
  double sum = 0.0;
  sum += 1.0 - 1.0 / 3.0 * matrix_trace(v[0]).real();
  sum += 1.0 - 1.0 / 3.0 * matrix_trace(v[1]).real();
  sum += 1.0 - 1.0 / 3.0 * matrix_trace(v[3]).real();
  return sum;
}

qacc double clf_topology_density(const CloverLeafField& clf,
                                 const Coordinate& xl)
// sum of the density of the topological charge Q
{
  const Vector<ColorMatrix> v = clf.get_elems_const(xl);
  array<ColorMatrix, 6> arr;
  for (int i = 0; i < 6; ++i) {
    arr[i] = (Complex)0.5 * (v[i] - matrix_adjoint(v[i]));
  }
  const double fac = -1.0 / (4.0 * PI * PI);
  double sum = 0.0;
  sum -= matrix_trace(arr[1] * arr[4]).real();
  sum += matrix_trace(arr[2] * arr[3]).real();
  sum += matrix_trace(arr[5] * arr[0]).real();
  return fac * sum;
}

inline void clf_plaq_action_field(FieldM<double, 1>& paf,
                                  const CloverLeafField& clf)
{
  TIMER("clf_plaq_action_field");
  const Geometry& geo = clf.geo();
  paf.init(geo);
  qassert(is_matching_geo(paf.geo(), geo));
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    paf.get_elem(xl) = clf_plaq_action_density(clf, xl);
  });
}

inline void clf_spatial_plaq_action_field(FieldM<double, 1>& spaf,
                                          const CloverLeafField& clf)
{
  TIMER("clf_spatial_plaq_action_field");
  const Geometry& geo = clf.geo();
  spaf.init(geo);
  qassert(is_matching_geo(spaf.geo(), geo));
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    spaf.get_elem(xl) = clf_spatial_plaq_action_density(clf, xl);
  });
}

inline void clf_topology_field(FieldM<double, 1>& topf,
                               const CloverLeafField& clf)
{
  TIMER("clf_topology_field");
  const Geometry& geo = clf.geo();
  topf.init(geo);
  qassert(is_matching_geo(topf.geo(), geo));
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    topf.get_elem(xl) = clf_topology_density(clf, xl);
  });
}

inline void clf_topology_field(FieldM<double, 1>& topf, const GaugeField& gf)
{
  TIMER("clf_topology_field");
  CloverLeafField clf;
  gf_clover_leaf_field(clf, gf);
  clf_topology_field(topf, clf);
}

inline void clf_topology_field_5(FieldM<double, 1>& topf,
                                 const CloverLeafField& clf1,
                                 const CloverLeafField& clf2,
                                 const CloverLeafField& clf3,
                                 const CloverLeafField& clf4,
                                 const CloverLeafField& clf5)
{
  TIMER("clf_topology_field_5");
  const Geometry& geo = clf1.geo();
  topf.init(geo);
  qassert(is_matching_geo(topf.geo(), geo));
  qassert(is_matching_geo(geo, clf2.geo()));
  qassert(is_matching_geo(geo, clf3.geo()));
  qassert(is_matching_geo(geo, clf4.geo()));
  qassert(is_matching_geo(geo, clf5.geo()));
  const double c5 = 1.0 / 20.0;
  const double c1 = (19.0 - 55.0 * c5) / 9.0;
  const double c2 = (1.0 - 64.0 * c5) / 9.0;
  const double c3 = (-64.0 + 640.0 * c5) / 45.0;
  const double c4 = 1.0 / 5.0 - 2.0 * c5;
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    topf.get_elem(xl) = c1 * clf_topology_density(clf1, xl);
    topf.get_elem(xl) += c2 / 16.0 * clf_topology_density(clf2, xl);
    topf.get_elem(xl) += c3 / 4.0 * clf_topology_density(clf3, xl);
    topf.get_elem(xl) += c4 / 9.0 * clf_topology_density(clf4, xl);
    topf.get_elem(xl) += c5 / 81.0 * clf_topology_density(clf5, xl);
  });
}

inline void clf_topology_field_5(FieldM<double, 1>& topf, const GaugeField& gf)
// https://arxiv.org/pdf/hep-lat/9701012v2.pdf
{
  TIMER("clf_topology_field_5(topf,gf)");
  CloverLeafField clf1, clf2, clf3, clf4, clf5;
  gf_clover_leaf_field_5(clf1, clf2, clf3, clf4, clf5, gf);
  clf_topology_field_5(topf, clf1, clf2, clf3, clf4, clf5);
}

inline void clf_topology_field_5_terms(FieldM<double, 5>& topf,
                                       const CloverLeafField& clf1,
                                       const CloverLeafField& clf2,
                                       const CloverLeafField& clf3,
                                       const CloverLeafField& clf4,
                                       const CloverLeafField& clf5)
{
  TIMER("clf_topology_field_5");
  const Geometry& geo = clf1.geo();
  topf.init(geo);
  qassert(is_matching_geo(topf.geo(), geo));
  qassert(is_matching_geo(geo, clf2.geo()));
  qassert(is_matching_geo(geo, clf3.geo()));
  qassert(is_matching_geo(geo, clf4.geo()));
  qassert(is_matching_geo(geo, clf5.geo()));
  const double c5 = 1.0 / 20.0;
  const double c1 = (19.0 - 55.0 * c5) / 9.0;
  const double c2 = (1.0 - 64.0 * c5) / 9.0;
  const double c3 = (-64.0 + 640.0 * c5) / 45.0;
  const double c4 = 1.0 / 5.0 - 2.0 * c5;
  qacc_for(index, geo.local_volume(), {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<double> v = topf.get_elems(xl);
    v[0] = c1 * clf_topology_density(clf1, xl);
    v[1] = c2 / 16.0 * clf_topology_density(clf2, xl);
    v[2] = c3 / 4.0 * clf_topology_density(clf3, xl);
    v[3] = c4 / 9.0 * clf_topology_density(clf4, xl);
    v[4] = c5 / 81.0 * clf_topology_density(clf5, xl);
  });
}

inline void clf_topology_field_5_terms(FieldM<double, 5>& topf, const GaugeField& gf)
// https://arxiv.org/pdf/hep-lat/9701012v2.pdf
{
  TIMER("clf_topology_field_5(topf,gf)");
  CloverLeafField clf1, clf2, clf3, clf4, clf5;
  gf_clover_leaf_field_5(clf1, clf2, clf3, clf4, clf5, gf);
  clf_topology_field_5_terms(topf, clf1, clf2, clf3, clf4, clf5);
}

inline double topology_charge_5(const GaugeField& gf)
{
  TIMER("topology_charge_5(gf)");
  FieldM<double, 1> topf;
  clf_topology_field_5(topf, gf);
  const Geometry& geo = topf.geo();
  qassert(geo.is_only_local);
  double sum = 0.0;
  for (long index = 0; index < geo.local_volume(); ++index) {
    sum += topf.get_elem(index);
  }
  return sum;
}

}  // namespace qlat
