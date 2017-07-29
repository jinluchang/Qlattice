#pragma once

#include <qlat/matrix.h>
#include <qlat/qcd.h>
#include <qlat/qcd-utils.h>

QLAT_START_NAMESPACE

struct CloverLeafField : FieldM<ColorMatrix,6>
{
  virtual const std::string& cname()
  {
    static const std::string s = "CloverLeafField";
    return s;
  }
};

inline ColorMatrix gf_clover_leaf_no_comm(const GaugeField& gf1, const Coordinate& xl, const int mu, const int nu)
{
  ColorMatrix m;
  set_zero(m);
  std::vector<int> path(4);
  path[0] = mu;
  path[1] = nu;
  path[2] = -mu-1;
  path[3] = -nu-1;
  m += gf_wilson_line_no_comm(gf1, xl, path);
  path[0] = -nu-1;
  path[1] = -mu-1;
  path[2] = nu;
  path[3] = mu;
  m += matrix_adjoint(gf_wilson_line_no_comm(gf1, xl, path));
  path[0] = nu;
  path[1] = -mu-1;
  path[2] = -nu-1;
  path[3] = mu;
  m += gf_wilson_line_no_comm(gf1, xl, path);
  path[0] = mu;
  path[1] = -nu-1;
  path[2] = -mu-1;
  path[3] = nu;
  m += matrix_adjoint(gf_wilson_line_no_comm(gf1, xl, path));
  return 0.25 * m;
}

inline void gf_clover_leaf_field_no_comm(CloverLeafField& clf, const GaugeField& gf1)
  // F_01, F_02, F_03, F_12, F_13, F_23
{
  TIMER("gf_clover_leaf_field_no_comm");
  const Geometry geo = geo_reform(gf1.geo, 6, 0);
  clf.init(geo);
  qassert(is_matching_geo_mult(clf.geo, geo));
#pragma omp parallel
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = clf.get_elems(xl);
    v[0] = gf_clover_leaf_no_comm(gf1, xl, 0, 1);
    v[1] = gf_clover_leaf_no_comm(gf1, xl, 0, 2);
    v[2] = gf_clover_leaf_no_comm(gf1, xl, 0, 3);
    v[3] = gf_clover_leaf_no_comm(gf1, xl, 1, 2);
    v[4] = gf_clover_leaf_no_comm(gf1, xl, 1, 3);
    v[5] = gf_clover_leaf_no_comm(gf1, xl, 2, 3);
  }
}

inline void gf_clover_leaf_field(CloverLeafField& clf, const GaugeField& gf)
{
  TIMER("gf_clover_leaf_field");
  GaugeField gf1;
  gf1.init(geo_resize(gf.geo, 1));
  gf1 = gf;
  refresh_expanded(gf1);
  gf_clover_leaf_field_no_comm(clf, gf1);
}

inline double clf_plaq_action_density(const CloverLeafField& clf, const Coordinate& xl)
  // \sum_P (1 - 1/3 * Re Tr U_P)
{
  const Vector<ColorMatrix> v = clf.get_elems_const(xl);
  double sum = 0.0;
  for (int i = 0; i < 6; ++i) {
    sum += 1.0 - 1.0/3.0 * matrix_trace(v[i]).real();
  }
  return sum;
}

inline double clf_topology_density(const CloverLeafField& clf, const Coordinate& xl)
  // sum of the density of the topological charge Q
{
  const Vector<ColorMatrix> v = clf.get_elems_const(xl);
  std::array<ColorMatrix,6> arr;
  for (int i = 0; i < 6; ++i) {
    arr[i] = 0.5 * (v[i] - matrix_adjoint(v[i]));
  }
  const double fac = -1.0 / (4.0 * PI * PI);
  double sum = 0.0;
  sum -= matrix_trace(arr[1] * arr[4]).real();
  sum += matrix_trace(arr[2] * arr[3]).real();
  sum += matrix_trace(arr[5] * arr[0]).real();
  return fac * sum;
}

inline void clf_plaq_action_field(FieldM<double,1>& paf, const CloverLeafField& clf)
{
  TIMER("clf_plaq_action_field");
  const Geometry& geo = clf.geo;
  paf.init(geo);
  qassert(is_matching_geo(paf.geo, geo));
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate& xl = geo.coordinate_from_index(index);
    paf.get_elem(xl) = clf_plaq_action_density(clf, xl);
  }
}

inline void clf_topology_field(FieldM<double,1>& topf, const CloverLeafField& clf)
{
  TIMER("clf_topology_field");
  const Geometry& geo = clf.geo;
  topf.init(geo);
  qassert(is_matching_geo(topf.geo, geo));
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate& xl = geo.coordinate_from_index(index);
    topf.get_elem(xl) = clf_topology_density(clf, xl);
  }
}

QLAT_END_NAMESPACE
