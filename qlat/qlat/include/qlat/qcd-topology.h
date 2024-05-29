#pragma once

#include <qlat-utils/matrix.h>
#include <qlat/qcd-utils.h>
#include <qlat/qcd.h>

namespace qlat
{  //

RealD topology_charge_5(const GaugeField& gf);

void clf_plaq_action_density_field(Field<RealD>& paf, const GaugeField& gf);

void clf_spatial_plaq_action_density_field(Field<RealD>& paf,
                                           const GaugeField& gf);

void clf_topology_field(Field<RealD>& topf, const GaugeField& gf);

void clf_topology_field_5(Field<RealD>& topf, const GaugeField& gf);

void clf_topology_field_5_terms(Field<RealD>& topf, const GaugeField& gf);

// ----------------------------------------------

struct CloverLeafField : FieldM<ColorMatrix, 6> {
};

void gf_clover_leaf_field_no_comm(CloverLeafField& clf, const GaugeField& gf1);

void gf_clover_leaf_field_m_n_no_comm(CloverLeafField& clf,
                                      const GaugeField& gf1, const int m,
                                      const int n);

void gf_clover_leaf_field(CloverLeafField& clf, const GaugeField& gf);

void clf_plaq_action_density_field(Field<RealD>& paf,
                                   const CloverLeafField& clf);

void clf_spatial_plaq_action_field(Field<RealD>& spaf,
                                   const CloverLeafField& clf);

void clf_topology_field(Field<RealD>& topf, const CloverLeafField& clf);

// ----------------------------------------------

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
  return (ComplexD)0.25 * m;
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
  return (ComplexD)0.25 * cm;
}

qacc RealD clf_plaq_action_density(const CloverLeafField& clf,
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

qacc RealD clf_spatial_plaq_action_density(const CloverLeafField& clf,
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

qacc RealD clf_topology_density(const CloverLeafField& clf,
                                const Coordinate& xl)
// sum of the density of the topological charge Q
{
  const Vector<ColorMatrix> v = clf.get_elems_const(xl);
  array<ColorMatrix, 6> arr;
  for (int i = 0; i < 6; ++i) {
    arr[i] = (ComplexD)0.5 * (v[i] - matrix_adjoint(v[i]));
  }
  const double fac = -1.0 / (4.0 * PI * PI);
  double sum = 0.0;
  sum -= matrix_trace(arr[1] * arr[4]).real();
  sum += matrix_trace(arr[2] * arr[3]).real();
  sum += matrix_trace(arr[5] * arr[0]).real();
  return fac * sum;
}

}  // namespace qlat
