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

}  // namespace qlat
