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
                                      const GaugeField& gf1, const Int m,
                                      const Int n);

void gf_clover_leaf_field(CloverLeafField& clf, const GaugeField& gf);

void clf_plaq_action_density_field(Field<RealD>& paf,
                                   const CloverLeafField& clf);

void clf_spatial_plaq_action_field(Field<RealD>& spaf,
                                   const CloverLeafField& clf);

void clf_topology_field(Field<RealD>& topf, const CloverLeafField& clf);

}  // namespace qlat
