#include <qlat/qcd.h>
#include <qlat/qcd-topology.h>

namespace qlat
{  //

RealD gf_avg_spatial_plaq(const GaugeField& gf)
{
  return gf_avg_spatial_plaq<RealD>(gf);
}

RealD gf_avg_plaq(const GaugeField& gf) { return gf_avg_plaq<RealD>(gf); }

RealD gf_avg_link_trace(const GaugeField& gf)
{
  return gf_avg_link_trace<RealD>(gf);
}

RealD topology_charge_5(const GaugeField& gf)
// interface function
{
  TIMER("topology_charge_5(gf)");
  FieldM<RealD, 1> topf;
  clf_topology_field_5(topf, gf);
  const Geometry& geo = topf.geo();
  qassert(geo.is_only_local);
  RealD sum = 0.0;
  for (Long index = 0; index < geo.local_volume(); ++index) {
    sum += topf.get_elem(index);
  }
  return sum;
}

void clf_plaq_action_density_field(Field<RealD>& paf, const GaugeField& gf)
// interface function
// \sum_P (1 - 1/3 * Re Tr U_P)
//
// Action = beta * total_volume() * action_density
// Single instanton action = 8 * sqr(PI) / g^2
// beta = 6/g^2
{
  TIMER("clf_plaq_action_density_field");
  CloverLeafField clf;
  gf_clover_leaf_field(clf, gf);
  clf_plaq_action_density_field(paf, clf);
}

void clf_spatial_plaq_action_density_field(Field<RealD>& paf, const GaugeField& gf)
// interface function
// \sum_P(spatial only) (1 - 1/3 * Re Tr U_P)
{
  TIMER("clf_spatial_plaq_action_density_field");
  CloverLeafField clf;
  gf_clover_leaf_field(clf, gf);
  clf_spatial_plaq_action_field(paf, clf);
}

void clf_topology_field(Field<RealD>& topf, const GaugeField& gf)
// interface function
{
  TIMER("clf_topology_field");
  CloverLeafField clf;
  gf_clover_leaf_field(clf, gf);
  clf_topology_field(topf, clf);
}

void clf_topology_field_5_terms(Field<RealD>& topf, const GaugeField& gf)
// interface function
// https://arxiv.org/pdf/hep-lat/9701012v2.pdf
// topf.geo.multiplicity == 5
{
  TIMER("clf_topology_field_5_terms(topf,gf)");
  CloverLeafField clf1, clf2, clf3, clf4, clf5;
  gf_clover_leaf_field_5(clf1, clf2, clf3, clf4, clf5, gf);
  clf_topology_field_5_terms(topf, clf1, clf2, clf3, clf4, clf5);
}

void clf_topology_field_5(Field<RealD>& topf, const GaugeField& gf)
// interface function
// https://arxiv.org/pdf/hep-lat/9701012v2.pdf
{
  TIMER("clf_topology_field_5(topf,gf)");
  CloverLeafField clf1, clf2, clf3, clf4, clf5;
  gf_clover_leaf_field_5(clf1, clf2, clf3, clf4, clf5, gf);
  clf_topology_field_5(topf, clf1, clf2, clf3, clf4, clf5);
}

}  // namespace qlat
