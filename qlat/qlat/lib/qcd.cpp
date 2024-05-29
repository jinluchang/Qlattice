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

// ------------------------------------

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

// ------------------------------------

void gf_clover_leaf_field_no_comm(CloverLeafField& clf, const GaugeField& gf1)
// F_01, F_02, F_03, F_12, F_13, F_23
{
  TIMER("gf_clover_leaf_field_no_comm");
  const Geometry geo = geo_resize(gf1.geo());
  const Int multiplicity = 6;
  clf.init(geo, multiplicity);
  qassert(is_matching_geo(clf.geo(), geo));
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = clf.geo();
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

void gf_clover_leaf_field_m_n_no_comm(CloverLeafField& clf,
                                      const GaugeField& gf1, const int m,
                                      const int n)
// F_01, F_02, F_03, F_12, F_13, F_23
{
  TIMER("gf_clover_leaf_field_m_n_no_comm");
  const Geometry geo = geo_resize(gf1.geo());
  const Int multiplicity = 6;
  clf.init(geo, multiplicity);
  qassert(is_matching_geo(clf.geo(), geo));
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = clf.geo();
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

void gf_clover_leaf_field(CloverLeafField& clf, const GaugeField& gf)
{
  TIMER("gf_clover_leaf_field");
  GaugeField gf1;
  gf1.init(geo_resize(gf.geo(), 1));
  gf1 = gf;
  refresh_expanded(gf1);
  gf_clover_leaf_field_no_comm(clf, gf1);
}

void gf_clover_leaf_field_5(CloverLeafField& clf1, CloverLeafField& clf2,
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

void clf_plaq_action_density_field(Field<RealD>& paf,
                                   const CloverLeafField& clf)
{
  TIMER("clf_plaq_action_density_field");
  const Geometry& geo = clf.geo();
  paf.init(geo, 1);
  qassert(is_matching_geo(paf.geo(), geo));
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = paf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    paf.get_elem(xl) = clf_plaq_action_density(clf, xl);
  });
}

void clf_spatial_plaq_action_field(Field<RealD>& spaf,
                                   const CloverLeafField& clf)
{
  TIMER("clf_spatial_plaq_action_field");
  const Geometry& geo = clf.geo();
  spaf.init(geo, 1);
  qassert(is_matching_geo(spaf.geo(), geo));
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = spaf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    spaf.get_elem(xl) = clf_spatial_plaq_action_density(clf, xl);
  });
}

void clf_topology_field(Field<RealD>& topf, const CloverLeafField& clf)
{
  TIMER("clf_topology_field");
  const Geometry& geo = clf.geo();
  topf.init(geo, 1);
  qassert(is_matching_geo(topf.geo(), geo));
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = clf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    topf.get_elem(xl) = clf_topology_density(clf, xl);
  });
}

void clf_topology_field_5(Field<RealD>& topf, const CloverLeafField& clf1,
                          const CloverLeafField& clf2,
                          const CloverLeafField& clf3,
                          const CloverLeafField& clf4,
                          const CloverLeafField& clf5)
{
  TIMER("clf_topology_field_5");
  const Geometry& geo = clf1.geo();
  topf.init(geo, 1);
  qassert(is_matching_geo(geo, topf.geo()));
  qassert(is_matching_geo(geo, clf1.geo()));
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

void clf_topology_field_5_terms(Field<RealD>& topf, const CloverLeafField& clf1,
                                const CloverLeafField& clf2,
                                const CloverLeafField& clf3,
                                const CloverLeafField& clf4,
                                const CloverLeafField& clf5)
// topf.geo.multiplicity == 5
{
  TIMER("clf_topology_field_5_terms");
  const Geometry& geo = clf1.geo();
  topf.init(geo, 5);
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
    const Geometry& geo = clf1.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<double> v = topf.get_elems(xl);
    v[0] = c1 * clf_topology_density(clf1, xl);
    v[1] = c2 / 16.0 * clf_topology_density(clf2, xl);
    v[2] = c3 / 4.0 * clf_topology_density(clf3, xl);
    v[3] = c4 / 9.0 * clf_topology_density(clf4, xl);
    v[4] = c5 / 81.0 * clf_topology_density(clf5, xl);
  });
}

}  // namespace qlat
