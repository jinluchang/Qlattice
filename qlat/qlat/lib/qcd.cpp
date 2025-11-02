#include <qlat/qcd-acc.h>
#include <qlat/qcd-topology.h>
#include <qlat/qcd.h>

namespace qlat
{  //

static qacc ColorMatrix gf_clover_leaf_m_n_no_comm(const GaugeField& gf1,
                                                   const Coordinate& xl,
                                                   const Int mu, const Int m,
                                                   const Int nu, const Int n)
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

static qacc RealD clf_plaq_action_density(const CloverLeafField& clf,
                                          const Coordinate& xl)
// \sum_P (1 - 1/3 * Re Tr U_P)
//
// 6 plaqs in the above sum (for each site).
//
// Action = beta * total_volume * 6 * (1 - avg_plaq)
// Action = beta * total_volume() * action_density
// Single instanton action = 8 * sqr(PI) / g^2
// beta = 6/g^2
{
  const Vector<ColorMatrix> v = clf.get_elems_const(xl);
  RealD sum = 0.0;
  for (Int i = 0; i < 6; ++i) {
    sum += 1.0 - 1.0 / 3.0 * matrix_trace(v[i]).real();
  }
  return sum;
}

static qacc RealD clf_spatial_plaq_action_density(const CloverLeafField& clf,
                                                  const Coordinate& xl)
// \sum_P(spatial only) (1 - 1/3 * Re Tr U_P)
{
  const Vector<ColorMatrix> v = clf.get_elems_const(xl);
  RealD sum = 0.0;
  sum += 1.0 - 1.0 / 3.0 * matrix_trace(v[0]).real();
  sum += 1.0 - 1.0 / 3.0 * matrix_trace(v[1]).real();
  sum += 1.0 - 1.0 / 3.0 * matrix_trace(v[3]).real();
  return sum;
}

static qacc RealD clf_topology_density(const CloverLeafField& clf,
                                       const Coordinate& xl)
// sum of the density of the topological charge Q
{
  const Vector<ColorMatrix> v = clf.get_elems_const(xl);
  array<ColorMatrix, 6> arr;
  for (Int i = 0; i < 6; ++i) {
    arr[i] = (ComplexD)0.5 * (v[i] - matrix_adjoint(v[i]));
  }
  const RealD fac = -1.0 / (4.0 * PI * PI);
  RealD sum = 0.0;
  sum -= matrix_trace(arr[1] * arr[4]).real();
  sum += matrix_trace(arr[2] * arr[3]).real();
  sum += matrix_trace(arr[5] * arr[0]).real();
  return fac * sum;
}

// ------------------------------------

template <class T>
qacc RealD gf_plaq_no_comm(const Vector<ColorMatrixT<T>>& v,
                           const array<Vector<ColorMatrixT<T>>, DIMN>& vms,
                           const Int m1, const Int m2)
{
  ColorMatrixT<T> cm =
      v[m1] * vms[m1][m2] * matrix_adjoint(v[m2] * vms[m2][m1]);
  return matrix_trace(cm).real() / NUM_COLOR;
}

template <class T>
static RealD gf_avg_plaq_no_comm(const GaugeFieldT<T>& gf)
// assume proper communication is done
{
  TIMER("gf_avg_plaq_no_comm");
  const Geometry geo = geo_resize(gf.geo());
  FieldM<RealD, 1> cf;
  cf.init(geo);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = cf.geo();
    Coordinate xl = geo.coordinate_from_index(index);
    const Vector<ColorMatrixT<T>> v = gf.get_elems_const(xl);
    array<Vector<ColorMatrixT<T>>, DIMN> vms;
    for (Int m = 0; m < DIMN; ++m) {
      xl[m] += 1;
      vms[m] = gf.get_elems_const(xl);
      xl[m] -= 1;
    }
    RealD sum = 0.0;
    sum += gf_plaq_no_comm(v, vms, 0, 1);
    sum += gf_plaq_no_comm(v, vms, 0, 2);
    sum += gf_plaq_no_comm(v, vms, 0, 3);
    sum += gf_plaq_no_comm(v, vms, 1, 2);
    sum += gf_plaq_no_comm(v, vms, 1, 3);
    sum += gf_plaq_no_comm(v, vms, 2, 3);
    sum /= DIMN * (DIMN - 1) / 2;
    if (std::isnan(sum)) {
      qerr(ssprintf("WARNING: isnan in gf_avg_plaq"));
    }
    cf.get_elem(index) = sum;
  });
  const std::vector<RealD> sum_vec = field_sum(cf);
  qassert(sum_vec.size() == 1);
  RealD sum = sum_vec[0];
  glb_sum(sum);
  sum /= geo.total_volume();
  return sum;
}

template <class T>
static RealD gf_avg_plaq(const GaugeFieldT<T>& gf)
{
  TIMER("gf_avg_plaq");
  GaugeFieldT<T> gf1;
  gf1.init(geo_resize(gf.geo(), Coordinate(0, 0, 0, 0), Coordinate(1, 1, 1, 1)));
  gf1 = gf;
  refresh_expanded(gf1);
  return gf_avg_plaq_no_comm(gf1);
}

template <class T>
static RealD gf_avg_spatial_plaq_no_comm(const GaugeFieldT<T>& gf)
// assume proper communication is done
{
  TIMER("gf_avg_spatial_plaq_no_comm");
  const Geometry& geo = gf.geo();
  std::vector<RealD> sums(omp_get_max_threads(), 0.0);
#pragma omp parallel
  {
    RealD sum_avg_plaq = 0.0;
#pragma omp for
    for (Long index = 0; index < geo.local_volume(); ++index) {
      Coordinate xl = geo.coordinate_from_index(index);
      const Vector<ColorMatrixT<T> > v = gf.get_elems_const(xl);
      std::vector<Vector<ColorMatrixT<T> > > vms(DIMN - 1);
      for (Int m = 0; m < DIMN - 1; ++m) {
        xl[m] += 1;
        vms[m] = gf.get_elems_const(xl);
        xl[m] -= 1;
      }
      RealD avg_plaq = 0.0;
      for (Int m1 = 1; m1 < 3; ++m1) {
        for (Int m2 = 0; m2 < m1; ++m2) {
          ColorMatrixT<T> cm =
              v[m1] * vms[m1][m2] * matrix_adjoint(v[m2] * vms[m2][m1]);
          avg_plaq += matrix_trace(cm).real() / NUM_COLOR;
          if (std::isnan(avg_plaq)) {
            fdisplayln(stdout, ssprintf("WARNING: isnan in gf_avg_plaq"));
            qassert(false);
          }
        }
      }
      avg_plaq /= (DIMN - 1) * (DIMN - 2) / 2;
      sum_avg_plaq += avg_plaq;
    }
    sums[omp_get_thread_num()] = sum_avg_plaq;
  }
  RealD sum = 0.0;
  for (size_t i = 0; i < sums.size(); ++i) {
    sum += sums[i];
  }
  glb_sum(sum);
  sum /= geo.total_volume();
  return sum;
}

template <class T>
static RealD gf_avg_spatial_plaq(const GaugeFieldT<T>& gf)
{
  TIMER("gf_avg_spatial_plaq(gf)");
  GaugeFieldT<T> gf1;
  gf1.init(geo_resize(gf.geo(), Coordinate(0, 0, 0, 0), Coordinate(1, 1, 1, 0)));
  gf1 = gf;
  refresh_expanded(gf1);
  return gf_avg_spatial_plaq_no_comm(gf1);
}

template <class T>
static RealD gf_avg_link_trace(const GaugeFieldT<T>& gf)
{
  TIMER("gf_avg_link_trace");
  const Geometry& geo = gf.geo();
  FieldM<RealD, 1> cf;
  cf.init(geo);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = cf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    const Vector<ColorMatrixT<T>> v = gf.get_elems_const(xl);
    RealD sum = 0;
    for (Int m = 0; m < v.size(); ++m) {
      sum += matrix_trace(v[m]).real() / NUM_COLOR;
    }
    sum /= v.size();
    cf.get_elem(index) = sum;
  });
  const std::vector<RealD> sum_vec = field_sum(cf);
  qassert(sum_vec.size() == 1);
  RealD sum = sum_vec[0];
  glb_sum(sum);
  sum /= geo.total_volume();
  return sum;
}

// ------------------------------------

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

template <class T>
static void gf_plaq_field_no_comm(Field<RealD>& f_plaq,
                                  const GaugeFieldT<T>& gf)
// assume proper communication is done
{
  TIMER("gf_plaq_field_no_comm");
  const Geometry geo = geo_resize(gf.get_geo());
  f_plaq.init(geo, 6);
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = f_plaq.geo();
    Coordinate xl = geo.coordinate_from_index(index);
    const Vector<ColorMatrixT<T>> v = gf.get_elems_const(xl);
    array<Vector<ColorMatrixT<T>>, DIMN> vms;
    for (Int m = 0; m < DIMN; ++m) {
      xl[m] += 1;
      vms[m] = gf.get_elems_const(xl);
      xl[m] -= 1;
    }
    Vector<RealD> pv = f_plaq.get_elems(index);
    pv[0] = gf_plaq_no_comm(v, vms, 0, 1);
    pv[1] = gf_plaq_no_comm(v, vms, 0, 2);
    pv[2] = gf_plaq_no_comm(v, vms, 0, 3);
    pv[3] = gf_plaq_no_comm(v, vms, 1, 2);
    pv[4] = gf_plaq_no_comm(v, vms, 1, 3);
    pv[5] = gf_plaq_no_comm(v, vms, 2, 3);
  });
}

template <class T>
static void gf_plaq_field(Field<RealD>& f_plaq, const GaugeFieldT<T>& gf)
{
  TIMER("gf_plaq_field");
  GaugeFieldT<T> gf1;
  gf1.init(
      geo_resize(gf.geo(), Coordinate(0, 0, 0, 0), Coordinate(1, 1, 1, 1)));
  gf1 = gf;
  refresh_expanded(gf1);
  gf_plaq_field_no_comm(f_plaq, gf1);
}

void gf_plaq_field(Field<RealD>& f_plaq, const GaugeField& gf)
{
  gf_plaq_field<Real>(f_plaq, gf);
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
// topf.multiplicity == 5
{
  TIMER("clf_topology_field_5_terms(topf,gf)");
  const Geometry& geo = gf.geo();
  topf.init(geo, 5);
  qassert(is_matching_geo(topf.geo(), geo));
  const RealD c5 = 1.0 / 20.0;
  const RealD c1 = (19.0 - 55.0 * c5) / 9.0;
  const RealD c2 = (1.0 - 64.0 * c5) / 9.0;
  const RealD c3 = (-64.0 + 640.0 * c5) / 45.0;
  const RealD c4 = 1.0 / 5.0 - 2.0 * c5;
  GaugeField gf1;
  gf1.init(geo_resize(gf.geo(), 3));
  gf1 = gf;
  refresh_expanded(gf1);
  CloverLeafField clf1, clf2, clf3, clf4, clf5;
  gf_clover_leaf_field_m_n_no_comm(clf1, gf1, 1, 1);
  qassert(is_matching_geo(geo, clf1.geo()));
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = gf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<RealD> v = topf.get_elems(xl);
    v[0] = c1 * clf_topology_density(clf1, xl);
  });
  clf1.init();
  gf_clover_leaf_field_m_n_no_comm(clf2, gf1, 2, 2);
  qassert(is_matching_geo(geo, clf2.geo()));
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = gf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<RealD> v = topf.get_elems(xl);
    v[1] = c2 / 16.0 * clf_topology_density(clf2, xl);
  });
  clf2.init();
  gf_clover_leaf_field_m_n_no_comm(clf3, gf1, 1, 2);
  qassert(is_matching_geo(geo, clf3.geo()));
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = gf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<RealD> v = topf.get_elems(xl);
    v[2] = c3 / 4.0 * clf_topology_density(clf3, xl);
  });
  clf3.init();
  gf_clover_leaf_field_m_n_no_comm(clf4, gf1, 1, 3);
  qassert(is_matching_geo(geo, clf4.geo()));
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = gf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<RealD> v = topf.get_elems(xl);
    v[3] = c4 / 9.0 * clf_topology_density(clf4, xl);
  });
  clf4.init();
  gf_clover_leaf_field_m_n_no_comm(clf5, gf1, 3, 3);
  qassert(is_matching_geo(geo, clf5.geo()));
  qacc_for(index, geo.local_volume(), {
    const Geometry& geo = gf.geo();
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<RealD> v = topf.get_elems(xl);
    v[4] = c5 / 81.0 * clf_topology_density(clf5, xl);
  });
  clf5.init();
}

void clf_topology_field_5(Field<RealD>& topf, const GaugeField& gf)
// interface function
// https://arxiv.org/pdf/hep-lat/9701012v2.pdf
// topf.multiplicity == 1
{
  TIMER("clf_topology_field_5(topf,gf)");
  const Geometry& geo = gf.geo();
  topf.init(geo, 1);
  qassert(is_matching_geo(geo, topf.geo()));
  Field<RealD> top_terms_f;
  clf_topology_field_5_terms(top_terms_f, gf);
  qassert(is_matching_geo(geo, top_terms_f.geo()));
  qacc_for(index, geo.local_volume(), {
    RealD& s = topf.get_elem(index);
    s = 0;
    Vector<RealD> v = top_terms_f.get_elems(index);
    for (Int m = 0; m < v.size(); ++m) {
      s += v[m];
    }
  });
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
                                      const GaugeField& gf1, const Int m,
                                      const Int n)
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
      for (Int i = 0; i < 6; ++i) {
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

}  // namespace qlat
