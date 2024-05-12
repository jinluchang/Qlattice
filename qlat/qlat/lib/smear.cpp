#define QLAT_INSTANTIATE_SMEAR

#include <qlat/qcd-smear.h>
#include <qlat/vector_utils/utils_smear_vecs.h>

namespace qlat
{  //

ColorMatrix gf_link_ape_smear_no_comm(const GaugeField& gf,
                                      const Coordinate& xl, const int mu,
                                      const double alpha)
{
  return color_matrix_su_projection(
      (ComplexD)(1.0 - alpha) * gf.get_elem(xl, mu) +
      (ComplexD)(alpha / 6.0) * gf_staple_no_comm(gf, xl, mu));
}

void gf_ape_smear_no_comm(GaugeField& gf, const GaugeField& gf0,
                          const double alpha)
{
  TIMER_VERBOSE("gf_ape_smear_no_comm");
  qassert(&gf != &gf0);
  const Geometry& geo = gf0.geo();
  gf.init(geo_resize(geo));
  qassert(is_matching_geo(geo, gf.geo()));
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = gf.get_elems(xl);
    for (int mu = 0; mu < DIMN; ++mu) {
      v[mu] = gf_link_ape_smear_no_comm(gf0, xl, mu, alpha);
    }
  }
}

void gf_ape_smear(GaugeField& gf, const GaugeField& gf0, const double alpha,
                  const Long steps)
{
  TIMER_VERBOSE("gf_ape_smear");
  gf = gf0;
  GaugeField gf1;
  gf1.init(geo_resize(gf0.geo(), 1));
  for (Long i = 0; i < steps; ++i) {
    gf1 = gf;
    refresh_expanded(gf1);
    gf_ape_smear_no_comm(gf, gf1, alpha);
  }
}

ColorMatrix gf_link_spatial_ape_smear_no_comm(const GaugeField& gf,
                                              const Coordinate& xl,
                                              const int mu, const double alpha)
{
  const double multi = mu == 3 ? 6.0 : 4.0;
  return color_matrix_su_projection(
      (ComplexD)(1.0 - alpha) * gf.get_elem(xl, mu) +
      (ComplexD)(alpha / multi) * gf_spatial_staple_no_comm(gf, xl, mu));
}

void gf_spatial_ape_smear_no_comm(GaugeField& gf, const GaugeField& gf0,
                                  const double alpha)
{
  TIMER_VERBOSE("gf_spatial_ape_smear_no_comm");
  qassert(&gf != &gf0);
  const Geometry& geo = gf0.geo();
  gf.init(geo_resize(geo));
  qassert(is_matching_geo(geo, gf.geo()));
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = gf.get_elems(xl);
    for (int mu = 0; mu < 3; ++mu) {
      // No need to smear the temperal link (mu == 3)
      v[mu] = gf_link_spatial_ape_smear_no_comm(gf0, xl, mu, alpha);
    }
  }
}

void gf_spatial_ape_smear(GaugeField& gf, const GaugeField& gf0,
                          const double alpha, const Long steps)
{
  TIMER_VERBOSE("gf_spatial_ape_smear");
  gf = gf0;
  const Coordinate expan_left(1, 1, 1, 0);
  const Coordinate expan_right(1, 1, 1, 0);
  GaugeField gf1;
  gf1.init(geo_resize(gf0.geo(), expan_left, expan_right));
  for (Long i = 0; i < steps; ++i) {
    gf1 = gf;
    refresh_expanded(gf1);
    gf_spatial_ape_smear_no_comm(gf, gf1, alpha);
  }
}

ColorMatrix gf_link_hyp_smear_3_no_comm(const GaugeField& gf,
                                        const Coordinate& xl, const int mu,
                                        const int nu, const int rho,
                                        const double alpha3)
{
  ColorMatrix ret;
  set_zero(ret);
  const Coordinate xl_mu = coordinate_shifts(xl, mu);
  for (int m = 0; m < DIMN; ++m) {
    if (mu != m && nu != m && rho != m) {
      ret += gf.get_elem(xl, m) * gf.get_elem(coordinate_shifts(xl, m), mu) *
             matrix_adjoint(gf.get_elem(xl_mu, m));
      ret += matrix_adjoint(gf.get_elem(coordinate_shifts(xl, -m - 1), m)) *
             gf.get_elem(coordinate_shifts(xl, -m - 1), mu) *
             gf.get_elem(coordinate_shifts(xl_mu, -m - 1), m);
    }
  }
  ret = (ComplexD)(1.0 - alpha3) * gf.get_elem(xl, mu) +
        (ComplexD)(alpha3 / 2.0) * ret;
  return color_matrix_su_projection(ret);
}

ColorMatrix gf_link_hyp_smear_2_no_comm(const GaugeField& gf,
                                        const Coordinate& xl, const int mu,
                                        const int nu, const double alpha2,
                                        const double alpha3)
{
  ColorMatrix ret;
  set_zero(ret);
  const Coordinate xl_mu = coordinate_shifts(xl, mu);
  for (int m = 0; m < DIMN; ++m) {
    if (mu != m && nu != m) {
      ret += gf_link_hyp_smear_3_no_comm(gf, xl, m, mu, nu, alpha3) *
             gf_link_hyp_smear_3_no_comm(gf, coordinate_shifts(xl, m), mu, m,
                                         nu, alpha3) *
             matrix_adjoint(
                 gf_link_hyp_smear_3_no_comm(gf, xl_mu, m, mu, nu, alpha3));
      ret += matrix_adjoint(gf_link_hyp_smear_3_no_comm(
                 gf, coordinate_shifts(xl, -m - 1), m, mu, nu, alpha3)) *
             gf_link_hyp_smear_3_no_comm(gf, coordinate_shifts(xl, -m - 1), mu,
                                         m, nu, alpha3) *
             gf_link_hyp_smear_3_no_comm(gf, coordinate_shifts(xl_mu, -m - 1),
                                         m, mu, nu, alpha3);
    }
  }
  ret = (ComplexD)(1.0 - alpha2) * gf.get_elem(xl, mu) +
        (ComplexD)(alpha2 / 4.0) * ret;
  return color_matrix_su_projection(ret);
}

ColorMatrix gf_link_hyp_smear_1_no_comm(const GaugeField& gf,
                                        const Coordinate& xl, const int mu,
                                        const double alpha1,
                                        const double alpha2,
                                        const double alpha3)
{
  ColorMatrix ret;
  set_zero(ret);
  const Coordinate xl_mu = coordinate_shifts(xl, mu);
  for (int m = 0; m < DIMN; ++m) {
    if (mu != m) {
      ret += gf_link_hyp_smear_2_no_comm(gf, xl, m, mu, alpha2, alpha3) *
             gf_link_hyp_smear_2_no_comm(gf, coordinate_shifts(xl, m), mu, m,
                                         alpha2, alpha3) *
             matrix_adjoint(
                 gf_link_hyp_smear_2_no_comm(gf, xl_mu, m, mu, alpha2, alpha3));
      ret += matrix_adjoint(gf_link_hyp_smear_2_no_comm(
                 gf, coordinate_shifts(xl, -m - 1), m, mu, alpha2, alpha3)) *
             gf_link_hyp_smear_2_no_comm(gf, coordinate_shifts(xl, -m - 1), mu,
                                         m, alpha2, alpha3) *
             gf_link_hyp_smear_2_no_comm(gf, coordinate_shifts(xl_mu, -m - 1),
                                         m, mu, alpha2, alpha3);
    }
  }
  ret = (ComplexD)(1.0 - alpha1) * gf.get_elem(xl, mu) +
        (ComplexD)(alpha1 / 6.0) * ret;
  return color_matrix_su_projection(ret);
}

ColorMatrix gf_link_hyp_smear_no_comm(const GaugeField& gf,
                                      const Coordinate& xl, const int mu,
                                      const double alpha1, const double alpha2,
                                      const double alpha3)
{
  return gf_link_hyp_smear_1_no_comm(gf, xl, mu, alpha1, alpha2, alpha3);
}

void gf_hyp_smear_no_comm(GaugeField& gf, const GaugeField& gf0,
                          const double alpha1, const double alpha2,
                          const double alpha3)
{
  TIMER_VERBOSE("gf_hyp_smear_no_comm");
  qassert(&gf != &gf0);
  const Geometry& geo = gf0.geo();
  gf.init(geo_resize(geo));
  qassert(is_matching_geo(geo, gf.geo()));
#pragma omp parallel for
  for (Long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = gf.get_elems(xl);
    for (int mu = 0; mu < DIMN; ++mu) {
      v[mu] = gf_link_hyp_smear_no_comm(gf0, xl, mu, alpha1, alpha2, alpha3);
    }
  }
}

void gf_hyp_smear(GaugeField& gf, const GaugeField& gf0, const double alpha1,
                  const double alpha2, const double alpha3)
// values in paper is 0.75 0.6 0.3
// 10.1103/PhysRevD.64.034504 Eq(4)
{
  TIMER_VERBOSE("gf_hyp_smear");
  GaugeField gf1;
  gf1.init(geo_resize(gf0.geo(), 1));
  gf1 = gf0;
  refresh_expanded(gf1);
  gf_hyp_smear_no_comm(gf, gf1, alpha1, alpha2, alpha3);
}

}  // namespace qlat
