#pragma once

#include <qlat/matrix.h>
#include <qlat/qcd-utils.h>
#include <qlat/qcd.h>

namespace qlat
{  //

inline ColorMatrix color_matrix_sub_invert(const ColorMatrix& x, const int ind)
// get su2 submatrix of x and return the su3 matrix that
// has the inverse of this matrix in the relevant row and column
{
  const int su2_index[][3] = {{0, 1, 2}, {0, 2, 1}, {1, 2, 0}};
  const int i1 = su2_index[ind][0];
  const int i2 = su2_index[ind][1];
  // const int zero_rc = su2_index[ind][2];
  // project onto SU(2)
  double p0 = x(i1, i1).real() + x(i2, i2).real();
  double p1 = x(i1, i2).imag() + x(i2, i1).imag();
  double p2 = x(i1, i2).real() - x(i2, i1).real();
  double p3 = x(i1, i1).imag() - x(i2, i2).imag();
  const double psqr = sqrt(p0 * p0 + p1 * p1 + p2 * p2 + p3 * p3);
  ColorMatrix y;
  set_unit(y);
  if (psqr == 0.0) {
    return y;
  } else {
    double ipsqr = 1.0 / psqr;
    p0 *= ipsqr;
    p1 *= ipsqr;
    p2 *= ipsqr;
    p3 *= ipsqr;
    // fill with inverse
    y(i1, i1) = Complex(p0, -p3);
    y(i2, i2) = Complex(p0, p3);
    y(i1, i2) = Complex(-p2, -p1);
    y(i2, i1) = Complex(p2, -p1);
    return y;
  }
}

inline ColorMatrix color_matrix_su_projection(const ColorMatrix& x,
                                              const double tolerance = 1.0e-8)
{
  // usually takes ~5 hits, so just exit
  // if hits the max, as something is
  // probably very wrong.
  const int max_iter = 10000;
  const ColorMatrix xdag = matrix_adjoint(x);
  ColorMatrix tmp = xdag;
  ColorMatrix y;
  set_unit(y);
  double old_tr = matrix_trace(xdag).real();
  for (int i = 0; i < max_iter; i++) {
    // loop over su2 subgroups
    double diff = 0.0;
    for (int j = 0; j < 3; j++) {
      const ColorMatrix inv = color_matrix_sub_invert(tmp, j);
      // y  .DotMEqual( inv, ycopy );
      y = inv * y;
      // tmp.DotMEqual( y, xdag );
      tmp = y * xdag;
      const double tr = matrix_trace(tmp).real();
      const double dtr = tr - old_tr;
      if (dtr > diff) {
        diff = dtr;
      }
      old_tr = tr;
    }
    // for single precision the difference seems
    // to never get below 1e-7 (not too suprising)
    if (diff < tolerance) {
      break;
    }
    qassert(i < max_iter - 1);
  }
  unitarize(y);
  return y;
}

inline ColorMatrix gf_link_ape_smear_no_comm(const GaugeField& gf,
                                             const Coordinate& xl, const int mu,
                                             const double alpha)
{
  return color_matrix_su_projection(
      (ComplexT)(1.0 - alpha) * gf.get_elem(xl, mu) +
      (ComplexT)(alpha / 6.0) * gf_staple_no_comm(gf, xl, mu));
}

inline void gf_ape_smear_no_comm(GaugeField& gf, const GaugeField& gf0,
                                 const double alpha)
{
  TIMER_VERBOSE("gf_ape_smear_no_comm");
  qassert(&gf != &gf0);
  const Geometry& geo = gf0.geo();
  gf.init(geo_resize(geo));
  qassert(is_matching_geo_mult(geo, gf.geo()));
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = gf.get_elems(xl);
    for (int mu = 0; mu < DIMN; ++mu) {
      v[mu] = gf_link_ape_smear_no_comm(gf0, xl, mu, alpha);
    }
  }
}

inline void gf_ape_smear(GaugeField& gf, const GaugeField& gf0,
                         const double alpha, const long steps = 1)
{
  TIMER_VERBOSE("gf_ape_smear");
  gf = gf0;
  GaugeField gf1;
  gf1.init(geo_resize(gf0.geo(), 1));
  for (long i = 0; i < steps; ++i) {
    gf1 = gf;
    refresh_expanded(gf1);
    gf_ape_smear_no_comm(gf, gf1, alpha);
  }
}

inline ColorMatrix gf_link_spatial_ape_smear_no_comm(const GaugeField& gf,
                                                     const Coordinate& xl,
                                                     const int mu,
                                                     const double alpha)
{
  const double multi = mu == 3 ? 6.0 : 4.0;
  return color_matrix_su_projection(
      (ComplexT)(1.0 - alpha) * gf.get_elem(xl, mu) +
      (ComplexT)(alpha / multi) * gf_spatial_staple_no_comm(gf, xl, mu));
}

inline void gf_spatial_ape_smear_no_comm(GaugeField& gf, const GaugeField& gf0,
                                         const double alpha)
{
  TIMER_VERBOSE("gf_spatial_ape_smear_no_comm");
  qassert(&gf != &gf0);
  const Geometry& geo = gf0.geo();
  gf.init(geo_resize(geo));
  qassert(is_matching_geo_mult(geo, gf.geo()));
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = gf.get_elems(xl);
    for (int mu = 0; mu < 3; ++mu) {
      // No need to smear the temperal link (mu == 3)
      v[mu] = gf_link_spatial_ape_smear_no_comm(gf0, xl, mu, alpha);
    }
  }
}

inline void gf_spatial_ape_smear(GaugeField& gf, const GaugeField& gf0,
                                 const double alpha, const long steps = 1)
{
  TIMER_VERBOSE("gf_spatial_ape_smear");
  gf = gf0;
  const Coordinate expan_left(1, 1, 1, 0);
  const Coordinate expan_right(1, 1, 1, 0);
  GaugeField gf1;
  gf1.init(geo_resize(gf0.geo(), expan_left, expan_right));
  for (long i = 0; i < steps; ++i) {
    gf1 = gf;
    refresh_expanded(gf1);
    gf_spatial_ape_smear_no_comm(gf, gf1, alpha);
  }
}

inline ColorMatrix gf_link_hyp_smear_3_no_comm(const GaugeField& gf,
                                               const Coordinate& xl,
                                               const int mu, const int nu,
                                               const int rho,
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
  ret = (ComplexT)(1.0 - alpha3) * gf.get_elem(xl, mu) +
        (ComplexT)(alpha3 / 2.0) * ret;
  return color_matrix_su_projection(ret);
}

inline ColorMatrix gf_link_hyp_smear_2_no_comm(const GaugeField& gf,
                                               const Coordinate& xl,
                                               const int mu, const int nu,
                                               const double alpha2,
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
  ret = (ComplexT)(1.0 - alpha2) * gf.get_elem(xl, mu) +
        (ComplexT)(alpha2 / 4.0) * ret;
  return color_matrix_su_projection(ret);
}

inline ColorMatrix gf_link_hyp_smear_1_no_comm(
    const GaugeField& gf, const Coordinate& xl, const int mu,
    const double alpha1, const double alpha2, const double alpha3)
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
  ret = (ComplexT)(1.0 - alpha1) * gf.get_elem(xl, mu) +
        (ComplexT)(alpha1 / 6.0) * ret;
  return color_matrix_su_projection(ret);
}

inline ColorMatrix gf_link_hyp_smear_no_comm(const GaugeField& gf,
                                             const Coordinate& xl, const int mu,
                                             const double alpha1,
                                             const double alpha2,
                                             const double alpha3)
{
  return gf_link_hyp_smear_1_no_comm(gf, xl, mu, alpha1, alpha2, alpha3);
}

inline void gf_hyp_smear_no_comm(GaugeField& gf, const GaugeField& gf0,
                                 const double alpha1, const double alpha2,
                                 const double alpha3)
{
  TIMER_VERBOSE("gf_hyp_smear_no_comm");
  qassert(&gf != &gf0);
  const Geometry& geo = gf0.geo();
  gf.init(geo_resize(geo));
  qassert(is_matching_geo_mult(geo, gf.geo()));
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = gf.get_elems(xl);
    for (int mu = 0; mu < DIMN; ++mu) {
      v[mu] = gf_link_hyp_smear_no_comm(gf0, xl, mu, alpha1, alpha2, alpha3);
    }
  }
}

inline void gf_hyp_smear(GaugeField& gf, const GaugeField& gf0,
                         const double alpha1, const double alpha2,
                         const double alpha3)
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

template <class T>
void smear_propagator(Propagator4dT<T>& prop, const GaugeFieldT<T>& gf1,
                      const double coef, const int step,
                      const CoordinateD& mom = CoordinateD(),
                      const bool smear_in_time_dir = false)
// gf1 is left_expanded and refreshed
// set_left_expanded_gauge_field(gf1, gf)
// prop is of normal size
{
  //TIMER_VERBOSE("==smear propagator");
  TIMER_FLOPS("smear_propagator");
  const int n_avg = smear_in_time_dir ? 8 : 6;
  const long long vGb = prop.geo().local_volume()*12*4;
  timer.flops += vGb * step * n_avg * (3*(3*6 + 2*2));
  if (0 == step) {
    return;
  }
  const Geometry& geo = prop.geo();
  const Geometry geo1 =
      smear_in_time_dir
          ? geo_resize(geo, 1)
          : geo_resize(geo, Coordinate(1, 1, 1, 0), Coordinate(1, 1, 1, 0));
  const int dir_limit = smear_in_time_dir ? 4 : 3;
  array<Complex, 8> mom_factors_v;
  box_acc<array<Complex, 8> > mom_factors(mom_factors_v); // (array<Complex, 8>());
  for (int i = 0; i < 8; ++i) {
    const int dir = i - 4;
    const double phase = dir >= 0 ? mom[dir] : -mom[-dir - 1];
    mom_factors()[i] = std::polar(coef / n_avg, -phase);
  }
  Propagator4dT<T> prop1;
  prop1.init(geo1);
  for (int i = 0; i < step; ++i) {
    prop1 = prop;
    refresh_expanded_1(prop1);
    qacc_for(index, geo.local_volume(), {
      const Coordinate xl = prop.geo().coordinate_from_index(index);
      WilsonMatrixT<T>& wm = prop.get_elem(xl);
      wm *= 1 - coef;
      for (int dir = -dir_limit; dir < dir_limit; ++dir) {
        const Coordinate xl1 = coordinate_shifts(xl, dir);
        ColorMatrixT<T> link =
            dir >= 0
                ? gf1.get_elem(xl, dir)
                : (ColorMatrixT<T>)matrix_adjoint(gf1.get_elem(xl1, -dir - 1));
        link *= mom_factors()[dir + 4];
        wm += link * prop1.get_elem(xl1);
      }
    });
  }
}

}  // namespace qlat
