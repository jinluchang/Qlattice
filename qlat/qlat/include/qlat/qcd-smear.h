#pragma once

#include <qlat-utils/matrix.h>
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
    y(i1, i1) = ComplexD(p0, -p3);
    y(i2, i2) = ComplexD(p0, p3);
    y(i1, i2) = ComplexD(-p2, -p1);
    y(i2, i1) = ComplexD(p2, -p1);
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

void gf_ape_smear(GaugeField& gf, const GaugeField& gf0, const double alpha,
                  const Long steps = 1);

void gf_spatial_ape_smear(GaugeField& gf, const GaugeField& gf0,
                          const double alpha, const Long steps = 1);

void gf_hyp_smear(GaugeField& gf, const GaugeField& gf0, const double alpha1,
                  const double alpha2, const double alpha3);

template <class T>
void prop_smear(Propagator4dT<T>& prop, const GaugeFieldT<T>& gf1,
                const double coef, const int step,
                const CoordinateD& mom = CoordinateD(),
                const bool smear_in_time_dir = false)
// gf1 is left_expanded and refreshed
// set_left_expanded_gauge_field(gf1, gf)
// prop is of normal size
{
  TIMER_FLOPS("prop_smear");
  const int n_avg = smear_in_time_dir ? 8 : 6;
  const Long vGb = prop.geo().local_volume() * 12 * 4;
  timer.flops += vGb * step * n_avg * (3 * (3 * 6 + 2 * 2));
  if (0 == step) {
    return;
  }
  const Geometry& geo = prop.geo();
  const Geometry geo1 =
      smear_in_time_dir
          ? geo_resize(geo, 1)
          : geo_resize(geo, Coordinate(1, 1, 1, 0), Coordinate(1, 1, 1, 0));
  const int dir_limit = smear_in_time_dir ? 4 : 3;
  array<ComplexD, 8> mom_factors_v;
  box_acc<array<ComplexD, 8>> mom_factors(
      mom_factors_v);  // (array<ComplexD, 8>());
  for (int i = 0; i < 8; ++i) {
    const int dir = i - 4;
    const double phase = dir >= 0 ? mom[dir] : -mom[-dir - 1];
    mom_factors()[i] = qpolar(coef / n_avg, -phase);
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

#ifdef QLAT_INSTANTIATE_SMEAR
#define QLAT_EXTERN
#else
#define QLAT_EXTERN extern
#endif

QLAT_EXTERN template void prop_smear<Real>(Propagator4d&, const GaugeField&,
                                           const double, const int,
                                           const CoordinateD&, const bool);

#undef QLAT_EXTERN

}  // namespace qlat
