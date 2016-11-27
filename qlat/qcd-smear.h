#pragma once

#include <qlat/matrix.h>
#include <qlat/qcd.h>
#include <qlat/qcd-utils.h>

QLAT_START_NAMESPACE

ColorMatrix color_matrix_sub_inverse(const ColorMatrix& x, const int ind)
  // get su2 submatrix of x and return the su3 matrix that
  // has the inverse of this matrix in the relevant row and column
{
  const int su2_index[][3] = {
    {0,1,2},
    {0,2,1},
    {1,2,0}
  };
  const int i1 = su2_index[ind][0];
  const int i2 = su2_index[ind][1];
  const int zero_rc = su2_index[ind][2];
  // project onto SU(2)
  double p0 = x(i1,i1).real() + x(i2,i2).real();
  double p1 = x(i1,i2).imag() + x(i2,i1).imag();
  double p2 = x(i1,i2).real() - x(i2,i1).real();
  double p3 = x(i1,i1).imag() - x(i2,i2).imag();
  const double psqr = sqrt(p0*p0 + p1*p1 + p2*p2 + p3*p3);
  double ipsqr;
  if (psqr == 0.0) {
    ipsqr = 1.0;
  } else {
    ipsqr = 1.0/psqr;
  }
  p0 *= ipsqr;
  p1 *= ipsqr;
  p2 *= ipsqr;
  p3 *= ipsqr;
  ColorMatrix y;
  set_unit(y);
  // fill with inverse
  y(i1,i1) = Complex( p0,-p3);
  y(i2,i2) = Complex( p0, p3);
  y(i1,i2) = Complex(-p2,-p1);
  y(i2,i1) = Complex( p2,-p1);
  return y;
}

ColorMatrix color_matrix_su_projection(const ColorMatrix& x, const double tolerance = 1.0e-8)
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
  for (int i=0;i<max_iter;i++) {
    // loop over su2 subgroups
    double diff = 0.0;
    for (int j=0;j<3;j++) {
      const ColorMatrix inv = color_matrix_sub_inverse(tmp, j);
      // y  .DotMEqual( inv, ycopy );
      y = inv * y;
      // tmp.DotMEqual( y, xdag );
      tmp = y * xdag;
      const double tr = matrix_trace(tmp).real();
      const double dtr = tr-old_tr;
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
  return y;
}

inline void gf_ape_smear_no_comm(GaugeField& gf, const GaugeField& gf0, const double alpha)
{
  TIMER_VERBOSE("gf_ape_smear_no_comm");
  qassert(&gf != &gf0);
  const Geometry& geo = gf0.geo;
  gf.init(geo_resize(geo));
  qassert(is_matching_geo_mult(geo, gf.geo));
#pragma omp parallel for
  for (long index = 0; index < geo.local_volume(); ++index) {
    const Coordinate xl = geo.coordinate_from_index(index);
    Vector<ColorMatrix> v = gf.get_elems(xl);
    const Vector<ColorMatrix> v0 = gf0.get_elems_const(xl);
    for (int mu = 0; mu < DIM; ++mu) {
      v[mu] = (1.0-alpha) * v0[mu] + alpha/6.0 * gf_staple_no_comm(gf0, xl, mu);
      v[mu] = color_matrix_su_projection(v[mu]);
    }
  }
}

inline void gf_ape_smear(GaugeField& gf, const GaugeField& gf0, const double alpha)
{
  TIMER_VERBOSE("gf_ape_smear");
  GaugeField gf1;
  gf1.init(geo_resize(gf0.geo, 1));
  gf1 = gf0;
  refresh_expanded(gf1);
  gf_ape_smear_no_comm(gf, gf1, alpha);
}

QLAT_END_NAMESPACE
