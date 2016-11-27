#pragma once

#include <qlat/matrix.h>

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

QLAT_END_NAMESPACE
