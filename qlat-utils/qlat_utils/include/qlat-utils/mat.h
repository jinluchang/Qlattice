#pragma once

#include <qlat-utils/mat-vec.h>
#include <qlat-utils/matrix.h>
#include <qlat-utils/vector.h>

namespace qlat
{  //

const SpinMatrixT<RealD>& get_gamma_matrix(const Int mu);

void benchmark_matrix_functions(const Long count = 128);

qacc WilsonMatrix g5_herm(const WilsonMatrix& m)
{
  const box<SpinMatrixConstants>& smc = get_spin_matrix_constants();
  const SpinMatrix& gamma5 = smc().gamma5;
  const WilsonMatrix ret = gamma5 * (WilsonMatrix)matrix_adjoint(m) * gamma5;
  return ret;
}

qacc void set_zero(ColorMatrix& x) { set_zero<3, ComplexD>(x); }
qacc void set_zero(SpinMatrix& x) { set_zero<4, ComplexD>(x); }
qacc void set_zero(WilsonMatrix& x) { set_zero<12, ComplexD>(x); }
qacc void set_zero(NonRelWilsonMatrix& x) { set_zero<6, ComplexD>(x); }
qacc void set_zero(IsospinMatrix& x) { set_zero<2, ComplexD>(x); }
qacc void set_zero(WilsonVector& x) { set_zero<12, ComplexD>(x); }

qacc SpinMatrix operator*(const SpinMatrix& m1, const SpinMatrix& m2)
{
  return operator*<4, ComplexD>(m1, m2);
}

qacc SpinMatrix operator*(const ComplexD& a, const SpinMatrix& m)
{
  return operator*<4, ComplexD>(a, m);
}

qacc SpinMatrix operator*(const SpinMatrix& m, const ComplexD& a)
{
  return operator*<4, ComplexD>(m, a);
}

qacc WilsonMatrix operator*(const WilsonMatrix& m1, const WilsonMatrix& m2)
{
  return operator*<12, ComplexD>(m1, m2);
}

qacc WilsonMatrix operator*(const ComplexD& a, const WilsonMatrix& m)
{
  return operator*<12, ComplexD>(a, m);
}

qacc WilsonMatrix operator*(const WilsonMatrix& m, const ComplexD& a)
{
  return operator*<12, ComplexD>(m, a);
}

qacc WilsonMatrix operator*(const SpinMatrix& m1, const WilsonMatrix& m2)
{
  return operator*<Real>(m1, m2);
}

qacc WilsonMatrix operator*(const WilsonMatrix& m1, const SpinMatrix& m2)
{
  return operator*<Real>(m1, m2);
}

qacc ColorMatrix matrix_adjoint(const ColorMatrix& m)
{
  return matrix_adjoint<3, ComplexD>(m);
}

qacc SpinMatrix matrix_adjoint(const SpinMatrix& m)
{
  return matrix_adjoint<4, ComplexD>(m);
}

qacc WilsonMatrix matrix_adjoint(const WilsonMatrix& m)
{
  return matrix_adjoint<12, ComplexD>(m);
}

qacc IsospinMatrix matrix_adjoint(const IsospinMatrix& m)
{
  return matrix_adjoint<2, ComplexD>(m);
}

qacc ComplexD epsilon_contraction(const Int v_s1, const Int b_s1,
                                  const Int v_s2, const Int b_s2,
                                  const Int v_s3, const Int b_s3,
                                  const WilsonMatrix& wm1,
                                  const WilsonMatrix& wm2,
                                  const WilsonMatrix& wm3)
{
  return epsilon_contraction<RealD>(v_s1, b_s1, v_s2, b_s2, v_s3, b_s3, wm1,
                                    wm2, wm3);
}

}  // namespace qlat
