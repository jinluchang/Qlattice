#pragma once

#include <qlat-utils/mat-vec.h>
#include <qlat-utils/matrix.h>
#include <qlat-utils/vector.h>

namespace qlat
{  //

const SpinMatrixT<>& get_gamma_matrix(const int mu);

qacc WilsonMatrix g5_herm(const WilsonMatrix& m)
{
  const box_acc<SpinMatrixConstants>& smc = get_spin_matrix_constants();
  const SpinMatrix& gamma5 = smc().gamma5;
  const WilsonMatrix ret = gamma5 * (WilsonMatrix)matrix_adjoint(m) * gamma5;
  return ret;
}

qacc void set_zero(ColorMatrix& x) { set_zero<3, Complex>(x); }
qacc void set_zero(SpinMatrix& x) { set_zero<4, Complex>(x); }
qacc void set_zero(WilsonMatrix& x) { set_zero<12, Complex>(x); }
qacc void set_zero(NonRelWilsonMatrix& x) { set_zero<6, Complex>(x); }
qacc void set_zero(IsospinMatrix& x) { set_zero<2, Complex>(x); }
qacc void set_zero(WilsonVector& x) { set_zero<12, Complex>(x); }

qacc Vector<Complex> get_data(const ColorMatrix& x)
{
  return get_data<3, Complex>(x);
}
qacc Vector<Complex> get_data(const SpinMatrix& x)
{
  return get_data<4, Complex>(x);
}
qacc Vector<Complex> get_data(const WilsonMatrix& x)
{
  return get_data<12, Complex>(x);
}
qacc Vector<Complex> get_data(const NonRelWilsonMatrix& x)
{
  return get_data<6, Complex>(x);
}
qacc Vector<Complex> get_data(const IsospinMatrix& x)
{
  return get_data<2, Complex>(x);
}
qacc Vector<Complex> get_data(const WilsonVector& x)
{
  return get_data<12, Complex>(x);
}

qacc Complex mat_tr(const WilsonMatrix& m) { return matrix_trace(m); }
qacc Complex mat_tr(const SpinMatrix& m) { return matrix_trace(m); }
qacc Complex mat_tr(const WilsonMatrix& m1, const WilsonMatrix& m2)
{
  return matrix_trace(m1, m2);
}
qacc Complex mat_tr(const WilsonMatrix& m1, const SpinMatrix& m2)
{
  return matrix_trace(m1, m2);
}
qacc Complex mat_tr(const SpinMatrix& m1, const WilsonMatrix& m2)
{
  return matrix_trace(m1, m2);
}
qacc Complex mat_tr(const SpinMatrix& m1, const SpinMatrix& m2)
{
  return matrix_trace(m1, m2);
}

qacc SpinMatrix operator*(const SpinMatrix& m1, const SpinMatrix& m2)
{
  return operator*<4, Complex>(m1, m2);
}

qacc SpinMatrix operator*(const Complex& a, const SpinMatrix& m)
{
  return operator*<4, Complex>(a, m);
}

qacc SpinMatrix operator*(const SpinMatrix& m, const Complex& a)
{
  return operator*<4, Complex>(m, a);
}

qacc WilsonMatrix operator*(const WilsonMatrix& m1, const WilsonMatrix& m2)
{
  return operator*<12, Complex>(m1, m2);
}

qacc WilsonMatrix operator*(const Complex& a, const WilsonMatrix& m)
{
  return operator*<12, Complex>(a, m);
}

qacc WilsonMatrix operator*(const WilsonMatrix& m, const Complex& a)
{
  return operator*<12, Complex>(m, a);
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
  return matrix_adjoint<3, Complex>(m);
}

qacc SpinMatrix matrix_adjoint(const SpinMatrix& m)
{
  return matrix_adjoint<4, Complex>(m);
}

qacc WilsonMatrix matrix_adjoint(const WilsonMatrix& m)
{
  return matrix_adjoint<12, Complex>(m);
}

qacc IsospinMatrix matrix_adjoint(const IsospinMatrix& m)
{
  return matrix_adjoint<2, Complex>(m);
}

}  // namespace qlat
