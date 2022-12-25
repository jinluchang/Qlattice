#include <qlat-utils/lib/mat.h>
#include <qlat-utils/matrix.h>

namespace qlat
{  //

WilsonMatrix g5_herm(const WilsonMatrix& m)
{
  const box_acc<SpinMatrixConstants>& smc = get_spin_matrix_constants();
  const SpinMatrix& gamma5 = smc().gamma5;
  const WilsonMatrix ret = gamma5 * (WilsonMatrix)matrix_adjoint(m) * gamma5;
  return ret;
}

Complex mat_tr(const WilsonMatrix& m) { return matrix_trace(m); }
Complex mat_tr(const SpinMatrix& m) { return matrix_trace(m); }
Complex mat_tr(const WilsonMatrix& m1, const WilsonMatrix& m2)
{
  return matrix_trace(m1, m2);
}
Complex mat_tr(const WilsonMatrix& m1, const SpinMatrix& m2)
{
  return matrix_trace(m1, m2);
}
Complex mat_tr(const SpinMatrix& m1, const WilsonMatrix& m2)
{
  return matrix_trace(m1, m2);
}
Complex mat_tr(const SpinMatrix& m1, const SpinMatrix& m2)
{
  return matrix_trace(m1, m2);
}

SpinMatrix operator*(const SpinMatrix& m1, const SpinMatrix& m2)
{
  return operator*<4, Complex >(m1, m2);
}

SpinMatrix operator*(const Complex& a, const SpinMatrix& m)
{
  return operator*<4, Complex >(a, m);
}

SpinMatrix operator*(const SpinMatrix& m, const Complex& a)
{
  return operator*<4, Complex >(m, a);
}

WilsonMatrix operator*(const WilsonMatrix& m1, const WilsonMatrix& m2)
{
  return operator*<12, Complex >(m1, m2);
}

WilsonMatrix operator*(const Complex& a, const WilsonMatrix& m)
{
  return operator*<12, Complex >(a, m);
}

WilsonMatrix operator*(const WilsonMatrix& m, const Complex& a)
{
  return operator*<12, Complex >(m, a);
}

WilsonMatrix operator*(const SpinMatrix& m1, const WilsonMatrix& m2)
{
  return operator*<Real>(m1, m2);
}

WilsonMatrix operator*(const WilsonMatrix& m1, const SpinMatrix& m2)
{
  return operator*<Real>(m1, m2);
}

}  // namespace qlat
