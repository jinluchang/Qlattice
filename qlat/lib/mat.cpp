#include <qlat/lib/mat.h>
#include <qlat/matrix.h>

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

}  // namespace qlat
