#include <qlat/vector_utils/utils_Matrix_prod.h>

namespace qlat
{

void matrix_prodP(qlat::ComplexT<double>** a, qlat::ComplexT<double>** b,
                  qlat::ComplexT<double>** c, const Long m, const Long n,
                  const Long w, const Long L, bool Conj, bool trans, bool GPU,
                  QBOOL dummy)
{
  matrix_prodPT(a, b, c, m, n, w, L, Conj, trans, GPU, dummy);
}

void matrix_prodP(qlat::ComplexT<RealF>** a, qlat::ComplexT<RealF>** b,
                  qlat::ComplexT<RealF>** c, const Long m, const Long n,
                  const Long w, const Long L, bool Conj, bool trans, bool GPU,
                  QBOOL dummy)
{
  matrix_prodPT(a, b, c, m, n, w, L, Conj, trans, GPU, dummy);
}

void matrix_prod(qlat::ComplexT<double>* A, qlat::ComplexT<double>* B,
                 qlat::ComplexT<double>* C, const Long m, const Long n,
                 const Long w, const Long L, bool Conj, bool trans, bool GPU,
                 QBOOL dummy)
{
  matrix_prodT(A, B, C, m, n, w, L, Conj, trans, GPU, dummy);
}

void matrix_prod(qlat::ComplexT<RealF>* A, qlat::ComplexT<RealF>* B,
                 qlat::ComplexT<RealF>* C, const Long m, const Long n,
                 const Long w, const Long L, bool Conj, bool trans, bool GPU,
                 QBOOL dummy)
{
  matrix_prodT(A, B, C, m, n, w, L, Conj, trans, GPU, dummy);
}

}  // namespace qlat
