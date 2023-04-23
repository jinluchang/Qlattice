#include <qlat-utils/mat.h>
#include <qlat-utils/matrix.h>
#include <qlat-utils/vector.h>

namespace qlat
{  //

const SpinMatrixT<>& get_gamma_matrix(const int mu)
// CPS's convention gamma matrices
{
  return SpinMatrixConstantsT<>::get_cps_gamma(mu);
}

}  // namespace qlat
