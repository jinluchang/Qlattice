#pragma once

#include "data-load.h"

namespace qlat
{  //

inline void set_local_va_current(Vector<WilsonMatrix> fc,
                                 const WilsonMatrix& prop_a,
                                 const WilsonMatrix& prop_b)
// ->- prop_a ->- op ->- inv_prop_b ->-
{
  const std::array<SpinMatrix, 8>& va_ms = get_va_matrices();
  const SpinMatrix& gamma5 = SpinMatrixConstants::get_gamma5();
  const WilsonMatrix inv_prop_b =
      gamma5 * (WilsonMatrix)matrix_adjoint(prop_b) * gamma5;
  for (int m = 0; m < 8; ++m) {
    fc[m] = inv_prop_b * va_ms[m] * prop_a;
  }
}

inline int tsep_op_wall_src(const std::string& job_tag)
// parameter
{
  if (job_tag == "24D" or job_tag == "32D" or job_tag == "24DH") {
    return 8;
  } else if (job_tag == "32Dfine") {
    return 10;
  } else if (job_tag == "48I") {
    return 12;
  } else if (job_tag == "64I") {
    return 18;
  } else {
    qassert(false);
  }
  return 8;
}

}  // namespace qlat
