#include <qlat/core.h>

namespace qlat
{  //

WilsonMatrix g5_herm(const WilsonMatrix& m);

Complex mat_tr(const WilsonMatrix& m);
Complex mat_tr(const SpinMatrix& m);
Complex mat_tr(const WilsonMatrix& m1, const WilsonMatrix& m2);
Complex mat_tr(const WilsonMatrix& m1, const SpinMatrix& m2);
Complex mat_tr(const SpinMatrix& m1, const WilsonMatrix& m2);
Complex mat_tr(const SpinMatrix& m1, const SpinMatrix& m2);

}  // namespace qlat
