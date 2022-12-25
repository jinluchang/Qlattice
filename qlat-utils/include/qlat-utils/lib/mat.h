#include <qlat-utils/matrix.h>

namespace qlat
{  //

WilsonMatrix g5_herm(const WilsonMatrix& m);

Complex mat_tr(const WilsonMatrix& m);
Complex mat_tr(const SpinMatrix& m);
Complex mat_tr(const WilsonMatrix& m1, const WilsonMatrix& m2);
Complex mat_tr(const WilsonMatrix& m1, const SpinMatrix& m2);
Complex mat_tr(const SpinMatrix& m1, const WilsonMatrix& m2);
Complex mat_tr(const SpinMatrix& m1, const SpinMatrix& m2);

SpinMatrix operator*(const SpinMatrix& m1, const SpinMatrix& m2);
SpinMatrix operator*(const Complex& a, const SpinMatrix& m);
SpinMatrix operator*(const SpinMatrix& m, const Complex& a);

WilsonMatrix operator*(const WilsonMatrix& m1, const WilsonMatrix& m2);
WilsonMatrix operator*(const Complex& a, const WilsonMatrix& m);
WilsonMatrix operator*(const WilsonMatrix& m, const Complex& a);
WilsonMatrix operator*(const SpinMatrix& m1, const WilsonMatrix& m2);
WilsonMatrix operator*(const WilsonMatrix& m1, const SpinMatrix& m2);

}  // namespace qlat
