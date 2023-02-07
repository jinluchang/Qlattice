#pragma once

#include <qlat-utils/mat-vec.h>

namespace qlat
{  //

const SpinMatrixT<>& get_gamma_matrix(const int mu);

WilsonMatrix g5_herm(const WilsonMatrix& m);

void set_zero(ColorMatrix& x);
void set_zero(SpinMatrix& x);
void set_zero(WilsonMatrix& x);
void set_zero(NonRelWilsonMatrix& x);
void set_zero(IsospinMatrix& x);
void set_zero(WilsonVector& x);

Vector<Complex> get_data(const ColorMatrix& x);
Vector<Complex> get_data(const SpinMatrix& x);
Vector<Complex> get_data(const WilsonMatrix& x);
Vector<Complex> get_data(const NonRelWilsonMatrix& x);
Vector<Complex> get_data(const IsospinMatrix& x);
Vector<Complex> get_data(const WilsonVector& x);

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
