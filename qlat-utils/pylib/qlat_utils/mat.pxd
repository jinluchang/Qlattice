from qlat_utils.lat_io cimport *

cdef extern from "qlat-utils/mat-vec.h" namespace "qlat":

    cdef cppclass ColorMatrix:
        ColorMatrix()
        const ColorMatrix& operator=(const ColorMatrix& m)
        Complex* data()

    cdef cppclass SpinMatrix:
        SpinMatrix()
        const SpinMatrix& operator=(const SpinMatrix& m)
        Complex* data()

    cdef cppclass WilsonMatrix:
        WilsonMatrix()
        const WilsonMatrix& operator=(const WilsonMatrix& m)
        Complex* data()

    cdef cppclass NonRelWilsonMatrix:
        NonRelWilsonMatrix()
        const NonRelWilsonMatrix& operator=(const NonRelWilsonMatrix& m)
        Complex* data()

    cdef cppclass IsospinMatrix:
        IsospinMatrix()
        const IsospinMatrix& operator=(const IsospinMatrix& m)
        Complex* data()

    cdef cppclass WilsonVector:
        WilsonVector()
        const WilsonVector& operator=(const WilsonVector& m)
        Complex* data()

cdef extern from "qlat-utils/lib/mat.h" namespace "qlat":

    void set_zero(ColorMatrix& x)
    void set_zero(SpinMatrix& x)
    void set_zero(WilsonMatrix& x)
    void set_zero(NonRelWilsonMatrix& x)
    void set_zero(IsospinMatrix& x)
    void set_zero(WilsonVector& x)

    Vector[Complex] get_data(const ColorMatrix& x)
    Vector[Complex] get_data(const SpinMatrix& x)
    Vector[Complex] get_data(const WilsonMatrix& x)
    Vector[Complex] get_data(const NonRelWilsonMatrix& x)
    Vector[Complex] get_data(const IsospinMatrix& x)
    Vector[Complex] get_data(const WilsonVector& x)

    Complex mat_tr(const SpinMatrix& m)
    Complex mat_tr(const WilsonMatrix& m)
    Complex mat_tr(const WilsonMatrix& m1, const WilsonMatrix& m2)
    Complex mat_tr(const WilsonMatrix& m1, const SpinMatrix& m2)
    Complex mat_tr(const SpinMatrix& m1, const WilsonMatrix& m2)
    Complex mat_tr(const SpinMatrix& m1, const SpinMatrix& m2)

    const SpinMatrix& get_gamma_matrix(const int mu)

    WilsonMatrix g5_herm(const WilsonMatrix& m)

    SpinMatrix operator*(const Complex& a, const SpinMatrix& m)
    SpinMatrix operator*(const SpinMatrix& m, const Complex& a)
    SpinMatrix operator*(const SpinMatrix& m1, const SpinMatrix& m2)

    WilsonMatrix operator*(const Complex& a, const WilsonMatrix& m)
    WilsonMatrix operator*(const WilsonMatrix& m, const Complex& a)
    WilsonMatrix operator*(const WilsonMatrix& m1, const WilsonMatrix& m2)
    WilsonMatrix operator*(const SpinMatrix& m1, const WilsonMatrix& m2)
    WilsonMatrix operator*(const WilsonMatrix& m1, const SpinMatrix& m2)
