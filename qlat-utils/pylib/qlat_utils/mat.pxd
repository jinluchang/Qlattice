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

cdef extern from "qlat-utils/mat.h" namespace "qlat":

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

    Complex matrix_trace(const SpinMatrix& m)
    Complex matrix_trace(const ColorMatrix& m)
    Complex matrix_trace(const WilsonMatrix& m)
    Complex matrix_trace(const WilsonMatrix& m1, const WilsonMatrix& m2)
    Complex matrix_trace(const WilsonMatrix& m1, const SpinMatrix& m2)
    Complex matrix_trace(const SpinMatrix& m1, const WilsonMatrix& m2)
    Complex matrix_trace(const SpinMatrix& m1, const SpinMatrix& m2)
    Complex matrix_trace(const WilsonMatrix& m1, const ColorMatrix& m2)
    Complex matrix_trace(const ColorMatrix& m1, const WilsonMatrix& m2)
    Complex matrix_trace(const ColorMatrix& m1, const ColorMatrix& m2)

    const SpinMatrix& get_gamma_matrix(const int mu)

    void benchmark_matrix_functions(const long count)

    WilsonMatrix g5_herm(const WilsonMatrix& m)

    SpinMatrix operator*(const Complex& a, const SpinMatrix& m)
    SpinMatrix operator*(const SpinMatrix& m, const Complex& a)
    SpinMatrix operator*(const SpinMatrix& m1, const SpinMatrix& m2)

    ColorMatrix operator*(const Complex& a, const ColorMatrix& m)
    ColorMatrix operator*(const ColorMatrix& m, const Complex& a)
    ColorMatrix operator*(const ColorMatrix& m1, const ColorMatrix& m2)

    WilsonMatrix operator*(const Complex& a, const WilsonMatrix& m)
    WilsonMatrix operator*(const WilsonMatrix& m, const Complex& a)
    WilsonMatrix operator*(const WilsonMatrix& m1, const WilsonMatrix& m2)
    WilsonMatrix operator*(const SpinMatrix& m1, const WilsonMatrix& m2)
    WilsonMatrix operator*(const WilsonMatrix& m1, const SpinMatrix& m2)
    WilsonMatrix operator*(const ColorMatrix& m1, const WilsonMatrix& m2)
    WilsonMatrix operator*(const WilsonMatrix& m1, const ColorMatrix& m2)
