from libcpp.complex import complex[double] as Complex

cdef extern from "qlat/matrix.h" namespace "qlat":

    cdef cppclass SpinMatrix:
        SpinMatrix()
        const SpinMatrix& operator=(const SpinMatrix& m)
        const SpinMatrix& operator*(const Complex& a, const SpinMatrix& m)
        const SpinMatrix& operator*(const SpinMatrix& m, const Complex& a)
        const SpinMatrix& operator*(const SpinMatrix& m1, const SpinMatrix& m2)

    cdef cppclass WilsonMatrix:
        WilsonMatrix()
        const WilsonMatrix& operator=(const WilsonMatrix& m)
        const WilsonMatrix& operator*(const Complex& a, const WilsonMatrix& m)
        const WilsonMatrix& operator*(const WilsonMatrix& m, const Complex& a)
        const WilsonMatrix& operator*(const WilsonMatrix& m1, const WilsonMatrix& m2)
        const WilsonMatrix& operator*(const SpinMatrix& m1, const WilsonMatrix& m2)
        const WilsonMatrix& operator*(const WilsonMatrix& m1, const SpinMatrix& m2)

cdef extern from "qlat/lib/mat.h" namespace "qlat":

    cdef void set_zero(SpinMatrix& m)
    cdef void set_zero(WilsonMatrix& m)

    cdef Complex mat_tr(const SpinMatrix& m)
    cdef Complex mat_tr(const WilsonMatrix& m)
    cdef Complex mat_tr(const WilsonMatrix& m1, const WilsonMatrix& m2)
    cdef Complex mat_tr(const WilsonMatrix& m1, const SpinMatrix& m2)
    cdef Complex mat_tr(const SpinMatrix& m1, const WilsonMatrix& m2)
    cdef Complex mat_tr(const SpinMatrix& m1, const SpinMatrix& m2)

    cdef WilsonMatrix g5_herm(const WilsonMatrix& m)
