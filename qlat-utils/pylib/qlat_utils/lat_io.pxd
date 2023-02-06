from .complex cimport *
from .vector cimport *
from .timer cimport *

cdef extern from "qlat-utils/lat-io.h" namespace "qlat":

    cdef cppclass LatDim:
        std_string name
        long size
        std_vector[std_string] indices
        LatDim()

    LatDim lat_dim_re_im()
    LatDim lat_dim_number(const std_string& name, const long start, const long end)
    LatDim lat_dim_number(const std_string& name, const long start, const long end, const long inc)
    LatDim lat_dim_string(const std_string& name, const std_vector[std_string]& indices)

    ctypedef std_vector[LatDim] LatInfo

    cdef cppclass LatData:
        LatInfo info
        std_vector[double] res
        LatData()
        const LatData& operator=(const LatData& ld)
        void load(const std_string& fn)
        void save(const std_string& fn)
        bool is_complex()
        int ndim()
        double* data()

    LatData operator*(const Complex& a, const LatData& ld)
    LatData operator*(const double a, const LatData& ld)
    LatData operator*(const LatData& ld, const Complex& a)
    LatData operator*(const LatData& ld, const double a)
    LatData operator+(const LatData& ld1, const LatData& ld2)
    LatData operator-(const LatData& ld1, const LatData& ld2)

    bool is_matching(const LatData& ld1, const LatData& ld2)
    long lat_data_size(LatData& ld)
    void lat_data_alloc(LatData& ld)

    Vector[double] get_data(const LatData& x)
    void set_zero(LatData& x)
    std_string show(const LatData& x)
    void clear(LatData& x)
    double qnorm(const LatData& x)
