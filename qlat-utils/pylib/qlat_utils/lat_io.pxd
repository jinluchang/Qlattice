from libcpp.string cimport string as std_string
from .complex cimport *
from .vector cimport *

cdef extern from "qlat-utils/lat-io.h" namespace "qlat":

    cdef cppclass LatData:
        LatData()
        void load(const std_string& fn)
        void save(const std_string& fn)
        const LatData& operator=(const LatData& ld)

    Vector[double] get_data(const LatData& ld)
