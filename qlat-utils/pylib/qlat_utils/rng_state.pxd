from libcpp.string cimport string as std_string
from .complex cimport *

cdef extern from "qlat-utils/rng-state.h" namespace "qlat":

    cdef cppclass RngState:
        RngState()
        RngState(const std_string& seed)
        RngState(const RngState& rs0, const std_string& sindex)
        const RngState& operator=(const RngState& rs)
        RngState split(const std_string& sindex)
        RngState newtype(const unsigned long type)

    uint64_t rand_gen(RngState& rs)
    double u_rand_gen(RngState& rs, const double upper, const double lower)
    double g_rand_gen(RngState& rs, const double center, const double sigma)
