from libcpp.vector cimport vector as std_vector
from libcpp.string cimport string as std_string
from cpython.ref cimport PyObject

from .rng_state cimport *

cdef extern from "qlat-utils/cache.h" namespace "qlat":

    void random_permute[T](std_vector[T]& vec, const RngState& rs)
    std_vector[std_string] get_all_caches_info()
    void clear_all_caches()
    void displayln_malloc_stats()
