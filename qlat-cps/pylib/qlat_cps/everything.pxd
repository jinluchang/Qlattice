from qlat.everything cimport *

cdef extern from "qlat-cps/qlat-cps.h" namespace "qlat":

    void begin_with_cps(const std_vector[std_string]& sargs,
                        const Coordinate& total_site)

    void end_with_cps(const bool is_preserving_cache)
