from qlat.everything cimport *

cdef extern from "qlat-cps/qlat-cps.h" namespace "qlat":

    void begin_with_cps(const std_vector[std_string]& sargs,
                        const Coordinate& total_site) except +

    void end_with_cps(const bool is_preserving_cache) except +

    void save_cps_prop_double(const Field[WilsonMatrix]& prop, const std_string& path) except +

    void load_cps_prop_double(Field[WilsonMatrix]& prop, const std_string& path) except +
