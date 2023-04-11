from qlat_utils.everything cimport *

cdef extern from "qlat-grid/qlat-grid.h" namespace "qlat":
    void begin_with_grid(const std_vector[std_string]& sargs,
                         const std_vector[Coordinate]& node_size_list)

    void end_with_grid(const bool is_preserving_cache)
