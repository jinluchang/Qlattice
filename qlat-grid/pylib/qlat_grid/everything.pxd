from libcpp cimport bool

cdef extern from "qlat-grid/qlat-grid.h" namespace "qlat":
    void grid_end()
    void grid_end(const bool is_preserving_cache)
