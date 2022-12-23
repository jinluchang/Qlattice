from libcpp cimport bool

cdef extern from "qlat-grid/qlat-grid.h" namespace "qlat":
    cdef void grid_end()
    cdef void grid_end(const bool is_preserving_cache)
