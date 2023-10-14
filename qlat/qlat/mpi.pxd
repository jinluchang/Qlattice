from qlat_utils.coordinate cimport *

cdef extern from "qlat/mpi.h" namespace "qlat":

    int begin(const int id_node, const Coordinate& size_node, const int color) except +
    int end(const bool is_preserving_cache) except +
