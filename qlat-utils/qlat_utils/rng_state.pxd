from . cimport everything as cqlat_utils

cdef class RngState:
    cdef cqlat_utils.RngState xx
    cdef readonly long cdata
