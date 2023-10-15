from . cimport everything as cqlat_utils

cdef class LatData:
    cdef cqlat_utils.LatData xx
    cdef readonly long cdata
