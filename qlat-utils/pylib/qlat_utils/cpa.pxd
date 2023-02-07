from . cimport everything as cqlat_utils

cdef class Coordinate:
    cdef cqlat_utils.Coordinate xx

cdef class RngState:
    cdef cqlat_utils.RngState xx
    cdef readonly long cdata

cdef class LatData:
    cdef cqlat_utils.LatData xx
    cdef readonly long cdata
