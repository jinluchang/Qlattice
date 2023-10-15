from . cimport everything as cqlat_utils

cdef class Coordinate:
    cdef cqlat_utils.Coordinate xx

cdef class CoordinateD:
    cdef cqlat_utils.CoordinateD xx
