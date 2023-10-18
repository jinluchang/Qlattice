from . cimport everything as cqlat

cdef class Geometry:

    cdef cqlat.Geometry xx

    cdef readonly long cdata
