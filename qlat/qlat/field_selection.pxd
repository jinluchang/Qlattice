from . cimport everything as cqlat
from .geometry cimport Geometry

cdef class PointsSelection:

    cdef cqlat.PointsSelection xx
    cdef readonly cqlat.Long cdata
    cdef readonly int view_count

cdef class FieldSelection:

    cdef cqlat.FieldSelection xx
    cdef readonly cqlat.Long cdata
    cdef readonly int view_count
