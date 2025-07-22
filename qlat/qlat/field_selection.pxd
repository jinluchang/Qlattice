from . cimport everything as cqlat
from .geometry cimport Geometry

cdef class SelectedShufflePlan:

    cdef cqlat.SelectedShufflePlan xx
    cdef public list psel_src_list
    cdef public list psel_dst_list

cdef class PointsSelection:

    cdef cqlat.PointsSelection xx
    cdef readonly cqlat.Long cdata
    cdef readonly int view_count

cdef class FieldSelection:

    cdef cqlat.FieldSelection xx
    cdef readonly cqlat.Long cdata
    cdef readonly int view_count
