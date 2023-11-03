from . cimport everything as cqlat
from .field_selection cimport FieldSelection, PointsSelection

cdef class FieldBase:

    cdef readonly cqlat.Long cdata
    cdef readonly int view_count

cdef class SelectedFieldBase:

    cdef readonly cqlat.Long cdata
    cdef readonly int view_count
    cdef public FieldSelection fsel

cdef class SelectedPointsBase:

    cdef readonly cqlat.Long cdata
    cdef readonly int view_count
    cdef public PointsSelection psel
