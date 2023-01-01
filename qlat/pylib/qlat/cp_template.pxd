from qlat_utils.cp cimport *
from . cimport everything as cqlat

cdef class Geometry:

    cdef cqlat.Geometry xx

    cdef readonly long cdata

cdef class ElemType:

    pass

cdef class Field:

    pass

cdef class SelectedField:

    pass

cdef class SelectedPoints:

    pass
