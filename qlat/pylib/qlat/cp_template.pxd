from qlat_utils.cp cimport *
from . cimport everything as cqlat

cdef class Geometry:

    cdef cqlat.Geometry xx

    cdef readonly long cdata

cdef class FieldBase:

    pass

cdef class SelectedField:

    pass

cdef class SelectedPoints:

    pass
