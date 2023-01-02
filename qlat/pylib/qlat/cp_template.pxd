from qlat_utils.cp cimport *
from . cimport everything as cqlat

cdef class Geometry:

    cdef cqlat.Geometry xx

    cdef readonly long cdata

cdef class PointSelection:

    cdef cqlat.PointSelection xx

    cdef readonly long cdata

cdef class FieldSelection:

    cdef cqlat.FieldSelection xx

    cdef readonly long cdata

cdef class FieldBase:

    pass

cdef class SelectedFieldBase:

    pass

cdef class SelectedPointsBase:

    pass
