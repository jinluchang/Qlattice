cdef class FieldTYPENAME(FieldBase):

    cdef cqlat.Field[cqlat.TYPENAME] xx

    cdef readonly long cdata

    cdef readonly int view_count

cdef class SelectedFieldTYPENAME(SelectedFieldBase):

    cdef cqlat.SelectedField[cqlat.TYPENAME] xx

    cdef readonly long cdata

    cdef public FieldSelection fsel

    cdef readonly int view_count

cdef class SelectedPointsTYPENAME(SelectedPointsBase):

    cdef cqlat.SelectedPoints[cqlat.TYPENAME] xx

    cdef readonly long cdata

    cdef public PointSelection psel

    cdef readonly int view_count
