cdef class FieldTYPENAME(FieldBase):

    cdef cqlat.Field[cqlat.TYPENAME] xx

    cdef readonly long cdata

cdef class SelectedFieldTYPENAME(SelectedFieldBase):

    cdef cqlat.SelectedField[cqlat.TYPENAME] xx

    cdef readonly long cdata

cdef class SelectedPointsTYPENAME(SelectedPointsBase):

    cdef cqlat.SelectedPoints[cqlat.TYPENAME] xx

    cdef readonly long cdata
