cdef class FieldTYPENAME(FieldBase):

    cdef cqlat.Field[cqlat.TYPENAME] xx

cdef class SelectedFieldTYPENAME(SelectedFieldBase):

    cdef cqlat.SelectedField[cqlat.TYPENAME] xx

cdef class SelectedPointsTYPENAME(SelectedPointsBase):

    cdef cqlat.SelectedPoints[cqlat.TYPENAME] xx
