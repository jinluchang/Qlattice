cdef class ElemTypeTYPENAME(ElemType):

    pass

cdef class FieldTYPENAME(Field):

    cdef cqlat.Field[cqlat.TYPENAME] xx

    cdef readonly long cdata
