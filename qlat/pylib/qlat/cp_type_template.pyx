cdef class ElemTypeTYPENAME(ElemType):

    name = "TYPENAME"

### -------------------------------------------------------------------

cdef class FieldTYPENAME(Field):

    ctype = ElemTypeTYPENAME
