cdef class ElemTypeTYPENAME(ElemType):

    name = "TYPENAME"

### -------------------------------------------------------------------

cdef class FieldTYPENAME(Field):

    ctype = ElemTypeTYPENAME

    def __cinit__(self):
        self.xx = cc.Field[cc.TYPENAME]()
        self.cdata = <long>&(self.xx)

### -------------------------------------------------------------------
