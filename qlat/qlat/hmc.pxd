from . cimport everything as cqlat
from .field_types cimport FieldColorMatrix

cdef class GaugeMomentum(FieldColorMatrix):

    cdef cqlat.Handle[cqlat.GaugeMomentum] xxx(self)
