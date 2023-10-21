from . cimport everything as cqlat
from .field_types cimport FieldColorMatrix

cdef class GaugeField(FieldColorMatrix):

    cdef cqlat.Handle[cqlat.GaugeField] xxx(self)

cdef class GaugeTransform(FieldColorMatrix):

    cdef cqlat.Handle[cqlat.GaugeTransform] xxx(self)
