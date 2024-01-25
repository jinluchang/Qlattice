from . cimport everything as cqlat

cdef class GaugeAction:

    cdef cqlat.GaugeAction xx

    cdef readonly cqlat.Long cdata
