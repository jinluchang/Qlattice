from . cimport everything as cqlat

cdef class CurrentMoments:

    cdef cqlat.CurrentMoments[cqlat.WilsonMatrix] xx
