from . cimport everything as cqlat

cdef class CommPlan:

    cdef cqlat.CommPlan xx

    cdef readonly cqlat.Long cdata
