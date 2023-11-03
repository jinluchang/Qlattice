from . cimport everything as cqlat

cdef class ShuffledFieldsWriter:

    cdef cqlat.ShuffledFieldsWriter xx
    cdef readonly cqlat.Long cdata

cdef class ShuffledFieldsReader:

    cdef cqlat.ShuffledFieldsReader xx
    cdef readonly set tags
    cdef readonly cqlat.Long cdata

cdef class ShuffledBitSet:

    cdef cqlat.ShuffledBitSet xx
    cdef readonly cqlat.Long cdata
