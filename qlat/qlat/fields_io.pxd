from . cimport everything as cqlat

cdef class ShuffledFieldsWriter:

    cdef cqlat.ShuffledFieldsWriter xx

cdef class ShuffledFieldsReader:

    cdef cqlat.ShuffledFieldsReader xx
    cdef readonly set tags

cdef class ShuffledBitSet:

    cdef cqlat.ShuffledBitSet xx
