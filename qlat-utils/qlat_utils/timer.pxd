from . cimport everything as cqlat_utils

cdef class Timer:
    cdef cqlat_utils.Timer xx
    cdef cqlat_utils.bool is_verbose

cdef class TimerNone:
    pass
