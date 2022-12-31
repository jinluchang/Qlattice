from . cimport everything as cqlat_utils

cdef class Types:

    pass

cdef class Timer:

    cdef cqlat_utils.Timer xx

    cdef cqlat_utils.bool is_verbose

cdef class TimerNone:

    pass

cdef class Buffer:

    cdef object obj

    cdef int ndim

    cdef Py_ssize_t itemsize

    cdef cqlat_utils.std_vector[Py_ssize_t] shape_strides # shape.size() == 2 * ndim

    cdef Py_ssize_t get_len(self)

    cdef void set_strides(self)

cdef class Coordinate:

    cdef cqlat_utils.Coordinate xx

cdef class RngState:

    cdef cqlat_utils.RngState xx

    cdef readonly long cdata

cdef class WilsonMatrix:

    cdef cqlat_utils.WilsonMatrix xx

    cdef readonly long cdata

cdef class SpinMatrix:

    cdef cqlat_utils.SpinMatrix xx

    cdef readonly long cdata

cdef cqlat_utils.SpinMatrix gamma_matrix_0
cdef cqlat_utils.SpinMatrix gamma_matrix_1
cdef cqlat_utils.SpinMatrix gamma_matrix_2
cdef cqlat_utils.SpinMatrix gamma_matrix_3
cdef cqlat_utils.SpinMatrix gamma_matrix_5
