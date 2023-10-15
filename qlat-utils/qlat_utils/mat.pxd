from . cimport everything as cqlat_utils

cdef class WilsonMatrix:
    cdef cqlat_utils.WilsonMatrix xx
    cdef readonly long cdata

cdef class SpinMatrix:
    cdef cqlat_utils.SpinMatrix xx
    cdef readonly long cdata

cdef class ColorMatrix:
    cdef cqlat_utils.ColorMatrix xx
    cdef readonly long cdata

cdef cqlat_utils.SpinMatrix gamma_matrix_0
cdef cqlat_utils.SpinMatrix gamma_matrix_1
cdef cqlat_utils.SpinMatrix gamma_matrix_2
cdef cqlat_utils.SpinMatrix gamma_matrix_3
cdef cqlat_utils.SpinMatrix gamma_matrix_5
