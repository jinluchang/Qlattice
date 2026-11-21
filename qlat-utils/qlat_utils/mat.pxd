from . cimport everything as cqlat_utils

# qlat::WilsonMatrix and qlat::SpinMatrix are over-aligned request types
# (up to 256 bytes), but CPython only guarantees 16-byte alignment for the
# memory of an object.  They are therefore kept in an over-aligned in-object
# storage and `xx` points into it (instead of holding the value inline, which
# used to rely on the unguaranteed object alignment).  qlat::ColorMatrix is
# only 16-byte aligned, so it stays a plain inline value.

cdef class WilsonMatrix:
    cdef cqlat_utils.QlatAlignedStorage[cqlat_utils.WilsonMatrix] storage
    cdef cqlat_utils.WilsonMatrix* xx

cdef class SpinMatrix:
    cdef cqlat_utils.QlatAlignedStorage[cqlat_utils.SpinMatrix] storage
    cdef cqlat_utils.SpinMatrix* xx

cdef class ColorMatrix:
    cdef cqlat_utils.ColorMatrix xx

cdef cqlat_utils.SpinMatrix gamma_matrix_0
cdef cqlat_utils.SpinMatrix gamma_matrix_1
cdef cqlat_utils.SpinMatrix gamma_matrix_2
cdef cqlat_utils.SpinMatrix gamma_matrix_3
cdef cqlat_utils.SpinMatrix gamma_matrix_5
