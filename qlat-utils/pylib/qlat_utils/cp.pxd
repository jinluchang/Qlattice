from . cimport everything as cqlat_utils

cdef class ElemType:

    pass

cdef class ElemTypeColorMatrix(ElemType):
    @staticmethod
    cdef char* buffer_format()
    @staticmethod
    cdef Py_ssize_t buffer_itemsize()
    @staticmethod
    cdef int buffer_ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] buffer_shape()
cdef class ElemTypeWilsonMatrix(ElemType):
    @staticmethod
    cdef char* buffer_format()
    @staticmethod
    cdef Py_ssize_t buffer_itemsize()
    @staticmethod
    cdef int buffer_ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] buffer_shape()
cdef class ElemTypeNonRelWilsonMatrix(ElemType):
    @staticmethod
    cdef char* buffer_format()
    @staticmethod
    cdef Py_ssize_t buffer_itemsize()
    @staticmethod
    cdef int buffer_ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] buffer_shape()
cdef class ElemTypeSpinMatrix(ElemType):
    @staticmethod
    cdef char* buffer_format()
    @staticmethod
    cdef Py_ssize_t buffer_itemsize()
    @staticmethod
    cdef int buffer_ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] buffer_shape()
cdef class ElemTypeWilsonVector(ElemType):
    @staticmethod
    cdef char* buffer_format()
    @staticmethod
    cdef Py_ssize_t buffer_itemsize()
    @staticmethod
    cdef int buffer_ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] buffer_shape()
cdef class ElemTypeComplex(ElemType):
    @staticmethod
    cdef char* buffer_format()
    @staticmethod
    cdef Py_ssize_t buffer_itemsize()
    @staticmethod
    cdef int buffer_ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] buffer_shape()
cdef class ElemTypeComplexF(ElemType):
    @staticmethod
    cdef char* buffer_format()
    @staticmethod
    cdef Py_ssize_t buffer_itemsize()
    @staticmethod
    cdef int buffer_ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] buffer_shape()
cdef class ElemTypeDouble(ElemType):
    @staticmethod
    cdef char* buffer_format()
    @staticmethod
    cdef Py_ssize_t buffer_itemsize()
    @staticmethod
    cdef int buffer_ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] buffer_shape()
cdef class ElemTypeFloat(ElemType):
    @staticmethod
    cdef char* buffer_format()
    @staticmethod
    cdef Py_ssize_t buffer_itemsize()
    @staticmethod
    cdef int buffer_ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] buffer_shape()
cdef class ElemTypeLong(ElemType):
    @staticmethod
    cdef char* buffer_format()
    @staticmethod
    cdef Py_ssize_t buffer_itemsize()
    @staticmethod
    cdef int buffer_ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] buffer_shape()
cdef class ElemTypeInt64t(ElemType):
    @staticmethod
    cdef char* buffer_format()
    @staticmethod
    cdef Py_ssize_t buffer_itemsize()
    @staticmethod
    cdef int buffer_ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] buffer_shape()
cdef class ElemTypeInt8t(ElemType):
    @staticmethod
    cdef char* buffer_format()
    @staticmethod
    cdef Py_ssize_t buffer_itemsize()
    @staticmethod
    cdef int buffer_ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] buffer_shape()
cdef class ElemTypeChar(ElemType):
    @staticmethod
    cdef char* buffer_format()
    @staticmethod
    cdef Py_ssize_t buffer_itemsize()
    @staticmethod
    cdef int buffer_ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] buffer_shape()

### -------------------------------------------------------------------

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
