from . cimport everything as cqlat_utils

cdef class Buffer:
    cdef object obj
    cdef char* buf # data pointer
    cdef int ndim
    cdef char* format
    cdef Py_ssize_t itemsize
    cdef cqlat_utils.std_vector[Py_ssize_t] shape_strides # shape.size() == 2 * ndim
    cdef Py_ssize_t get_len(self)
    cdef void set_dim_size(self, int dim, Py_ssize_t size)
    cdef void update_strides_from_shape(self)
    cdef void set_buffer(self, Py_buffer *buffer, int flags)

### -------------------------------------------------------------------

cdef class ElemType:
    pass

cdef class ElemTypeColorMatrix(ElemType):
    @staticmethod
    cdef char* format()
    @staticmethod
    cdef Py_ssize_t itemsize()
    @staticmethod
    cdef int ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] shape()
    @staticmethod
    cdef Py_ssize_t size()

cdef class ElemTypeWilsonMatrix(ElemType):
    @staticmethod
    cdef char* format()
    @staticmethod
    cdef Py_ssize_t itemsize()
    @staticmethod
    cdef int ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] shape()
    @staticmethod
    cdef Py_ssize_t size()

cdef class ElemTypeNonRelWilsonMatrix(ElemType):
    @staticmethod
    cdef char* format()
    @staticmethod
    cdef Py_ssize_t itemsize()
    @staticmethod
    cdef int ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] shape()
    @staticmethod
    cdef Py_ssize_t size()

cdef class ElemTypeIsospinMatrix(ElemType):
    @staticmethod
    cdef char* format()
    @staticmethod
    cdef Py_ssize_t itemsize()
    @staticmethod
    cdef int ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] shape()
    @staticmethod
    cdef Py_ssize_t size()

cdef class ElemTypeSpinMatrix(ElemType):
    @staticmethod
    cdef char* format()
    @staticmethod
    cdef Py_ssize_t itemsize()
    @staticmethod
    cdef int ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] shape()
    @staticmethod
    cdef Py_ssize_t size()

cdef class ElemTypeWilsonVector(ElemType):
    @staticmethod
    cdef char* format()
    @staticmethod
    cdef Py_ssize_t itemsize()
    @staticmethod
    cdef int ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] shape()
    @staticmethod
    cdef Py_ssize_t size()

cdef class ElemTypeComplexD(ElemType):
    @staticmethod
    cdef char* format()
    @staticmethod
    cdef Py_ssize_t itemsize()
    @staticmethod
    cdef int ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] shape()
    @staticmethod
    cdef Py_ssize_t size()

cdef class ElemTypeComplexF(ElemType):
    @staticmethod
    cdef char* format()
    @staticmethod
    cdef Py_ssize_t itemsize()
    @staticmethod
    cdef int ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] shape()
    @staticmethod
    cdef Py_ssize_t size()

cdef class ElemTypeRealD(ElemType):
    @staticmethod
    cdef char* format()
    @staticmethod
    cdef Py_ssize_t itemsize()
    @staticmethod
    cdef int ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] shape()
    @staticmethod
    cdef Py_ssize_t size()

cdef class ElemTypeRealF(ElemType):
    @staticmethod
    cdef char* format()
    @staticmethod
    cdef Py_ssize_t itemsize()
    @staticmethod
    cdef int ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] shape()
    @staticmethod
    cdef Py_ssize_t size()

cdef class ElemTypeLong(ElemType):
    @staticmethod
    cdef char* format()
    @staticmethod
    cdef Py_ssize_t itemsize()
    @staticmethod
    cdef int ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] shape()
    @staticmethod
    cdef Py_ssize_t size()

cdef class ElemTypeInt(ElemType):
    @staticmethod
    cdef char* format()
    @staticmethod
    cdef Py_ssize_t itemsize()
    @staticmethod
    cdef int ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] shape()
    @staticmethod
    cdef Py_ssize_t size()

cdef class ElemTypeInt64t(ElemType):
    @staticmethod
    cdef char* format()
    @staticmethod
    cdef Py_ssize_t itemsize()
    @staticmethod
    cdef int ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] shape()
    @staticmethod
    cdef Py_ssize_t size()

cdef class ElemTypeInt32t(ElemType):
    @staticmethod
    cdef char* format()
    @staticmethod
    cdef Py_ssize_t itemsize()
    @staticmethod
    cdef int ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] shape()
    @staticmethod
    cdef Py_ssize_t size()

cdef class ElemTypeInt8t(ElemType):
    @staticmethod
    cdef char* format()
    @staticmethod
    cdef Py_ssize_t itemsize()
    @staticmethod
    cdef int ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] shape()
    @staticmethod
    cdef Py_ssize_t size()

cdef class ElemTypeChar(ElemType):
    @staticmethod
    cdef char* format()
    @staticmethod
    cdef Py_ssize_t itemsize()
    @staticmethod
    cdef int ndim()
    @staticmethod
    cdef cqlat_utils.std_vector[Py_ssize_t] shape()
    @staticmethod
    cdef Py_ssize_t size()
