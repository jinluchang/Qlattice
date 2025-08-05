# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from . cimport everything as cc
from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

from .timer cimport *

### -------------------------------------------------------------------

cdef class Buffer:

    def __cinit__(self,
                  object obj=None,
                  int ndim=1,
                  int itemsize=1,
                  char* fmt=NULL,
                  char* buf=NULL,
                  ):
        self.obj = obj
        self.ndim = ndim
        self.shape_strides.resize(ndim * 2)
        self.format = fmt
        self.itemsize = itemsize
        self.buf = buf
        # need to set shape with self.set_dim_size(dim, dim_size)
        # need to call self.update_strides_from_shape()

    cdef void set_buffer(self, Py_buffer *buffer, int flags):
        cdef Py_ssize_t* shape = NULL
        cdef Py_ssize_t* strides = NULL
        if self.ndim != 0:
            shape = &self.shape_strides[0]
            strides = &self.shape_strides[self.ndim]
        buffer.buf = self.buf
        if flags & PyBUF_FORMAT:
            buffer.format = self.format
        else:
            buffer.format = NULL
        buffer.internal = NULL
        buffer.itemsize = self.itemsize
        buffer.len = self.get_len()
        buffer.ndim = self.ndim
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = shape
        buffer.strides = strides
        buffer.suboffsets = NULL

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        self.set_buffer(buffer, flags)

    def __releasebuffer__(self, Py_buffer *buffer):
        self.obj.release_buffer(self)

    cdef void set_dim_size(self, int dim, Py_ssize_t dim_size):
        self.shape_strides[dim] = dim_size

    cdef void update_strides_from_shape(self):
        cdef Py_ssize_t* shape = NULL
        cdef Py_ssize_t* strides = NULL
        if self.ndim != 0:
            shape = &self.shape_strides[0]
            strides = &self.shape_strides[self.ndim]
        cdef int i
        cdef Py_ssize_t stride = self.itemsize
        for i in range(self.ndim - 1, -1, -1):
            strides[i] = stride
            stride *= shape[i]

    cdef Py_ssize_t get_len(self):
        cdef int i
        cdef Py_ssize_t ret = self.itemsize
        for i in range(self.ndim):
            ret *= self.shape_strides[i]
        return ret

### -------------------------------------------------------------------

# ColorMatrix WilsonMatrix NonRelWilsonMatrix SpinMatrix WilsonVector ComplexD ComplexF RealD RealF Long Int Int64t Int32t Int8t Char

cdef class ElemType:
    name = ""

cdef class ElemTypeColorMatrix(ElemType):
    name = "ColorMatrix"
    sizeof_m = sizeof(cc.ColorMatrix)
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.ComplexD)
    @staticmethod
    cdef int ndim():
        return 2
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t](2, 3)
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.ColorMatrix)

cdef class ElemTypeWilsonMatrix(ElemType):
    name = "WilsonMatrix"
    sizeof_m = sizeof(cc.WilsonMatrix)
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.ComplexD)
    @staticmethod
    cdef int ndim():
        return 2
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t](2, 12)
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.WilsonMatrix)

cdef class ElemTypeNonRelWilsonMatrix(ElemType):
    name = "NonRelWilsonMatrix"
    sizeof_m = sizeof(cc.NonRelWilsonMatrix)
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.ComplexD)
    @staticmethod
    cdef int ndim():
        return 2
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t](2, 6)
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.NonRelWilsonMatrix)

cdef class ElemTypeIsospinMatrix(ElemType):
    name = "IsospinMatrix"
    sizeof_m = sizeof(cc.IsospinMatrix)
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.ComplexD)
    @staticmethod
    cdef int ndim():
        return 2
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t](2, 2)
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.IsospinMatrix)

cdef class ElemTypeSpinMatrix(ElemType):
    name = "SpinMatrix"
    sizeof_m = sizeof(cc.SpinMatrix)
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.ComplexD)
    @staticmethod
    cdef int ndim():
        return 2
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t](2, 4)
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.SpinMatrix)

cdef class ElemTypeWilsonVector(ElemType):
    name = "WilsonVector"
    sizeof_m = sizeof(cc.WilsonVector)
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.ComplexD)
    @staticmethod
    cdef int ndim():
        return 1
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t](1, 12)
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.WilsonVector)

cdef class ElemTypeComplexD(ElemType):
    name = "ComplexD"
    sizeof_m = sizeof(cc.ComplexD)
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.ComplexD)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.ComplexD)

cdef class ElemTypeComplexF(ElemType):
    name = "ComplexF"
    sizeof_m = sizeof(cc.ComplexF)
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zf'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.ComplexF)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.ComplexF)

cdef class ElemTypeRealD(ElemType):
    name = "RealD"
    sizeof_m = sizeof(cc.RealD)
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.RealD)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.RealD)

cdef class ElemTypeRealF(ElemType):
    name = "RealF"
    sizeof_m = sizeof(cc.RealF)
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'f'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.RealF)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.RealF)

cdef class ElemTypeLong(ElemType):
    name = "Long"
    sizeof_m = sizeof(cc.Long)
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'q'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Long)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.Long)

cdef class ElemTypeInt(ElemType):
    name = "Int"
    sizeof_m = sizeof(cc.Int)
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'i'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Int)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.Int)

cdef class ElemTypeChar(ElemType):
    name = "Char"
    sizeof_m = sizeof(cc.Char)
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'b'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Char)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.Char)

cdef class ElemTypeInt64t(ElemType):
    name = "Int64t"
    sizeof_m = sizeof(cc.Int64t)
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'q'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Int64t)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.Int64t)

cdef class ElemTypeInt32t(ElemType):
    name = "Int32t"
    sizeof_m = sizeof(cc.Int32t)
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'i'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Int32t)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.Int32t)

cdef class ElemTypeInt8t(ElemType):
    name = "Int8t"
    sizeof_m = sizeof(cc.Int8t)
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'b'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Int8t)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.Int8t)

### -------------------------------------------------------------------
