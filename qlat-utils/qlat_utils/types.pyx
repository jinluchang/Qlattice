# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat_utils.types``
=================================\n
Element type descriptors and buffer protocol support for lattice QCD data types.\n
Documentation: ``docs/qlat-utils/qlat_types.md``\n
.. note:: Update the documentation when updating this source file.
"""

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
    @staticmethod
    def py_format():
        return ElemTypeColorMatrix.format()
    @staticmethod
    def py_itemsize():
        return ElemTypeColorMatrix.itemsize()
    @staticmethod
    def py_ndim():
        return ElemTypeColorMatrix.ndim()
    @staticmethod
    def py_shape():
        return tuple(ElemTypeColorMatrix.shape())
    @staticmethod
    def py_size():
        return ElemTypeColorMatrix.size()

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
    @staticmethod
    def py_format():
        return ElemTypeWilsonMatrix.format()
    @staticmethod
    def py_itemsize():
        return ElemTypeWilsonMatrix.itemsize()
    @staticmethod
    def py_ndim():
        return ElemTypeWilsonMatrix.ndim()
    @staticmethod
    def py_shape():
        return tuple(ElemTypeWilsonMatrix.shape())
    @staticmethod
    def py_size():
        return ElemTypeWilsonMatrix.size()

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
    @staticmethod
    def py_format():
        return ElemTypeNonRelWilsonMatrix.format()
    @staticmethod
    def py_itemsize():
        return ElemTypeNonRelWilsonMatrix.itemsize()
    @staticmethod
    def py_ndim():
        return ElemTypeNonRelWilsonMatrix.ndim()
    @staticmethod
    def py_shape():
        return tuple(ElemTypeNonRelWilsonMatrix.shape())
    @staticmethod
    def py_size():
        return ElemTypeNonRelWilsonMatrix.size()

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
    @staticmethod
    def py_format():
        return ElemTypeIsospinMatrix.format()
    @staticmethod
    def py_itemsize():
        return ElemTypeIsospinMatrix.itemsize()
    @staticmethod
    def py_ndim():
        return ElemTypeIsospinMatrix.ndim()
    @staticmethod
    def py_shape():
        return tuple(ElemTypeIsospinMatrix.shape())
    @staticmethod
    def py_size():
        return ElemTypeIsospinMatrix.size()

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
    @staticmethod
    def py_format():
        return ElemTypeSpinMatrix.format()
    @staticmethod
    def py_itemsize():
        return ElemTypeSpinMatrix.itemsize()
    @staticmethod
    def py_ndim():
        return ElemTypeSpinMatrix.ndim()
    @staticmethod
    def py_shape():
        return tuple(ElemTypeSpinMatrix.shape())
    @staticmethod
    def py_size():
        return ElemTypeSpinMatrix.size()

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
    @staticmethod
    def py_format():
        return ElemTypeWilsonVector.format()
    @staticmethod
    def py_itemsize():
        return ElemTypeWilsonVector.itemsize()
    @staticmethod
    def py_ndim():
        return ElemTypeWilsonVector.ndim()
    @staticmethod
    def py_shape():
        return tuple(ElemTypeWilsonVector.shape())
    @staticmethod
    def py_size():
        return ElemTypeWilsonVector.size()

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
    @staticmethod
    def py_format():
        return ElemTypeComplexD.format()
    @staticmethod
    def py_itemsize():
        return ElemTypeComplexD.itemsize()
    @staticmethod
    def py_ndim():
        return ElemTypeComplexD.ndim()
    @staticmethod
    def py_shape():
        return tuple(ElemTypeComplexD.shape())
    @staticmethod
    def py_size():
        return ElemTypeComplexD.size()

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
    @staticmethod
    def py_format():
        return ElemTypeComplexF.format()
    @staticmethod
    def py_itemsize():
        return ElemTypeComplexF.itemsize()
    @staticmethod
    def py_ndim():
        return ElemTypeComplexF.ndim()
    @staticmethod
    def py_shape():
        return tuple(ElemTypeComplexF.shape())
    @staticmethod
    def py_size():
        return ElemTypeComplexF.size()

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
    @staticmethod
    def py_format():
        return ElemTypeRealD.format()
    @staticmethod
    def py_itemsize():
        return ElemTypeRealD.itemsize()
    @staticmethod
    def py_ndim():
        return ElemTypeRealD.ndim()
    @staticmethod
    def py_shape():
        return tuple(ElemTypeRealD.shape())
    @staticmethod
    def py_size():
        return ElemTypeRealD.size()

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
    @staticmethod
    def py_format():
        return ElemTypeRealF.format()
    @staticmethod
    def py_itemsize():
        return ElemTypeRealF.itemsize()
    @staticmethod
    def py_ndim():
        return ElemTypeRealF.ndim()
    @staticmethod
    def py_shape():
        return tuple(ElemTypeRealF.shape())
    @staticmethod
    def py_size():
        return ElemTypeRealF.size()

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
    @staticmethod
    def py_format():
        return ElemTypeLong.format()
    @staticmethod
    def py_itemsize():
        return ElemTypeLong.itemsize()
    @staticmethod
    def py_ndim():
        return ElemTypeLong.ndim()
    @staticmethod
    def py_shape():
        return tuple(ElemTypeLong.shape())
    @staticmethod
    def py_size():
        return ElemTypeLong.size()

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
    @staticmethod
    def py_format():
        return ElemTypeInt.format()
    @staticmethod
    def py_itemsize():
        return ElemTypeInt.itemsize()
    @staticmethod
    def py_ndim():
        return ElemTypeInt.ndim()
    @staticmethod
    def py_shape():
        return tuple(ElemTypeInt.shape())
    @staticmethod
    def py_size():
        return ElemTypeInt.size()

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
    @staticmethod
    def py_format():
        return ElemTypeChar.format()
    @staticmethod
    def py_itemsize():
        return ElemTypeChar.itemsize()
    @staticmethod
    def py_ndim():
        return ElemTypeChar.ndim()
    @staticmethod
    def py_shape():
        return tuple(ElemTypeChar.shape())
    @staticmethod
    def py_size():
        return ElemTypeChar.size()

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
    @staticmethod
    def py_format():
        return ElemTypeInt64t.format()
    @staticmethod
    def py_itemsize():
        return ElemTypeInt64t.itemsize()
    @staticmethod
    def py_ndim():
        return ElemTypeInt64t.ndim()
    @staticmethod
    def py_shape():
        return tuple(ElemTypeInt64t.shape())
    @staticmethod
    def py_size():
        return ElemTypeInt64t.size()

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
    @staticmethod
    def py_format():
        return ElemTypeInt32t.format()
    @staticmethod
    def py_itemsize():
        return ElemTypeInt32t.itemsize()
    @staticmethod
    def py_ndim():
        return ElemTypeInt32t.ndim()
    @staticmethod
    def py_shape():
        return tuple(ElemTypeInt32t.shape())
    @staticmethod
    def py_size():
        return ElemTypeInt32t.size()

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
    @staticmethod
    def py_format():
        return ElemTypeInt8t.format()
    @staticmethod
    def py_itemsize():
        return ElemTypeInt8t.itemsize()
    @staticmethod
    def py_ndim():
        return ElemTypeInt8t.ndim()
    @staticmethod
    def py_shape():
        return tuple(ElemTypeInt8t.shape())
    @staticmethod
    def py_size():
        return ElemTypeInt8t.size()

### -------------------------------------------------------------------
