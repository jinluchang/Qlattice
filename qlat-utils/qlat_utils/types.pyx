# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from . cimport everything as cc
from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

from .timer cimport *

### -------------------------------------------------------------------

cdef class Buffer:

    def __cinit__(self, object obj = None, int ndim = 1, int itemsize = 1):
        self.obj = obj
        self.ndim = ndim
        self.itemsize = itemsize
        self.shape_strides.resize(ndim * 2)

    cdef Py_ssize_t get_len(self):
        cdef int i
        cdef Py_ssize_t ret = 1
        for i in range(self.ndim):
            ret *= self.shape_strides[i]
        return ret

    cdef void set_strides(self):
        cdef Py_ssize_t* shapes = &self.shape_strides[0]
        cdef Py_ssize_t* strides = &self.shape_strides[self.ndim]
        cdef int i
        cdef Py_ssize_t stride = self.itemsize
        for i in range(self.ndim - 1, -1, -1):
            strides[i] = stride
            stride *= shapes[i]

### -------------------------------------------------------------------

# ColorMatrix WilsonMatrix NonRelWilsonMatrix SpinMatrix WilsonVector Complex ComplexF Double Float Long Int64t Int8t Char

cdef class ElemType:
    name = ""

cdef class ElemTypeColorMatrix(ElemType):
    name = "ColorMatrix"
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Complex)
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
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Complex)
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
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Complex)
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
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Complex)
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
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Complex)
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
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Complex)
    @staticmethod
    cdef int ndim():
        return 1
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t](1, 12)
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.WilsonVector)

cdef class ElemTypeComplex(ElemType):
    name = "Complex"
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'Zd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Complex)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.Complex)

cdef class ElemTypeComplexF(ElemType):
    name = "ComplexF"
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

cdef class ElemTypeDouble(ElemType):
    name = "Double"
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'd'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Double)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.Double)

cdef class ElemTypeFloat(ElemType):
    name = "Float"
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'f'
        return fmt
    @staticmethod
    cdef Py_ssize_t itemsize():
        return sizeof(cc.Float)
    @staticmethod
    cdef int ndim():
        return 0
    @staticmethod
    cdef cc.std_vector[Py_ssize_t] shape():
        return cc.std_vector[Py_ssize_t]()
    @staticmethod
    cdef Py_ssize_t size():
        return sizeof(cc.Float)

cdef class ElemTypeLong(ElemType):
    name = "Long"
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'q'
        if not sizeof(cc.Long) == 8:
            assert sizeof(cc.Long) == 4
            fmt = 'l'
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

cdef class ElemTypeInt64t(ElemType):
    name = "Int64t"
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

cdef class ElemTypeInt8t(ElemType):
    name = "Int8t"
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

cdef class ElemTypeChar(ElemType):
    name = "Char"
    @staticmethod
    cdef char* format():
        cdef char* fmt = 'c'
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

### -------------------------------------------------------------------

cdef class WilsonMatrix:

    def __cinit__(self):
        self.cdata = <long>&(self.xx)

    def __imatmul__(self, WilsonMatrix v1):
        self.xx = v1.xx
        return self

    def copy(self, cc.bool is_copying_data = True):
        cdef WilsonMatrix x = WilsonMatrix()
        if is_copying_data:
            x.xx = self.xx
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def set_zero(self):
        cc.set_zero(self.xx)

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef Buffer buf = Buffer(self, ElemTypeWilsonMatrix.ndim(), ElemTypeWilsonMatrix.itemsize())
        cdef char* fmt = ElemTypeWilsonMatrix.format()
        cdef cc.std_vector[Py_ssize_t] vec = ElemTypeWilsonMatrix.shape()
        assert vec.size() == buf.ndim
        cdef Py_ssize_t* shape = &buf.shape_strides[0]
        cdef Py_ssize_t* strides = &buf.shape_strides[buf.ndim]
        cdef int i
        for i in range(buf.ndim):
            shape[i] = vec[i]
        buf.set_strides()
        buffer.buf = <char*>(self.xx.data())
        if flags & PyBUF_FORMAT:
            buffer.format = fmt
        else:
            buffer.format = NULL
        buffer.internal = NULL
        buffer.itemsize = buf.itemsize
        buffer.len = buf.get_len()
        buffer.ndim = buf.ndim
        buffer.obj = buf
        buffer.readonly = 0
        buffer.shape = shape
        buffer.strides = strides
        buffer.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buffer):
        pass

    def g5_herm(self):
        self.xx = cc.g5_herm(self.xx)

### -------------------------------------------------------------------

cdef class SpinMatrix:

    def __cinit__(self):
        self.cdata = <long>&(self.xx)

    def __imatmul__(self, SpinMatrix v1):
        self.xx = v1.xx
        return self

    def copy(self, cc.bool is_copying_data = True):
        cdef SpinMatrix x = SpinMatrix()
        if is_copying_data:
            x.xx = self.xx
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def set_zero(self):
        cc.set_zero(self.xx)

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef Buffer buf = Buffer(self, ElemTypeSpinMatrix.ndim(), ElemTypeSpinMatrix.itemsize())
        cdef char* fmt = ElemTypeWilsonMatrix.format()
        cdef cc.std_vector[Py_ssize_t] vec = ElemTypeSpinMatrix.shape()
        assert vec.size() == buf.ndim
        cdef Py_ssize_t* shape = &buf.shape_strides[0]
        cdef Py_ssize_t* strides = &buf.shape_strides[buf.ndim]
        cdef int i
        for i in range(buf.ndim):
            shape[i] = vec[i]
        buf.set_strides()
        buffer.buf = <char*>(self.xx.data())
        if flags & PyBUF_FORMAT:
            buffer.format = fmt
        else:
            buffer.format = NULL
        buffer.internal = NULL
        buffer.itemsize = buf.itemsize
        buffer.len = buf.get_len()
        buffer.ndim = buf.ndim
        buffer.obj = buf
        buffer.readonly = 0
        buffer.shape = shape
        buffer.strides = strides
        buffer.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buffer):
        pass

### -------------------------------------------------------------------

cdef class ColorMatrix:

    def __cinit__(self):
        self.cdata = <long>&(self.xx)

    def __imatmul__(self, ColorMatrix v1):
        self.xx = v1.xx
        return self

    def copy(self, cc.bool is_copying_data = True):
        cdef ColorMatrix x = ColorMatrix()
        if is_copying_data:
            x.xx = self.xx
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def set_zero(self):
        cc.set_zero(self.xx)

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef Buffer buf = Buffer(self, ElemTypeColorMatrix.ndim(), ElemTypeColorMatrix.itemsize())
        cdef char* fmt = ElemTypeWilsonMatrix.format()
        cdef cc.std_vector[Py_ssize_t] vec = ElemTypeColorMatrix.shape()
        assert vec.size() == buf.ndim
        cdef Py_ssize_t* shape = &buf.shape_strides[0]
        cdef Py_ssize_t* strides = &buf.shape_strides[buf.ndim]
        cdef int i
        for i in range(buf.ndim):
            shape[i] = vec[i]
        buf.set_strides()
        buffer.buf = <char*>(self.xx.data())
        if flags & PyBUF_FORMAT:
            buffer.format = fmt
        else:
            buffer.format = NULL
        buffer.internal = NULL
        buffer.itemsize = buf.itemsize
        buffer.len = buf.get_len()
        buffer.ndim = buf.ndim
        buffer.obj = buf
        buffer.readonly = 0
        buffer.shape = shape
        buffer.strides = strides
        buffer.suboffsets = NULL

    def __releasebuffer__(self, Py_buffer *buffer):
        pass

### -------------------------------------------------------------------

def mat_tr_sm(SpinMatrix v):
    return cc.matrix_trace(v.xx)

def mat_tr_cm(ColorMatrix v):
    return cc.matrix_trace(v.xx)

def mat_tr_wm(WilsonMatrix v):
    return cc.matrix_trace(v.xx)

def mat_tr_wm_wm(WilsonMatrix v1, WilsonMatrix v2):
    return cc.matrix_trace(v1.xx, v2.xx)

def mat_tr_wm_sm(WilsonMatrix v1, SpinMatrix v2):
    return cc.matrix_trace(v1.xx, v2.xx)

def mat_tr_sm_wm(SpinMatrix v1, WilsonMatrix v2):
    return cc.matrix_trace(v1.xx, v2.xx)

def mat_tr_sm_sm(SpinMatrix v1, SpinMatrix v2):
    return cc.matrix_trace(v1.xx, v2.xx)

def mat_tr_wm_cm(WilsonMatrix v1, ColorMatrix v2):
    return cc.matrix_trace(v1.xx, v2.xx)

def mat_tr_cm_wm(ColorMatrix v1, WilsonMatrix v2):
    return cc.matrix_trace(v1.xx, v2.xx)

def mat_tr_cm_cm(ColorMatrix v1, ColorMatrix v2):
    return cc.matrix_trace(v1.xx, v2.xx)

### -------------------------------------------------------------------

def mat_mul_a_wm(cc.Complex v1, WilsonMatrix v2):
    cdef WilsonMatrix x = WilsonMatrix()
    x.xx = v1 * v2.xx
    return x

def mat_mul_a_sm(cc.Complex v1, SpinMatrix v2):
    cdef SpinMatrix x = SpinMatrix()
    x.xx = v1 * v2.xx
    return x

def mat_mul_a_cm(cc.Complex v1, ColorMatrix v2):
    cdef ColorMatrix x = ColorMatrix()
    x.xx = v1 * v2.xx
    return x

def mat_mul_wm_wm(WilsonMatrix v1, WilsonMatrix v2):
    cdef WilsonMatrix x = WilsonMatrix()
    x.xx = v1.xx * v2.xx
    return x

def mat_mul_sm_wm(SpinMatrix v1, WilsonMatrix v2):
    cdef WilsonMatrix x = WilsonMatrix()
    x.xx = v1.xx * v2.xx
    return x

def mat_mul_wm_sm(WilsonMatrix v1, SpinMatrix v2):
    cdef WilsonMatrix x = WilsonMatrix()
    x.xx = v1.xx * v2.xx
    return x

def mat_mul_sm_sm(SpinMatrix v1, SpinMatrix v2):
    cdef SpinMatrix x = SpinMatrix()
    x.xx = v1.xx * v2.xx
    return x

def mat_mul_cm_wm(ColorMatrix v1, WilsonMatrix v2):
    cdef WilsonMatrix x = WilsonMatrix()
    x.xx = v1.xx * v2.xx
    return x

def mat_mul_wm_cm(WilsonMatrix v1, ColorMatrix v2):
    cdef WilsonMatrix x = WilsonMatrix()
    x.xx = v1.xx * v2.xx
    return x

def mat_mul_cm_cm(ColorMatrix v1, ColorMatrix v2):
    cdef ColorMatrix x = ColorMatrix()
    x.xx = v1.xx * v2.xx
    return x

### -------------------------------------------------------------------

def get_gamma_matrix(int mu):
    cdef SpinMatrix x = SpinMatrix()
    x.xx = cc.get_gamma_matrix(mu)
    return x

gamma_matrix_0 = cc.get_gamma_matrix(0)
gamma_matrix_1 = cc.get_gamma_matrix(1)
gamma_matrix_2 = cc.get_gamma_matrix(2)
gamma_matrix_3 = cc.get_gamma_matrix(3)
gamma_matrix_5 = cc.get_gamma_matrix(5)

### -------------------------------------------------------------------

def benchmark_matrix_functions(int count):
    cc.benchmark_matrix_functions(count)

### -------------------------------------------------------------------

def as_wilson_matrix(x):
    cdef WilsonMatrix wm
    if isinstance(x, WilsonMatrix):
        return x
    elif x == 0:
        wm = WilsonMatrix()
        cc.set_zero(wm.xx)
        return wm

def wilson_matrix_g5_herm(WilsonMatrix x):
    cdef WilsonMatrix wm = WilsonMatrix()
    wm.xx = cc.g5_herm(x.xx)
    return wm

def as_wilson_matrix_g5_herm(x):
    cdef WilsonMatrix wm = WilsonMatrix()
    if isinstance(x, WilsonMatrix):
        wm.xx = cc.g5_herm((<WilsonMatrix>x).xx)
    elif x == 0:
        cc.set_zero(wm.xx)
    return wm
