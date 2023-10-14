cdef class ElemType:

    name = ""

### -------------------------------------------------------------------

# ColorMatrix WilsonMatrix NonRelWilsonMatrix SpinMatrix WilsonVector Complex ComplexF Double Float Long Int64t Int8t Char

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
