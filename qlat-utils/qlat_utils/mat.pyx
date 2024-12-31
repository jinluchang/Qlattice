# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from . cimport everything as cc
from .types cimport *

from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

import numpy as np

### -------------------------------------------------------------------

cdef class WilsonMatrix:

    def __cinit__(self):
        self.cdata = <cc.Long>&(self.xx)

    def __imatmul__(self, WilsonMatrix v1):
        self.xx = v1.xx
        return self

    def copy(self, cc.Bool is_copying_data = True):
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
        cdef int ndim = ElemTypeWilsonMatrix.ndim()
        cdef Buffer buf = Buffer(self, ndim)
        buf.format = ElemTypeWilsonMatrix.format()
        buf.itemsize = ElemTypeWilsonMatrix.itemsize()
        buf.buf = <char*>(self.xx.data())
        cdef cc.std_vector[Py_ssize_t] vec = ElemTypeWilsonMatrix.shape()
        cdef int dim
        for dim in range(ndim):
            buf.set_dim_size(dim, vec[dim])
        buf.update_strides_from_shape()
        buf.set_buffer(buffer, flags)

    def release_buffer(self, Buffer buf):
        pass

    def __setitem__(self, idx, val):
        """
        Implemented in terms of ``np.asarray``
        """
        np.asarray(self)[idx] = val

    def __getitem__(self, idx):
        """
        Implemented in terms of ``np.asarray``
        """
        return np.asarray(self)[idx]

    def g5_herm(self):
        self.xx = cc.g5_herm(self.xx)

    def qnorm(self):
        return cc.qnorm(self.xx)

    def __eq__(self, WilsonMatrix v1):
        return cc.qnorm(self.xx - v1.xx) == 0.0

    def __iadd__(self, WilsonMatrix v1):
        self.xx += v1.xx
        return self

    def __isub__(self, WilsonMatrix v1):
        self.xx -= v1.xx
        return self

    def __imul__(self, cc.PyComplexD v1):
        self.xx *= cc.ccpy_d(v1)
        return self

    def __add__(self, WilsonMatrix v1):
        cdef WilsonMatrix x = WilsonMatrix()
        x.xx = self.xx + v1.xx
        return x

    def __sub__(self, WilsonMatrix v1):
        cdef WilsonMatrix x = WilsonMatrix()
        x.xx = self.xx - v1.xx
        return x

    def __mul__(self, cc.PyComplexD v1):
        cdef WilsonMatrix x = WilsonMatrix()
        x.xx = cc.ccpy_d(v1) * self.xx
        return self

    def conjugate(self):
        cdef WilsonMatrix x = WilsonMatrix()
        x.xx = cc.matrix_conjugate(self.xx)
        return x

    def transpose(self):
        cdef WilsonMatrix x = WilsonMatrix()
        x.xx = cc.matrix_transpose(self.xx)
        return x

    def adjoint(self):
        cdef WilsonMatrix x = WilsonMatrix()
        x.xx = cc.matrix_adjoint(self.xx)
        return x

    def __repr__(self):
        return f"WilsonMatrix({self[:].tolist()!r})"

### -------------------------------------------------------------------

cdef class SpinMatrix:

    def __cinit__(self):
        self.cdata = <cc.Long>&(self.xx)

    def __imatmul__(self, SpinMatrix v1):
        self.xx = v1.xx
        return self

    def copy(self, cc.Bool is_copying_data = True):
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
        cdef int ndim = ElemTypeSpinMatrix.ndim()
        cdef Buffer buf = Buffer(self, ndim)
        buf.format = ElemTypeSpinMatrix.format()
        buf.itemsize = ElemTypeSpinMatrix.itemsize()
        buf.buf = <char*>(self.xx.data())
        cdef cc.std_vector[Py_ssize_t] vec = ElemTypeSpinMatrix.shape()
        cdef int dim
        for dim in range(ndim):
            buf.set_dim_size(dim, vec[dim])
        buf.update_strides_from_shape()
        buf.set_buffer(buffer, flags)

    def release_buffer(self, Buffer buf):
        pass

    def __setitem__(self, idx, val):
        """
        Implemented in terms of ``np.asarray``
        """
        np.asarray(self)[idx] = val

    def __getitem__(self, idx):
        """
        Implemented in terms of ``np.asarray``
        """
        return np.asarray(self)[idx]

    def qnorm(self):
        return cc.qnorm(self.xx)

    def __eq__(self, SpinMatrix v1):
        return cc.qnorm(self.xx - v1.xx) == 0.0

    def __iadd__(self, SpinMatrix v1):
        self.xx += v1.xx
        return self

    def __isub__(self, SpinMatrix v1):
        self.xx -= v1.xx
        return self

    def __imul__(self, cc.PyComplexD v1):
        self.xx *= cc.ccpy_d(v1)
        return self

    def __add__(self, SpinMatrix v1):
        cdef SpinMatrix x = SpinMatrix()
        x.xx = self.xx + v1.xx
        return x

    def __sub__(self, SpinMatrix v1):
        cdef SpinMatrix x = SpinMatrix()
        x.xx = self.xx - v1.xx
        return x

    def __mul__(self, cc.PyComplexD v1):
        cdef SpinMatrix x = SpinMatrix()
        x.xx = cc.ccpy_d(v1) * self.xx
        return self

    def conjugate(self):
        cdef SpinMatrix x = SpinMatrix()
        x.xx = cc.matrix_conjugate(self.xx)
        return x

    def transpose(self):
        cdef SpinMatrix x = SpinMatrix()
        x.xx = cc.matrix_transpose(self.xx)
        return x

    def adjoint(self):
        cdef SpinMatrix x = SpinMatrix()
        x.xx = cc.matrix_adjoint(self.xx)
        return x

    def __repr__(self):
        return f"SpinMatrix({self[:].tolist()!r})"

### -------------------------------------------------------------------

cdef class ColorMatrix:

    def __cinit__(self):
        self.cdata = <cc.Long>&(self.xx)

    def __imatmul__(self, ColorMatrix v1):
        self.xx = v1.xx
        return self

    def copy(self, cc.Bool is_copying_data = True):
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
        cdef int ndim = ElemTypeColorMatrix.ndim()
        cdef Buffer buf = Buffer(self, ndim)
        buf.format = ElemTypeColorMatrix.format()
        buf.itemsize = ElemTypeColorMatrix.itemsize()
        buf.buf = <char*>(self.xx.data())
        cdef cc.std_vector[Py_ssize_t] vec = ElemTypeColorMatrix.shape()
        cdef int dim
        for dim in range(ndim):
            buf.set_dim_size(dim, vec[dim])
        buf.update_strides_from_shape()
        buf.set_buffer(buffer, flags)

    def release_buffer(self, Buffer buf):
        pass

    def __setitem__(self, idx, val):
        """
        Implemented in terms of ``np.asarray``
        """
        np.asarray(self)[idx] = val

    def __getitem__(self, idx):
        """
        Implemented in terms of ``np.asarray``
        """
        return np.asarray(self)[idx]

    def qnorm(self):
        return cc.qnorm(self.xx)

    def __eq__(self, ColorMatrix v1):
        return cc.qnorm(self.xx - v1.xx) == 0.0

    def __iadd__(self, ColorMatrix v1):
        self.xx += v1.xx
        return self

    def __isub__(self, ColorMatrix v1):
        self.xx -= v1.xx
        return self

    def __imul__(self, cc.PyComplexD v1):
        self.xx *= cc.ccpy_d(v1)
        return self

    def __add__(self, ColorMatrix v1):
        cdef ColorMatrix x = ColorMatrix()
        x.xx = self.xx + v1.xx
        return x

    def __sub__(self, ColorMatrix v1):
        cdef ColorMatrix x = ColorMatrix()
        x.xx = self.xx - v1.xx
        return x

    def __mul__(self, cc.PyComplexD v1):
        cdef ColorMatrix x = ColorMatrix()
        x.xx = cc.ccpy_d(v1) * self.xx
        return self

    def conjugate(self):
        cdef ColorMatrix x = ColorMatrix()
        x.xx = cc.matrix_conjugate(self.xx)
        return x

    def transpose(self):
        cdef ColorMatrix x = ColorMatrix()
        x.xx = cc.matrix_transpose(self.xx)
        return x

    def adjoint(self):
        cdef ColorMatrix x = ColorMatrix()
        x.xx = cc.matrix_adjoint(self.xx)
        return x

    def __repr__(self):
        return f"ColorMatrix({self[:].tolist()!r})"

### -------------------------------------------------------------------

def mat_tr_sm(SpinMatrix v):
    return cc.pycc_d(cc.matrix_trace(v.xx))

def mat_tr_cm(ColorMatrix v):
    return cc.pycc_d(cc.matrix_trace(v.xx))

def mat_tr_wm(WilsonMatrix v):
    return cc.pycc_d(cc.matrix_trace(v.xx))

def mat_tr_wm_wm(WilsonMatrix v1, WilsonMatrix v2):
    return cc.pycc_d(cc.matrix_trace(v1.xx, v2.xx))

def mat_tr_wm_sm(WilsonMatrix v1, SpinMatrix v2):
    return cc.pycc_d(cc.matrix_trace(v1.xx, v2.xx))

def mat_tr_sm_wm(SpinMatrix v1, WilsonMatrix v2):
    return cc.pycc_d(cc.matrix_trace(v1.xx, v2.xx))

def mat_tr_sm_sm(SpinMatrix v1, SpinMatrix v2):
    return cc.pycc_d(cc.matrix_trace(v1.xx, v2.xx))

def mat_tr_wm_cm(WilsonMatrix v1, ColorMatrix v2):
    return cc.pycc_d(cc.matrix_trace(v1.xx, v2.xx))

def mat_tr_cm_wm(ColorMatrix v1, WilsonMatrix v2):
    return cc.pycc_d(cc.matrix_trace(v1.xx, v2.xx))

def mat_tr_cm_cm(ColorMatrix v1, ColorMatrix v2):
    return cc.pycc_d(cc.matrix_trace(v1.xx, v2.xx))

### -------------------------------------------------------------------

def mat_mul_a_wm(cc.PyComplexD v1, WilsonMatrix v2):
    cdef WilsonMatrix x = WilsonMatrix()
    x.xx = cc.ccpy_d(v1) * v2.xx
    return x

def mat_mul_a_sm(cc.PyComplexD v1, SpinMatrix v2):
    cdef SpinMatrix x = SpinMatrix()
    x.xx = cc.ccpy_d(v1) * v2.xx
    return x

def mat_mul_a_cm(cc.PyComplexD v1, ColorMatrix v2):
    cdef ColorMatrix x = ColorMatrix()
    x.xx = cc.ccpy_d(v1) * v2.xx
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

def mat_add_wm_wm(WilsonMatrix v1, WilsonMatrix v2):
    cdef WilsonMatrix x = WilsonMatrix()
    x.xx = v1.xx + v2.xx
    return x

def mat_add_sm_sm(SpinMatrix v1, SpinMatrix v2):
    cdef SpinMatrix x = SpinMatrix()
    x.xx = v1.xx + v2.xx
    return x

def mat_add_cm_cm(ColorMatrix v1, ColorMatrix v2):
    cdef ColorMatrix x = ColorMatrix()
    x.xx = v1.xx + v2.xx
    return x

### -------------------------------------------------------------------

def mat_epsilon_contraction_wm_wm_wm(int v_s1, int b_s1, int v_s2, int b_s2, int v_s3, int b_s3,
                                     WilsonMatrix wm1, WilsonMatrix wm2, WilsonMatrix wm3):
    cdef cc.ComplexD val
    val = cc.epsilon_contraction(v_s1, b_s1, v_s2, b_s2, v_s3, b_s3, wm1.xx, wm2.xx, wm3.xx)
    return cc.pycc_d(val)

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
