# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from . cimport everything as cc

from .mat cimport *
from .coordinate cimport *
from .lat_data cimport *

from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT
cimport numpy

from .timer import timer

import numpy as np
import functools

cdef class RngState:

    def __cinit__(self):
        self.cdata = <cc.Long>&(self.xx)

    def __init__(self, x=None, y=None):
        cdef cc.std_string seed
        if x is None:
            assert y is None
            # make a new rng
            self.xx = cc.RngState()
        elif isinstance(x, RngState):
            if y is None:
                # make a copy of x
                self.xx = (<RngState>x).xx
            else:
                # split into a new rng
                seed = str(y)
                self.xx = cc.RngState((<RngState>x).xx, seed)
        else:
            assert y is None
            # seed a new rng
            seed = str(x)
            self.xx = cc.RngState(seed)

    def __imatmul__(self, RngState v1):
        self.xx = v1.xx
        return self

    def copy(self, cc.Bool is_copying_data=True):
        cdef RngState x = RngState()
        if is_copying_data:
            x.xx = self.xx
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def split(self, const cc.std_string& seed):
        cdef RngState x = RngState()
        x.xx = self.xx.split(seed)
        return x

    def rand_gen(self):
        """Generate a uniformly distributed random integer ranges from 0 up to 2**64 - 1"""
        return cc.rand_gen(self.xx)

    def u_rand_gen(self, double upper=1.0, double lower=0.0):
        return cc.u_rand_gen(self.xx, upper, lower)

    def g_rand_gen(self, double center=0.0, double sigma=1.0):
        return cc.g_rand_gen(self.xx, center, sigma)

    def c_rand_gen(self, Coordinate size):
        """``size`` can be ``total_site`` of the lattice"""
        cdef Coordinate x = Coordinate()
        x.xx = cc.c_rand_gen(self.xx, size.xx)
        return x

    def select(self, list l):
        ri = self.rand_gen() % len(l)
        return l[ri]

    def rand_arr(self, shape):
        """
        return numpy array with `shape` and dtype=np.uint64
        Uniformly distributed random integer ranges from 0 up to 2**64 - 1
        """
        arr = np.empty(shape, dtype=np.uint64)
        arr_ravel = arr.ravel()
        self.rand_fill_u_long(arr.ravel())
        return arr

    def u_rand_arr(self, shape):
        """
        return numpy array with `shape` and dtype=np.float64
        Uniform distribution from 0.0 to 1.0
        """
        arr = np.empty(shape, dtype=np.float64)
        self.u_rand_fill_real_d(arr.ravel(), 1.0, 0.0)
        return arr

    def g_rand_arr(self, shape):
        """
        return numpy array with `shape` and dtype=np.float64
        Gaussian distribution center=0.0 with sigma=1.0
        """
        arr = np.empty(shape, dtype=np.float64)
        self.g_rand_fill_real_d(arr.ravel(), 0.0, 1.0)
        return arr

    def u_rand_fill(self, arr, double upper=1.0, double lower=0.0):
        """
        obsolete
        """
        arr = arr.ravel()
        assert arr.base is not None
        arr = arr.view(np.float64)
        assert arr.base is not None
        return self.u_rand_fill_real_d(arr, upper, lower)

    def g_rand_fill(self, arr, double center=0.0, double sigma=1.0):
        """
        obsolete
        """
        arr = arr.ravel()
        assert arr.base is not None
        arr = arr.view(np.float64)
        assert arr.base is not None
        return self.g_rand_fill_real_d(arr, center, sigma)

    def rand_fill_u_long(self, cc.ULong[:] arr):
        cdef cc.Long size = arr.size
        cdef cc.Long i
        for i in range(size):
            arr[i] = cc.rand_gen(self.xx)

    def u_rand_fill_real_d(self, double[:] arr, double upper=1.0, double lower=0.0):
        cdef cc.Long size = arr.size
        cdef cc.Long i
        for i in range(size):
            arr[i] = cc.u_rand_gen(self.xx, upper, lower)

    def g_rand_fill_real_d(self, double[:] arr, double center=0.0, double sigma=1.0):
        cdef cc.Long size = arr.size
        cdef cc.Long i
        for i in range(size):
            arr[i] = cc.g_rand_gen(self.xx, center, sigma)

### -------------------------------------------------------------------

@timer
def get_data_sig(x, RngState rs):
    """
    Return a signature (a floating point number, real or complex) of data viewed as a 1-D array of numbers.\n
    Result only depends on the value of the data, not the structure.
    ``x`` can be an instance of ``LatData``, ``np.ndarray``, etc.
    """
    if isinstance(x, np.ndarray):
        rs = rs.copy()
        arr = x.ravel()
        arr_rand = rs.u_rand_arr(arr.shape) * 2.0 - 1.0
        return np.sum(arr * arr_rand)
    elif isinstance(x, (LatData, LatDataInt, LatDataRealF, LatDataLong)):
        return get_data_sig(np.asarray(x), rs)
    elif isinstance(x, (SpinMatrix, ColorMatrix, WilsonMatrix,)):
        return get_data_sig(np.asarray(x), rs)
    elif isinstance(x, (int, float,)):
        return float(x)
    elif isinstance(x, (complex,)):
        return complex(x)
    elif hasattr(x, "get_data_sig"):
        return x.get_data_sig(rs)
    else:
        return None

@timer
def random_permute(list l, RngState rs):
    """
    Do not change `l`.
    Return a new permuted list.
    """
    cdef cc.Long size = len(l)
    cdef cc.std_vector[cc.PyObject*] vec
    vec.resize(size)
    cdef cc.Long i = 0
    for i in range(size):
        vec[i] = <cc.PyObject*>l[i]
    cc.random_permute(vec, rs.xx)
    cdef list l_new = []
    for i in range(size):
        l_new.append(<object>vec[i])
    return l_new
