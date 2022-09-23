import qlat.cqlat as c

import numpy as np

class WilsonMatrix:

    # self.cdata
    # self.base

    def __init__(self, *, cdata = None, base = None):
        # self.cdata will be newly allocated if not assigned
        if cdata is None:
            assert base is None
            self.cdata = c.mk_wilson_matrix()
            self.base = None
        else:
            assert base is not None
            self.cdata = cdata
            self.base = base

    def __del__(self):
        # only free is base is None
        assert isinstance(self.cdata, int)
        if self.base is None:
            c.free_wilson_matrix(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, WilsonMatrix)
        c.set_wilson_matrix(self, v1)
        return self

    def copy(self, is_copying_data = True):
        x = WilsonMatrix()
        if is_copying_data:
            x @= self
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def __getstate__(self):
        # will not pickle to file correct if base is not None:
        if self.base is not None:
            return self.cdata
        return self.get_value()

    def __setstate__(self, arr):
        # will not pickle to file correct if base is not None:
        if isinstance(arr, int):
            return self.__init__(cdata = arr, base = "foreign")
        self.__init__()
        return self.set_value(arr)

    def set_zero(self):
        return c.set_zero_wilson_matrix(self)

    def set_value(self, value_as_numpy_array):
        return c.set_state_wm(self, value_as_numpy_array.tobytes())

    def get_value(self):
        # return a 12x12 2-D numpy array contain the data
        arr = np.frombuffer(c.get_state_wm(self), dtype = complex).reshape(12, 12)
        return arr

    def g5_herm(self):
        wm = WilsonMatrix()
        c.set_g5_herm_wilson_matrix(wm, self)
        return wm

    def __mul__(self, other):
        if isinstance(other, (int, float, complex,)):
            return mat_mul_a_sc(other, self)
        elif isinstance(other, SpinMatrix):
            return mat_mul_sc_s(self, other)
        elif isinstance(other, WilsonMatrix):
            return mat_mul_sc_sc(self, other)
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, (int, float, complex,)):
            return mat_mul_a_sc(other, self)
        elif isinstance(other, SpinMatrix):
            return mat_mul_s_sc(other, self)
        elif isinstance(other, WilsonMatrix):
            return mat_mul_sc_sc(other, self)
        else:
            return NotImplemented

###

class SpinMatrix:

    # self.cdata

    def __init__(self):
        self.cdata = c.mk_spin_matrix()

    def __del__(self):
        assert isinstance(self.cdata, int)
        c.free_spin_matrix(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, SpinMatrix)
        c.set_spin_matrix(self, v1)
        return self

    def copy(self, is_copying_data = True):
        x = SpinMatrix()
        if is_copying_data:
            x @= self
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def set_zero(self):
        return c.set_zero_spin_matrix(self)

    def set_value(self, value_as_numpy_array):
        return c.set_state_sm(self, value_as_numpy_array.tobytes())

    def get_value(self):
        # return a 4x4 2-D numpy array contain the data
        arr = np.frombuffer(c.get_state_sm(self), dtype = complex).reshape(4, 4)
        return arr

    def __mul__(self, other):
        if isinstance(other, (int, float, complex,)):
            return mat_mul_a_s(other, self)
        elif isinstance(other, SpinMatrix):
            return mat_mul_s_s(self, other)
        else:
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, (int, float, complex,)):
            return mat_mul_a_s(other, self)
        elif isinstance(other, SpinMatrix):
            return mat_mul_s_s(other, self)
        else:
            return NotImplemented

###

def mat_mul_sc_sc(x, y):
    wm = WilsonMatrix()
    c.set_wm_mul_wm_wm(wm, x, y)
    return wm

def mat_mul_sc_s(x, y):
    wm = WilsonMatrix()
    c.set_wm_mul_wm_sm(wm, x, y)
    return wm

def mat_mul_s_sc(x, y):
    wm = WilsonMatrix()
    c.set_wm_mul_sm_wm(wm, x, y)
    return wm

def mat_mul_s_s(x, y):
    sm = SpinMatrix()
    c.set_sm_mul_sm_sm(sm, x, y)
    return sm

def mat_mul_a_s(coef, x):
    sm = SpinMatrix()
    c.set_sm_mul_a_sm(sm, coef, x)
    return sm

def mat_mul_a_sc(coef, x):
    wm = SpinMatrix()
    c.set_wm_mul_a_wm(wm, coef, x)
    return wm

def mat_sc_trace(x):
    return c.trace_wm(x)

def mat_sc_sc_trace(x, y):
    return c.trace_wm_wm(x, y)

def mat_sc_s_trace(x, y):
    return c.trace_sm_wm(y, x)

def mat_s_sc_trace(x, y):
    return c.trace_sm_wm(x, y)
