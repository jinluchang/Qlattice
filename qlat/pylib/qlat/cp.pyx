# cython: c_string_type=unicode, c_string_encoding=utf8

from . cimport everything as cp

def qremove_info(path):
    return cp.qremove_info(path)

def qremove_all_info(path):
    return cp.qremove_all_info(path)

cdef class SpinMatrix:
    cdef cp.SpinMatrix sm

    def __cinit__(self):
        self.sm = cp.SpinMatrix()

    cpdef SpinMatrix mul_sm(self, SpinMatrix other):
        cdef SpinMatrix ret = SpinMatrix()
        ret.sm = self.sm * other.sm
        return ret
