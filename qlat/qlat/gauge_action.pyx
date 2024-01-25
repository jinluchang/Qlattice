# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc

import qlat_utils as q

cdef class GaugeAction:

    def __cinit__(self):
        self.xx = cc.GaugeAction()
        self.cdata = <cc.Long>&(self.xx)

    def __init__(self, cc.RealD beta, cc.RealD c1=0.0):
        self.xx.beta = beta
        self.xx.c1 = c1

    def __imatmul__(self, GaugeAction v1):
        cc.assign_direct(self.xx, v1.xx)
        return self

    def beta(self):
        return self.xx.beta

    def c1(self):
        return self.xx.c1
