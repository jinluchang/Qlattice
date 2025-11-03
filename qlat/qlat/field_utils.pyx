# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc
from .geometry cimport Geometry
from .field_base cimport (
        FieldBase,
        SelectedFieldBase,
        SelectedPointsBase,
        )
from .field_types cimport (
        FieldInt8t,
        FieldRealD,
        )
from .selected_field_types cimport SelectedFieldRealD
from .selected_points_types cimport SelectedPointsRealD

from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

import cqlat as c
import qlat_utils as q
import numpy as np
import math

from .geometry import geo_resize
from .field_base import (
        Field,
        SelectedField,
        SelectedPoints,
        )

def field_expanded(f, expansion_left, expansion_right):
    geo = f.geo
    multiplicity = f.multiplicity
    geo_e = geo_resize(geo, expansion_left, expansion_right)
    f_e = type(f)(geo_e, multiplicity)
    f_e @= f
    return f_e

def refresh_expanded(field, comm_plan=None):
    if comm_plan is None:
        return c.refresh_expanded_field(field)
    else:
        return c.refresh_expanded_field(field, comm_plan)

def refresh_expanded_1(field):
    return c.refresh_expanded_1_field(field)


cdef class CommMarks(FieldInt8t):

    def __init__(self, Geometry geo=None, int multiplicity=0):
        super().__init__(geo, multiplicity)

###

cdef class CommPlan:

    def __cinit__(self):
        self.cdata = <cc.Long>&(self.xx)

    def __imatmul__(self, CommPlan v1):
        self.xx = v1.xx
        return self

    def copy(self, is_copying_data=True):
        x = type(self)()
        if is_copying_data:
            x @= self
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

###

def make_field_expand_comm_plan(CommMarks comm_marks):
    """
    comm_marks is of type Field(ElemTypeInt8t)
    """
    cp = CommPlan()
    c.make_field_expand_comm_plan(cp, comm_marks)
    return cp

def mk_phase_field(Geometry geo, lmom):
    """
    lmom is in lattice momentum unit
    exp(i * 2*pi/L * lmom \cdot xg )
    """
    f = Field(ElemTypeComplexD, geo, 1)
    c.set_phase_field(f, lmom)
    return f

class FastFourierTransform:

    def __init__(self, fft_infos, *, is_normalizing=False, mode_fft=1):
        # mode_fft in [ 0, 1, ]
        # fft_infos = [ ( fft_dir, is_forward, ), ... ]
        self.fft_infos = fft_infos
        self.is_normalizing = is_normalizing
        self.mode_fft = mode_fft

    def copy(self):
        return self.__copy__()

    def __mul__(self, fields):
        if isinstance(fields, FieldBase):
            return (self * [ fields, ])[0]
        assert isinstance(fields, list)
        for f in fields:
            assert isinstance(f, FieldBase)
        fields = [ f.copy() for f in fields ]
        fft_dirs, fft_is_forwards = zip(*self.fft_infos)
        c.fft_fields(fields, fft_dirs, fft_is_forwards, self.mode_fft)
        if self.is_normalizing and self.fft_infos:
            for field in fields:
                total_site = field.total_site
                scale_factor = 1
                for fft_dir, is_forward in self.fft_infos:
                    scale_factor *= total_site[fft_dir]
                scale_factor = 1.0 / math.sqrt(scale_factor)
                field *= scale_factor
        return fields

###

@q.timer
def mk_fft(is_forward, *, is_only_spatial=False, is_normalizing=False, mode_fft=1):
    if is_only_spatial:
        fft_infos = [
                (0, is_forward,),
                (1, is_forward,),
                (2, is_forward,),
                ]
        return FastFourierTransform(fft_infos, is_normalizing=is_normalizing, mode_fft=mode_fft)
    else:
        fft_infos = [
                (0, is_forward,),
                (1, is_forward,),
                (2, is_forward,),
                (3, is_forward,),
                ]
        return FastFourierTransform(fft_infos, is_normalizing=is_normalizing, mode_fft=mode_fft)

###

@q.timer
def qnorm_field(f):
    if isinstance(f, (FieldBase, SelectedFieldBase, SelectedPointsBase,)):
        f_n = f.qnorm_field()
    else:
        q.displayln_info("qnorm_field:", type(f))
        assert False
    return f_n

@q.timer
def sqrt_selected_points_real_d(SelectedPointsRealD f):
    cdef SelectedPointsRealD f_ret = f.copy(is_copying_data=False)
    cc.set_sqrt_field(f_ret.xx, f.xx)
    return f_ret

@q.timer
def sqrt_selected_field_real_d(SelectedFieldRealD f):
    cdef SelectedFieldRealD f_ret = f.copy(is_copying_data=False)
    cc.set_sqrt_field(f_ret.xx, f.xx)
    return f_ret

@q.timer
def sqrt_field_real_d(FieldRealD f):
    cdef FieldRealD f_ret = f.copy(is_copying_data=False)
    cc.set_sqrt_field(f_ret.xx, f.xx)
    return f_ret

@q.timer
def sqrt_field(f):
    if isinstance(f, FieldRealD):
        f_ret = sqrt_field_real_d(f)
    elif isinstance(f, SelectedFieldRealD):
        f_ret = sqrt_selected_field_real_d(f)
    elif isinstance(f, SelectedPointsRealD):
        f_ret = sqrt_selected_points_real_d(f)
    else:
        q.displayln_info("sqrt_field:", type(f))
        assert False
    return f_ret

###