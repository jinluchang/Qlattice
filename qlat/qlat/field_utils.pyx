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
        )

from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

import cqlat as c
import qlat_utils as q
import numpy as np
import math

from .geometry import geo_reform
from .field_base import (
        Field,
        SelectedField,
        SelectedPoints,
        )

def field_expanded(f, expansion_left, expansion_right):
    geo = f.geo()
    multiplicity = geo.multiplicity()
    geo_e = geo_reform(geo, multiplicity, expansion_left, expansion_right)
    f_e = type(f)(geo = geo_e)
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
        self.cdata = <long>&(self.xx)

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
    f = Field(ElemTypeComplex, geo, 1)
    c.set_phase_field(f, lmom)
    return f

class FastFourierTransform:

    def __init__(self, fft_infos, *, is_normalizing=False, mode_fft=1):
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
                total_site = field.total_site()
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
    if isinstance(f, FieldBase):
        f_n = Field(ElemTypeDouble)
        c.qnorm_field_field(f_n, f)
    elif isinstance(f, SelectedFieldBase):
        fsel = f.fsel
        f_n = SelectedField(ElemTypeDouble, fsel)
        c.qnorm_field_sfield(f_n, f)
    elif isinstance(f, SelectedPointsBase):
        psel = f.psel
        f_n = SelectedPoints(ElemTypeDouble, psel)
        c.qnorm_field_spfield(f_n, f)
    else:
        q.displayln_info("qnorm_field:", type(f))
        assert False
    return f_n

@q.timer
def sqrt_double_field(f):
    if isinstance(f, FieldBase):
        assert f.ctype is ElemTypeDouble
        f_ret = Field(ElemTypeDouble)
        c.set_sqrt_double_field(f_ret, f)
    elif isinstance(f, SelectedFieldBase):
        assert f.ctype == ElemTypeDouble
        fsel = f.fsel
        f_ret = SelectedField(ElemTypeDouble, fsel)
        c.set_sqrt_double_sfield(f_ret, f)
    elif isinstance(f, SelectedPointsBase):
        assert f.ctype == ElemTypeDouble
        psel = f.psel
        f_ret = SelectedPoints(ElemTypeDouble, psel)
        c.set_sqrt_double_spfield(f_ret, f)
    else:
        q.displayln_info("sqrt_double_field:", type(f))
        assert False
    return f_ret

###
