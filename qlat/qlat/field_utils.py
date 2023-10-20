from qlat_utils import *
from .c import *

import numpy as np
import math

from . import c

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

class FieldExpandCommPlan:

    # self.cdata

    def __init__(self):
        self.cdata = c.mk_field_expand_comm_plan()

    def __del__(self):
        assert isinstance(self.cdata, int)
        c.free_field_expand_comm_plan(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, FieldExpandCommPlan)
        c.set_field_expand_comm_plan(self, v1)
        return self

    def copy(self, is_copying_data=True):
        x = FieldExpandCommPlan()
        if is_copying_data:
            x @= self
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

###

def make_field_expand_comm_plan(comm_marks):
    """
    comm_marks is of type Field(ElemTypeInt8t)
    """
    cp = FieldExpandCommPlan()
    c.make_field_expand_comm_plan(cp, comm_marks)
    return cp

def mk_phase_field(geo: Geometry, lmom):
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

@timer
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

@timer
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
        displayln_info("qnorm_field:", type(f))
        assert False
    return f_n

@timer
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
        displayln_info("sqrt_double_field:", type(f))
        assert False
    return f_ret

###

@timer
def set_selected_points(sp, f):
    # deprecated use @=
    from qlat.selected_field import SelectedField
    assert isinstance(sp, SelectedPointsBase)
    if isinstance(f, FieldBase):
        c.set_spfield_field(sp, f)
    elif isinstance(f, SelectedFieldBase):
        c.set_spfield_sfield(sp, f)
    else:
        raise Exception("set_selected_points")

@timer
def set_selected_field(sf, f):
    # deprecated use @=
    displayln_info("set_selected_field: deprecated")
    assert isinstance(sf, SelectedFieldBase)
    if isinstance(f, FieldBase):
        c.set_sfield_field(sf, f)
    elif isinstance(f, SelectedFieldBase):
        c.set_sfield_sfield(sf, f)
    else:
        raise Exception("set_selected_field")
