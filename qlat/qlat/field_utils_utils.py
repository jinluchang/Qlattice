"""
Module ``qlat.field_utils_utils``
==================================\n
Pure-Python field utilities that do not call C++ functions directly: field
expansion/refresh, the Fourier transform factory (``FastFourierTransform``
and ``mk_fft``), norms, and the element-wise square root dispatch.\n
"""

import math

import qlat_utils as q

from .geometry import geo_resize
from .field_base import (
        FieldBase,
        SelectedFieldBase,
        SelectedPointsBase,
        )
from .field_types import FieldRealD
from .selected_field_types import SelectedFieldRealD
from .selected_points_types import SelectedPointsRealD
from .field_utils import (
        sqrt_field_real_d,
        sqrt_selected_field_real_d,
        sqrt_selected_points_real_d,
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
        field._cc_refresh_expanded()
    else:
        field._cc_refresh_expanded_plan(comm_plan)

def refresh_expanded_1(field):
    field._cc_refresh_expanded_1()

def refresh_expanded_field(field, comm_plan=None):
    """
    cqlat-compatible name for ``refresh_expanded``.
    """
    refresh_expanded(field, comm_plan)

def refresh_expanded_1_field(field):
    """
    cqlat-compatible name for ``refresh_expanded_1``.
    """
    refresh_expanded_1(field)

def merge_fields_ms_field(f, fs, ms):
    """
    cqlat-compatible name for ``Field._cc_merge_fields_ms``.
    """
    assert isinstance(f, FieldBase)
    f._cc_merge_fields_ms(fs, ms)

### -------------------------------------------------------------------

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
        fields[0]._cc_fft(fields, fft_dirs, fft_is_forwards, self.mode_fft)
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
