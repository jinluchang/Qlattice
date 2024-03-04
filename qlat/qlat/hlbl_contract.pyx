# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from . cimport everything as cc
from qlat_utils.all cimport *
cimport numpy
import numpy as np
import qlat_utils as q
from .field_selection cimport FieldSelection
from .selected_field_types cimport SelectedFieldRealD

cdef object convert_from_sl_table(const cc.SlTable& x):
    pass

def set_m_z_field_tag(
        FieldSelection fsel,
        Coordinate xg_x,
        Coordinate xg_y,
        const cc.RealD a,
        const cc.Int tag
        ):
    cdef SelectedFieldRealD smf_d = SelectedFieldRealD(fsel)
    cc.set_m_z_field_tag(smf_d.xx, fsel.xx, xg_x.xx, xg_y.xx, a, tag)
    return smf_d
