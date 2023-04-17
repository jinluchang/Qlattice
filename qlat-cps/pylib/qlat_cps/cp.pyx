# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat.all cimport *
from . cimport everything as cc

import sys
import qlat as q

def begin_with_cps(total_site):
    cdef Coordinate total_site_ = q.Coordinate(total_site)
    cc.begin_with_cps(sys.argv, total_site_.xx)

def end_with_cps(is_preserving_cache = False):
    if not is_preserving_cache:
        q.clean_cache()
    cc.end_with_cps(is_preserving_cache)

def save_cps_prop_double(prop, path):
    assert isinstance(prop, q.Prop)
    assert isinstance(path, str)
    cdef FieldWilsonMatrix field_wm = prop
    cc.save_cps_prop_double(field_wm.xx, path)

def load_cps_prop_double(prop, path):
    assert isinstance(prop, q.Prop)
    assert isinstance(path, str)
    cdef FieldWilsonMatrix field_wm = prop
    cc.load_cps_prop_double(field_wm.xx, path)
