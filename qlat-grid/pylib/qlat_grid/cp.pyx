# cython: c_string_type=unicode, c_string_encoding=utf8

from . cimport everything as cc

def end_with_grid(is_preserving_cache = False):
    cc.grid_end(is_preserving_cache)
