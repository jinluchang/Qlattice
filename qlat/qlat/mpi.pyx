# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

cimport qlat_utils.everything
from qlat_utils.all cimport *
from . cimport everything as cc
from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

import qlat_utils as q
import numpy as np

### -------------------------------------------------------------------

def begin(int id_node, Coordinate size_node, int color=0):
    cc.begin(id_node, size_node.xx, color)

def end(cc.bool is_preserving_cache=False):
    if not is_preserving_cache:
        q.clean_cache()
    cc.end(is_preserving_cache)

### -------------------------------------------------------------------

def sync_node():
    cc.sync_node()

def get_size_node():
    cdef Coordinate x = Coordinate()
    # cc.assign_direct(x.xx, cc.get_size_node())
    x.xx = cc.get_size_node()
    return x

def get_coor_node():
    cdef Coordinate x = Coordinate()
    # cc.assign_direct(x.xx, cc.get_coor_node())
    x.xx = cc.get_coor_node()
    return x

### -------------------------------------------------------------------

def bcast_long(cc.Long x, int root=0):
    cc.bcast(x, root)
    return x

def bcast_double(double x, int root=0):
    cc.bcast(x, root)
    return x

def bcast_complex(cc.PyComplexD x, int root=0):
    cc.bcast(cc.ccpy_d(x), root)
    return x

def bcast_lat_data_in_place(LatData ld, int root=0):
    cc.bcast(ld.xx, root)
    return ld

def bcast_lat_data(LatData ld, int root=0):
    cdef LatData ld1 = ld.copy()
    return bcast_lat_data_in_place(ld1, root)

### -------------------------------------------------------------------

def glb_sum_long(cc.Long x):
    cc.glb_sum(x)
    return x

def glb_sum_double(double x):
    cc.glb_sum(x)
    return x

def glb_sum_complex(cc.PyComplexD x):
    cc.glb_sum(cc.ccpy_d(x))
    return x

def glb_sum_lat_data_in_place(LatData ld):
    cc.glb_sum(ld.xx)
    return ld

def glb_sum_lat_data(LatData ld):
    cdef LatData ld1 = ld.copy()
    return glb_sum_lat_data_in_place(ld1)

### -------------------------------------------------------------------
