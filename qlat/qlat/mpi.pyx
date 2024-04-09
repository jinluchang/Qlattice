# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

cimport qlat_utils.everything
from qlat_utils.all cimport *
from . cimport everything as cc
from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

import qlat_utils as q
import numpy as np

### -------------------------------------------------------------------

def mpi_level_count():
    return cc.mpi_level_count()

def begin(int id_node, Coordinate size_node, int color=0):
    cc.begin(id_node, size_node.xx, color)

def end(cc.Bool is_preserving_cache=False):
    if not is_preserving_cache:
        q.clean_cache()
    cc.end(is_preserving_cache)
    if mpi_level_count() == 0:
        q.set_display_method()

### -------------------------------------------------------------------

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

def glb_sum_double(cc.RealD x):
    cc.glb_sum(x)
    return x

def glb_sum_complex(cc.PyComplexD x):
    cdef cc.ComplexD xx = cc.ccpy_d(x)
    cc.glb_sum(xx)
    x = cc.pycc_d(xx)
    return x

def glb_sum_lat_data_in_place(LatData ld):
    cc.glb_sum(ld.xx)
    return ld

def glb_sum_lat_data(LatData ld):
    cdef LatData ld1 = ld.copy()
    return glb_sum_lat_data_in_place(ld1)

@q.timer
def glb_sum_np(x):
    """
    x does NOT change
    """
    shape = x.shape
    dtype = x.dtype
    l = list(x.flatten())
    ld = LatData()
    if dtype == np.dtype('float64'):
        ld.from_list(l, is_complex=False)
    elif dtype == np.dtype('int64'):
        ld.from_list(list(map(float, l)), is_complex=False)
    elif dtype == np.dtype('complex128'):
        ld.from_list(l, is_complex=True)
    else:
        q.displayln(dtype)
        assert False
    ld.glb_sum_in_place()
    return np.array(ld.to_list(), dtype=dtype).reshape(shape)

@q.timer
def glb_sum(x):
    """
    x does NOT change
    """
    if isinstance(x, float):
        return glb_sum_double(x)
    elif isinstance(x, complex):
        return glb_sum_complex(x)
    elif isinstance(x, (int, np.int64)):
        return glb_sum_long(x)
    elif isinstance(x, np.ndarray):
        return glb_sum_np(x)
    elif isinstance(x, list):
        return [ glb_sum(x_i) for x_i in x ]
    elif isinstance(x, tuple):
        return tuple(glb_sum(x_i) for x_i in x)
    else:
        # possible types: q.LatData
        return x.glb_sum()

### -------------------------------------------------------------------
