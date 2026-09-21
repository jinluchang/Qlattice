# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.mpi``
===================\n
Low-level MPI utilities for lattice QCD simulations: MPI communicator begin/end,
node layout queries, broadcasting, and global sum reductions.\n
Documentation: ``docs/qlat/qlat_mpi.md``\n
.. note:: Update the documentation when updating this source file.
"""

cimport qlat_utils.everything
from qlat_utils.all cimport *
from . cimport everything as cc
from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

class q:
    from qlat_utils import (
        clean_cache,
        set_display_method,
    )

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

def is_initialized():
    """
    Return whether the global geometry node (`geon`) is initialized, i.e.
    whether `q.begin_with_mpi()` or `q.begin(id_node, size_node)` has been
    called.
    """
    return cc.get_geometry_node().initialized

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
    cdef cc.ComplexD xx = cc.ccpy_d(x)
    cc.bcast(xx, root)
    x = cc.pycc_d(xx)
    return x

def bcast_lat_data_in_place(LatData ld, int root=0):
    cc.bcast(ld.xx, root)
    return ld

def bcast_lat_data(LatData ld, int root=0):
    cdef LatData ld1 = ld.copy()
    return bcast_lat_data_in_place(ld1, root)

### -------------------------------------------------------------------

def get_id_node_list_for_shuffle():
    """
    return ``list``
    ``list`` is a list of ``id_node``s, which are intended to perform IO.
    ``list[id_node_in_shuffle] = id_node``
    """
    cdef list id_node_list = []
    cdef cc.std_vector[cc.Int] id_node_vec = cc.get_id_node_list_for_shuffle()
    for id_node in id_node_vec:
        id_node_list.append(id_node)
    return id_node_list

def get_id_node_in_shuffle_list():
    """
    return ``list``
    ``list[id_node] = id_node_in_shuffle``
    Related to ``get_id_node_list_for_shuffle()``
    """
    cdef list id_node_in_shuffle_list = []
    cdef cc.std_vector[cc.Int] id_node_in_shuffle_vec = cc.get_id_node_in_shuffle_list()
    for id_node in id_node_in_shuffle_vec:
        id_node_in_shuffle_list.append(id_node)
    return id_node_in_shuffle_list

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

### -------------------------------------------------------------------

# Re-exported from ``qlat.mpi_utils`` (where the pure-Python ``glb_sum``
# dispatch now lives) to keep ``qlat.mpi.glb_sum`` working.
from .mpi_utils import (
    glb_sum,
    glb_sum_np,
    )

