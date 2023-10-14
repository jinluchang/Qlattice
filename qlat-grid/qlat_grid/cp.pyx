# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat.all cimport *
from . cimport everything as cc

import sys
import qlat as q

def begin_with_grid(size_node_list = None):
    if size_node_list is None:
        size_node_list = []
    else:
        assert isinstance(size_node_list, list)
        size_node_list = list(map(q.Coordinate, size_node_list))
    size_node_list = size_node_list + q.default_size_node_list
    cdef long vec_size = len(size_node_list)
    cdef cc.std_vector[cc.Coordinate] size_node_vec = cc.std_vector[cc.Coordinate](vec_size)
    cdef long i
    cdef Coordinate size_node
    for i in range(vec_size):
        size_node = Coordinate(size_node_list[i])
        size_node_vec[i] = size_node.xx
    cc.begin_with_grid(sys.argv, size_node_vec)

def end_with_grid(is_preserving_cache = False):
    if not is_preserving_cache:
        q.clean_cache()
    cc.end_with_grid(is_preserving_cache)

def save_grid_prop_float(prop, path):
    assert isinstance(prop, q.Prop)
    assert isinstance(path, str)
    cdef FieldWilsonMatrix field_wm = prop
    cc.save_grid_prop_float(field_wm.xx, path)

def save_grid_prop_double(prop, path):
    assert isinstance(prop, q.Prop)
    assert isinstance(path, str)
    cdef FieldWilsonMatrix field_wm = prop
    cc.save_grid_prop_double(field_wm.xx, path)

def load_grid_prop_float(prop, path):
    assert isinstance(prop, q.Prop)
    assert isinstance(path, str)
    cdef FieldWilsonMatrix field_wm = prop
    cc.load_grid_prop_float(field_wm.xx, path)

def load_grid_prop_double(prop, path):
    assert isinstance(prop, q.Prop)
    assert isinstance(path, str)
    cdef FieldWilsonMatrix field_wm = prop
    cc.load_grid_prop_double(field_wm.xx, path)
