# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc
from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

import cqlat as c
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

cdef class Geometry:

    def __cinit__(self):
        self.xx = cc.Geometry()
        self.cdata = <long>&(self.xx)

    def __init__(self, total_site=None, int multiplicity=1):
        """
        if total_site is None: create geo uninitialized
        elif multiplicity is None: create geo with multiplicity = 1
        """
        self.xx.init(Coordinate(total_site).xx, multiplicity)

    def __imatmul__(self, Geometry v1):
        self.xx = v1.xx
        return self

    def copy(self, is_copying_data=True):
        x = Geometry()
        if is_copying_data:
            x @= self
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def total_site(self):
        return c.get_total_site_geo(self)

    def total_volume(self):
        return c.get_total_volume_geo(self)

    def local_site(self):
        return c.get_node_site_geo(self)

    def local_volume(self):
        return c.get_local_volume_geo(self)

    def multiplicity(self):
        return c.get_multiplicity_geo(self)

    def eo(self):
        return c.get_eo_geo(self)

    def expansion_left(self):
        return c.get_expansion_left_geo(self)

    def expansion_right(self):
        return c.get_expansion_right_geo(self)

    def id_node(self):
        return c.get_id_node_geo(self)

    def num_node(self):
        return c.get_num_node_geo(self)

    def coor_node(self):
        return c.get_coor_node_geo(self)

    def size_node(self):
        return c.get_size_node_geo(self)

    def __str__(self):
        return self.show()

    def __repr__(self):
        return self.show()

    def show(self):
        total_site = self.total_site()
        multiplicity = self.multiplicity()
        expan_left = self.expansion_left()
        expan_right = self.expansion_right()
        eo = self.eo()
        zero = [ 0, 0, 0, 0, ]
        if expan_left == zero and expan_right == zero and eo == 0:
            return f"Geometry({total_site}, {multiplicity})"
        else:
            return f"Geometry({total_site}, {multiplicity}, expan_left={expan_left}, expan_right={expan_right}, eo={eo})"

    def coordinate_g_from_l(self, xl):
        return c.coordinate_g_from_l_geo(self, xl)

    def coordinate_l_from_g(self, xg):
        return c.coordinate_l_from_g_geo(self, xg)

    def is_local(self, xl):
        return c.is_local_geo(self, xl)

    def is_local_xg(self, xg):
        """
        return a global coordinate is inside the local volume or not
        """
        return c.is_local_xg_geo(self, xg)

    def xg_list(self):
        """
        return xg for all local sites
        """
        return c.get_xg_list(self)

### -------------------------------------------------------------------

def geo_reform(Geometry geo, int multiplicity=1, expansion_left=None, expansion_right=None):
    cdef Coordinate el, er
    if expansion_left is None:
        el = Coordinate()
    elif isinstance(expansion_left, int):
        e = expansion_left
        el = Coordinate([ e, e, e, e, ])
    elif isinstance(expansion_left, list):
        el = Coordinate(expansion_left)
    else:
        assert isinstance(expansion_left, Coordinate)
        el = <Coordinate>expansion_left
    if expansion_right is None:
        er = Coordinate()
    elif isinstance(expansion_right, int):
        e = expansion_right
        er = Coordinate([ e, e, e, e, ])
    elif isinstance(expansion_right, list):
        er = Coordinate(expansion_right)
    else:
        assert isinstance(expansion_left, Coordinate)
        er = <Coordinate>expansion_left
    cdef Geometry geo_new = Geometry()
    geo_new.xx = cc.geo_reform(geo.xx, multiplicity, el.xx, er.xx)
    return geo_new

def geo_eo(Geometry geo, int eo=0):
    cdef Geometry geo_new = Geometry()
    geo_new.xx = cc.geo_eo(geo.xx, eo)
    return geo_new

### -------------------------------------------------------------------


