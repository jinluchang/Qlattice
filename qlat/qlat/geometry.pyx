# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc
from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

import qlat_utils as q
import numpy as np

### -------------------------------------------------------------------

cdef class Geometry:

    def __cinit__(self):
        self.xx = cc.Geometry()
        self.cdata = <cc.Long>&(self.xx)

    def __init__(self, *args):
        """
        Geometry()
        Geometry(total_site)
        Geometry(id_node, size_node, node_site)
        Geometry(coor_node, size_node, node_site)
        """
        n_args = len(args)
        if n_args == 0:
            self.xx.init()
        elif n_args == 1:
            if isinstance(args[0], Coordinate):
                self.init_from_total_site(*args)
            else:
                fname = q.get_fname()
                raise Exception(f"{fname}: args={args}")
        elif n_args == 3:
            if isinstance(args[0], Coordinate):
                self.init_from_coor_node_site(*args)
            elif isinstance(args[0], int):
                self.init_from_id_node_site(*args)
            else:
                fname = q.get_fname()
                raise Exception(f"{fname}: args={args}")
        else:
            fname = q.get_fname()
            raise Exception(f"{fname}: args={args}")

    def init_from_total_site(self, Coordinate total_site):
        """
        Use `total_site` to initialize `Geometry`.
        Will use the global `geon` to determine `id_node` and `size_node`.
        """
        self.xx.init(total_site.xx)

    def init_from_id_node_site(self, int id_node, Coordinate size_node, Coordinate node_site):
        """
        initialize Geometry without call MPI functions.
        #
        corr_node = q.coordinate_from_index(id_node, size_node)
        id_node = q.index_from_coordinate(corr_node, size_node)
        #
        total_site = size_node * node_site
        """
        self.xx.init(id_node, size_node.xx, node_site.xx)

    def init_from_coor_node_site(self, Coordinate coor_node, Coordinate size_node, Coordinate node_site):
        """
        initialize Geometry without call MPI functions.
        #
        corr_node = q.coordinate_from_index(id_node, size_node)
        id_node = q.index_from_coordinate(corr_node, size_node)
        #
        total_site = size_node * node_site
        """
        self.xx.init(coor_node.xx, size_node.xx, node_site.xx)

    def __imatmul__(self, Geometry v1):
        cc.assign_direct(self.xx, v1.xx)
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

    def __eq__(self, Geometry other):
        return self.xx == other.xx

    @property
    def total_site(self):
        cdef Coordinate x = Coordinate()
        cc.assign_direct(x.xx, self.xx.total_site())
        return x

    @property
    def total_volume(self):
        return self.xx.total_volume()

    @property
    def local_site(self):
        cdef Coordinate x = Coordinate()
        cc.assign_direct(x.xx, self.xx.local_site())
        return x

    @property
    def local_volume(self):
        return self.xx.local_volume()

    @property
    def local_volume_expanded(self):
        return self.xx.local_volume_expanded()

    @property
    def eo(self):
        return self.xx.eo

    @property
    def expansion_left(self):
        cdef Coordinate x = Coordinate()
        cc.assign_direct(x.xx, self.xx.expansion_left)
        return x

    @property
    def expansion_right(self):
        cdef Coordinate x = Coordinate()
        cc.assign_direct(x.xx, self.xx.expansion_right)
        return x

    @property
    def id_node(self):
        return self.xx.geon.id_node

    @property
    def num_node(self):
        return self.xx.geon.num_node

    @property
    def coor_node(self):
        cdef Coordinate x = Coordinate()
        cc.assign_direct(x.xx, self.xx.geon.coor_node)
        return x

    @property
    def size_node(self):
        cdef Coordinate x = Coordinate()
        cc.assign_direct(x.xx, self.xx.geon.size_node)
        return x

    @property
    def node_site(self):
        cdef Coordinate x = Coordinate()
        cc.assign_direct(x.xx, self.xx.node_site)
        return x

    @property
    def is_only_local(self):
        return self.xx.is_only_local

    def __str__(self):
        return self.show()

    def __repr__(self):
        return self.show()

    def show(self):
        cdef Coordinate total_site = self.total_site
        cdef Coordinate expan_left = self.expansion_left
        cdef Coordinate expan_right = self.expansion_right
        cdef int eo = self.eo
        cdef Coordinate zero = Coordinate()
        if expan_left == zero and expan_right == zero and eo == 0:
            return f"Geometry({total_site.to_list()})"
        else:
            return f"Geometry({total_site.to_list()}, expan_left={expan_left.to_list()}, expan_right={expan_right.to_list()}, eo={eo})"

    def coordinate_g_from_l(self, Coordinate xl not None):
        cdef Coordinate xg = Coordinate()
        cc.assign_direct(xg.xx, self.xx.coordinate_g_from_l(xl.xx))
        return xg

    def coordinate_l_from_g(self, Coordinate xg not None):
        cdef Coordinate xl = Coordinate()
        cc.assign_direct(xl.xx, self.xx.coordinate_l_from_g(xg.xx))
        return xl

    def index_from_coordinate(self, Coordinate xl not None):
        return self.xx.index_from_coordinate(xl.xx)

    def coordinate_from_index(self, cc.Long index):
        cdef Coordinate xl = Coordinate()
        cc.assign_direct(xl.xx, self.xx.coordinate_from_index(index))
        return xl

    def index_from_g_coordinate(self, Coordinate xg not None):
        return self.xx.index_from_g_coordinate(xg.xx)

    def g_coordinate_from_index(self, cc.Long index):
        cdef Coordinate xg = Coordinate()
        cc.assign_direct(xg.xx, self.xx.g_coordinate_from_index(index))
        return xg

    def is_local(self, Coordinate xl not None):
        return self.xx.is_local(xl.xx)

    def is_local_xg(self, Coordinate xg not None):
        """
        return a global coordinate is inside the local volume or not
        """
        cdef cc.Coordinate xl_xx = self.xx.coordinate_l_from_g(xg.xx)
        return self.xx.is_local(xl_xx)

    def xg_arr(self):
        """
        return xg for all local sites
        shape = (geo.local_volume, 4,)
        """
        from .field_selection import mk_xg_field
        f_xg = mk_xg_field(self)
        xg_arr = np.asarray(f_xg, dtype=np.int32)
        cdef cc.Long local_volume = self.local_volume
        return xg_arr.reshape((local_volume, 4,))

    def __getstate__(self):
        """
        State of `Geometry` is different on different nodes.
        """
        size_node = self.size_node
        coor_node = self.coor_node
        node_site = self.node_site
        expan_left = self.expansion_left
        expan_right = self.expansion_right
        return [ size_node, coor_node, node_site, expan_left, expan_right, ]

    def __setstate__(self, state):
        """
        State of `Geometry` is different on different nodes.
        """
        self.__init__()
        cdef Coordinate size_node
        cdef Coordinate coor_node
        cdef Coordinate node_site
        cdef Coordinate expan_left
        cdef Coordinate expan_right
        [ size_node, coor_node, node_site, expan_left, expan_right, ] = state
        self.xx.init(coor_node.xx, size_node.xx, node_site.xx)
        self.xx.expansion_left = expan_left.xx
        self.xx.expansion_right = expan_right.xx
        self.xx.reset_node_site_expanded()

### -------------------------------------------------------------------

def geo_resize(Geometry geo, expansion_left=None, expansion_right=None):
    """
    expansion_left can be None, or int, or Coordinate
    expansion_right can be None, or int, or Coordinate
    Default is zero.
    """
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
        assert isinstance(expansion_right, Coordinate)
        er = <Coordinate>expansion_right
    cdef Geometry geo_new = Geometry()
    cc.assign_direct(geo_new.xx, cc.geo_resize(geo.xx, el.xx, er.xx))
    return geo_new

def geo_eo(Geometry geo, int eo=0):
    cdef Geometry geo_new = Geometry()
    cc.assign_direct(geo_new.xx, cc.geo_eo(geo.xx, eo))
    return geo_new

### -------------------------------------------------------------------
