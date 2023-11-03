# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc
from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

from .geometry cimport *
from .field_types cimport *

import cqlat as c
import qlat_utils as q
import numpy as np

### -------------------------------------------------------------------

cdef class PointsSelection:

    def __cinit__(self):
        self.cdata = <cc.Long>&(self.xx)
        self.geo = None
        self.view_count = 0

    def __init__(self, xg_arr=None, Geometry geo=None):
        """
        PointsSelection()
        PointsSelection(n_points, geo)
        PointsSelection(xg, geo)
        PointsSelection(xg_arr, geo)
        PointsSelection(xg_list, geo)
        """
        self.set_xg_arr(xg_arr, geo)

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef int ndim = 2
        cdef Buffer buf = Buffer(self, ndim)
        buf.format = 'i'
        buf.itemsize = sizeof(cc.Int)
        buf.buf = <char*>(self.xx.data())
        buf.set_dim_size(0, self.xx.size())
        buf.set_dim_size(1, 4)
        buf.update_strides_from_shape()
        buf.set_buffer(buffer, flags)
        self.view_count += 1

    def release_buffer(self, Buffer buf):
        assert buf.obj is self
        self.view_count -= 1

    def __imatmul__(self, PointsSelection v1 not None):
        if self.view_count > 0:
            raise Exception("PointsSelection.__imatmul__: self.view_count>0")
        self.geo = v1.geo
        cc.assign_direct(self.xx, v1.xx)
        return self

    def copy(self):
        x = type(self)()
        x @= self
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def set_rand(self, RngState rs not None, Coordinate total_site not None, cc.Long n_points):
        if self.view_count > 0:
            raise ValueError("can't re-init while being viewed")
        self.geo = Geometry(total_site)
        cc.assign_direct(self.xx, cc.mk_random_point_selection(total_site.xx, n_points, rs.xx))

    def save(self, const cc.std_string& path):
        cc.save_point_selection_info(self.xx, path)

    def load(self, const cc.std_string& path, Geometry geo):
        if self.view_count > 0:
            raise ValueError("can't re-init while being viewed")
        self.geo = geo
        cc.assign_direct(self.xx, cc.load_point_selection_info(path))

    def n_points(self):
        return self.xx.size()

    def set_n_points(self, const cc.Long n_points):
        self.xx.resize(n_points)

    def xg_arr(self):
        """
        return xg for all selected points
        shape = (psel.n_points(), 4,)
        """
        return np.asarray(self, dtype=np.int32)

    @q.timer
    def set_xg_arr(self, xg_arr=None, Geometry geo=None):
        """
        psel.set_xg_arr()
        psel.set_xg_arr(n_points, geo)
        psel.set_xg_arr(xg, geo)
        psel.set_xg_arr(xg_arr, geo)
        psel.set_xg_arr(xg_list, geo)
        """
        if self.view_count > 0:
            raise ValueError("can't re-init while being viewed")
        self.geo = geo
        cdef cc.Long n_points
        cdef cc.Long i
        cdef Coordinate xg
        cdef cc.Coordinate *p_xg
        if xg_arr is None:
            self.xx = cc.PointsSelection()
        elif isinstance(xg_arr, int):
            n_points = xg_arr
            self.set_n_points(n_points)
        elif isinstance(xg_arr, Coordinate):
            xg = xg_arr
            n_points = 1
            self.xx = cc.PointsSelection(n_points)
            cc.assign_direct(self.xx[0], xg.xx)
        elif isinstance(xg_arr, np.ndarray):
            n_points = len(xg_arr)
            self.xx = cc.PointsSelection(n_points)
            assert xg_arr.shape == (n_points, 4,)
            np.asarray(self, dtype=np.int32)[:] = xg_arr
        elif isinstance(xg_arr, list):
            n_points = len(xg_arr)
            self.xx = cc.PointsSelection(n_points)
            for i in range(n_points):
                self.xx[i] = Coordinate(xg_arr[i]).xx
        else:
            raise Exception(f"PointsSelection.set_xg_arr({xg_arr},{geo})")

    def set_geo(self, Geometry geo):
        self.geo = geo

    def __getitem__(self, cc.Long idx):
        cdef Coordinate xg = Coordinate()
        cc.assign_direct(xg.xx, self.xx[idx])
        return xg

    def __setitem__(self, cc.Long idx, Coordinate xg not None):
        cc.assign_direct(self.xx[idx], xg.xx)

    def __iter__(self):
        cdef cc.Long idx
        cdef cc.Long n_points = self.n_points()
        for idx in range(n_points):
            yield self[idx]

    def __len__(self):
        return self.n_points()

### -------------------------------------------------------------------

cdef class FieldSelection:

    def __cinit__(self):
        self.cdata = <cc.Long>&(self.xx)

    def __init__(self, Coordinate total_site=None, cc.Long n_per_tslice=-1, RngState rs=None, PointsSelection psel=None):
        if total_site is not None:
            assert rs is not None
            self.set_rand(rs, total_site, n_per_tslice)
            if psel is not None:
                self.add_psel(psel)
            self.update()

    def __imatmul__(self, FieldSelection v1):
        cc.assign_direct(self.xx, v1.xx)
        return self

    def copy(self):
        x = type(self)()
        x @= self
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def set_uniform(self, Coordinate total_site, val=0):
        """
        default (val = 0) select every sites
        val = -1 deselection everything
        """
        cc.mk_field_selection(self.xx.f_rank, total_site.xx, val)
        self.update()

    def set_rand(self, RngState rs, Coordinate total_site, cc.Long n_per_tslice):
        cc.mk_field_selection(self.xx.f_rank, total_site.xx, n_per_tslice, rs.xx)
        self.update()

    def add_psel(self, PointsSelection psel, cc.Long rank_psel=1024 * 1024 * 1024 * 1024 * 1024):
        """
        Add psel points to the selection, with the rank specified as rank_psel.
        If the point is already selected with lower rank, the rank is unchanged.
        """
        cc.add_field_selection(self.xx.f_rank, psel.xx, rank_psel)
        self.update()

    def update(self):
        """
        update various indices based on f_rank
        """
        cc.update_field_selection(self.xx)

    def to_psel(self):
        cdef PointsSelection psel = PointsSelection(None, self.geo())
        cc.assign_direct(psel.xx, cc.psel_from_fsel(self.xx))
        return psel

    def to_psel_local(self):
        cdef PointsSelection psel = PointsSelection(None, self.geo())
        cc.assign_direct(psel.xx, cc.psel_from_fsel_local(self.xx))
        return psel

    def save(self, const cc.std_string& path):
        return cc.write_field_selection(self.xx, path)

    def load(self, const cc.std_string& path):
        return cc.read_field_selection(self.xx, path)

    def geo(self):
        cdef Geometry geo = Geometry()
        geo.xx = self.xx.f_rank.get_geo()
        return geo

    def total_site(self):
        cdef Coordinate total_site = Coordinate()
        cc.assign_direct(total_site.xx, self.xx.f_rank.get_geo().total_site())
        return total_site

    def n_elems(self):
        return self.xx.n_elems

    def __getitem__(self, cc.Long idx):
        cdef Coordinate xg = Coordinate()
        cdef cc.Long index = self.xx.indices[idx]
        cc.assign_direct(xg.xx, self.xx.f_local_idx.get_geo().coordinate_from_index(index))
        return xg

    def idx_from_coordinate(self, Coordinate xg not None):
        cdef cc.Coordinate xl_xx = self.xx.f_local_idx.get_geo().coordinate_l_from_g(xg.xx)
        cdef cc.Long idx = self.xx.f_local_idx.get_elem(xl_xx)
        return idx

    def coordinate_from_idx(self, cc.Long idx):
        cdef Coordinate xg = Coordinate()
        cdef cc.Long index = self.xx.indices[idx]
        cc.assign_direct(xg.xx, self.xx.f_local_idx.get_geo().coordinate_from_index(index))
        return xg

### -------------------------------------------------------------------

cache_point_selection = q.mk_cache("point_selection")

@q.timer
def mk_xg_field(Geometry geo):
    cdef FieldInt f = FieldInt()
    cc.set_xg_field(f.xx, geo.xx)
    return f

def get_psel_single(Coordinate total_site, Coordinate xg=None):
    """
    [ xg, ]
    need total_site to set the psel.geo property
    """
    if xg is None:
        xg = Coordinate([ -1, -1, -1, -1, ])
    param_tuple = (total_site[0], total_site[1], total_site[2], total_site[3], xg[0], xg[1], xg[2], xg[3],)
    if param_tuple not in cache_point_selection:
        psel = PointsSelection([ xg, ], Geometry(total_site))
        cache_point_selection[param_tuple] = psel
    return cache_point_selection[param_tuple]

def get_psel_tslice(Coordinate total_site, *, int t_dir=3):
    """
    if t_dir = 3, then [ [-1,-1,-1,-1,], [-1,-1,-1,1,], ..., [-1,-1,-1,total_site[3]-1], ]
    if t_dir = 2, then [ [-1,-1,-1,-1,], [-1,-1,1,-1,], ..., [-1,-1,total_site[2]-1,-1], ]
    need total_site to set the psel.geo property
    """
    assert 0 <= t_dir and t_dir < 4
    param_tuple = (total_site[0], total_site[1], total_site[2], total_site[3], t_dir,)
    if param_tuple not in cache_point_selection:
        psel = PointsSelection(None, Geometry(total_site))
        cc.assign_direct(psel.xx, cc.mk_tslice_point_selection(total_site[t_dir], t_dir))
        cache_point_selection[param_tuple] = psel
    return cache_point_selection[param_tuple]

def is_matching_fsel(FieldSelection fsel1, FieldSelection fsel2):
    return cc.is_matching_fsel(fsel1.xx, fsel2.xx)

### -------------------------------------------------------------------
