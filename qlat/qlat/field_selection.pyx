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
        self.cdata = <long>&(self.xx)
        self.geo = None
        self.view_count = 0

    def __init__(self, xg_arr=None, Geometry geo=None):
        """
        PointsSelection()
        PointsSelection(n_points, geo)
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

    def set_rand(self, RngState rs not None, Coordinate total_site not None, long n_points):
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

    def set_n_points(self, const long n_points):
        self.xx.resize(n_points)

    def xg_arr(self):
        """
        return xg for all selected points
        shape = (psel.n_points(), 4,)
        """
        return np.asarray(self, dtype=np.int32)

    def set_xg_arr(self, xg_arr=None, Geometry geo=None):
        """
        psel.set_xg_arr()
        psel.set_xg_arr(n_points, geo)
        psel.set_xg_arr(xg_arr, geo)
        psel.set_xg_arr(xg_list, geo)
        """
        if self.view_count > 0:
            raise ValueError("can't re-init while being viewed")
        self.geo = geo
        cdef long n_points
        cdef long i
        if xg_arr is None:
            self.xx = cc.PointsSelection()
        elif isinstance(xg_arr, int):
            n_points = xg_arr
            self.set_n_points(n_points)
        elif isinstance(xg_arr, np.ndarray):
            n_points = len(xg_arr)
            self.xx = cc.PointsSelection(n_points)
            assert xg_arr.shape == (n_points, 4,)
            np.asarray(self, dtype=np.int32)[:] = xg_arr
        elif isinstance(xg_arr, list):
            n_points = len(xg_arr)
            self.xx = cc.PointsSelection(n_points)
            for i in range(n_points):
                cc.assign_direct(self.xx[i], Coordinate(xg_arr[i]).xx)
        else:
            raise Exception(f"PointsSelection.set_xg_arr({xg_arr},{geo})")

    def set_geo(self, Geometry geo):
        self.geo = geo

    def __getitem__(self, long idx):
        cdef Coordinate xg = Coordinate()
        cc.assign_direct(xg.xx, self.xx[idx])
        return xg

    def __setitem__(self, long idx, Coordinate xg not None):
        cc.assign_direct(self.xx[idx], xg.xx)

    def __iter__(self):
        cdef long idx
        cdef long n_points = self.n_points()
        for idx in range(n_points):
            yield self[idx]

    def __len__(self):
        return self.n_points()

    # def coordinate_from_idx(self, long idx):
    #     print("CHECK: WARNING: PointsSelection.coordinate_from_idx")
    #     return self[idx]

    # def to_list(self):
    #     print("CHECK: WARNING: PointsSelection.to_list")
    #     return [ self.coordinate_from_idx(idx).to_list() for idx in range(self.n_points()) ]

    # def from_list(self, xg_arr=None, Geometry geo=None):
    #     print("CHECK: WARNING: PointsSelection.from_list")
    #     self.set_xg_arr(xg_arr, geo)

### -------------------------------------------------------------------

cdef class FieldSelection:

    def __cinit__(self):
        self.cdata = <long>&(self.xx)

    def __init__(self, Coordinate total_site=None, long n_per_tslice=-1, RngState rs=None, PointsSelection psel=None):
        if total_site is not None:
            assert rs is not None
            c.set_rand_fsel(self, rs, total_site.to_list(), n_per_tslice)
            if psel is not None:
                c.add_psel_fsel(self, psel)
            self.update()
            self.update(n_per_tslice)

    def __imatmul__(self, FieldSelection v1):
        self.xx = v1.xx
        return self

    def copy(self):
        x = FieldSelection()
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
        c.set_uniform_fsel(self, total_site.to_list(), val)

    def set_rand(self, RngState rs, Coordinate total_site, long n_per_tslice):
        assert isinstance(rs, RngState)
        c.set_rand_fsel(self, rs, total_site.to_list(), n_per_tslice)
        self.update()
        self.update(n_per_tslice)

    def add_psel(self, PointsSelection psel, long rank_psel=1024 * 1024 * 1024 * 1024 * 1024):
        """
        Add psel points to the selection, with the rank specified as rank_psel.
        If the point is already selected with lower rank, the rank is unchanged.
        """
        c.add_psel_fsel(self, psel, rank_psel)
        self.update()

    def update(self, long n_per_tslice=-1):
        """
        if n_per_tslice < 0: only update various indices
        if n_per_tslice >= 0: only update parameters (n_per_tslice and prob)
        """
        c.update_fsel(self, n_per_tslice)

    def select_rank_range(self, long rank_start=0, long rank_stop=-1):
        """
        return new fsel with selected points that
        rank_start <= rank and (rank < rank_stop or rank_stop == -1)
        Does NOT change the n_per_tslice parameter for the new fsel
        """
        fsel = FieldSelection()
        c.select_rank_range_fsel(fsel, self, rank_start, rank_stop)
        fsel.update()
        fsel.update(self.n_per_tslice())
        return fsel

    def select_t_range(self, long rank_start=0, long rank_stop=-1):
        """
        return new fsel with selected points that
        t_start <= t and (t < t_stop or t_stop == -1)
        rank_start <= rank < rank_stop (rank_stop = -1 implies unlimited)
        Does NOT change the n_per_tslice parameter for the new fsel
        """
        fsel = FieldSelection()
        c.select_rank_range_fsel(fsel, self, rank_start, rank_stop)
        fsel.update()
        fsel.update(self.n_per_tslice())
        return fsel

    def to_psel(self):
        psel = PointsSelection(None, self.geo())
        c.set_psel_fsel(psel, self)
        return psel

    def to_psel_local(self):
        psel = PointsSelection(None, self.geo())
        c.set_psel_fsel_local(psel, self)
        return psel

    def save(self, str path):
        return c.save_fsel(self, path)

    def load(self, str path, long n_per_tslice):
        return c.load_fsel(self, path, n_per_tslice)

    def geo(self):
        geo = Geometry()
        c.set_geo_fsel(geo, self)
        return geo

    def total_site(self):
        return c.get_total_site_fsel(self)

    def n_elems(self):
        return c.get_n_elems_fsel(self)

    def n_per_tslice(self):
        return c.get_n_per_tslice_fsel(self)

    def prob(self):
        """
        return fsel.prob
        n_per_tslice / spatial_volume
        """
        return c.get_prob_fsel(self)

    def idx_from_coordinate(self, Coordinate xg not None):
        return c.get_idx_from_coordinate_fsel(self, xg)

    def coordinate_from_idx(self, long idx):
        return c.get_coordinate_from_idx_fsel(self, idx)

### -------------------------------------------------------------------

@q.timer
def mk_xg_field(Geometry geo):
    cdef FieldInt f = FieldInt()
    cc.set_xg_field(f.xx, geo.xx)
    return f

### -------------------------------------------------------------------
