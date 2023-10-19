# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc
from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

import cqlat as c
import qlat_utils as q
import numpy as np

### -------------------------------------------------------------------

cdef class PointsSelection:

    def __cinit__(self):
        self.cdata = <long>&(self.xx)
        self.geo = None

    def __init__(self, coordinate_list=None, Geometry geo=None):
        cdef long n_points
        cdef long i
        self.geo = geo
        if coordinate_list is None:
            return
        if self.view_count > 0:
            raise ValueError("can't re-init while being viewed")
        n_points = len(coordinate_list)
        self.xx = cc.PointsSelection(n_points)
        for i in range(n_points):
            c = coordinate_list[i]
            if not isinstance(c, Coordinate):
                c = Coordinate(c)
            self.xx[i] = (<Coordinate>c).xx

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef long n_points = self.xx.size()
        cdef int ndim = 2
        cdef char* fmt = "i"
        cdef Buffer buf = Buffer(self, ndim, sizeof(cc.Int32t))
        cdef Py_ssize_t* shape = &buf.shape_strides[0]
        cdef Py_ssize_t* strides = &buf.shape_strides[buf.ndim]
        shape[0] = n_points
        shape[1] = 4
        buf.set_strides()
        buffer.buf = <char*>(self.xx.data())
        if flags & PyBUF_FORMAT:
            buffer.format = fmt
        else:
            buffer.format = NULL
        buffer.internal = NULL
        buffer.itemsize = buf.itemsize
        buffer.len = buf.get_len()
        buffer.ndim = buf.ndim
        buffer.obj = buf
        buffer.readonly = 0
        buffer.shape = shape
        buffer.strides = strides
        buffer.suboffsets = NULL
        self.view_count += 1

    def __releasebuffer__(self, Py_buffer *buffer):
        self.view_count -= 1

    def __imatmul__(self, PointsSelection v1):
        self.geo = v1.geo
        self.xx = v1.xx
        return self

    def copy(self):
        x = PointsSelection()
        x @= self
        return x

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        return self.copy()

    def set_rand(self, rs, total_site, n_points):
        if self.view_count > 0:
            raise ValueError("can't re-init while being viewed")
        c.set_rand_psel(self, rs, total_site, n_points)
        self.geo = Geometry(total_site)

    def save(self, path):
        c.save_psel(self, path)

    def load(self, path, geo):
        if self.view_count > 0:
            raise ValueError("can't re-init while being viewed")
        c.load_psel(self, path)
        self.geo = geo

    def n_points(self):
        return c.get_n_points_psel(self)

    def to_list(self):
        return c.mk_list_psel(self)

    def from_list(self, coordinate_list, geo=None):
        if self.view_count > 0:
            raise ValueError("can't re-init while being viewed")
        c.set_list_psel(self, coordinate_list)
        self.geo = geo
        return self

    def coordinate_from_idx(self, idx):
        return c.get_coordinate_from_idx_psel(self, idx)

### -------------------------------------------------------------------

cdef class FieldSelection:

    def __cinit__(self):
        self.cdata = <long>&(self.xx)

    def __init__(self, total_site=None, n_per_tslice=-1, rs=None, psel=None):
        if total_site is not None:
            assert isinstance(rs, RngState)
            assert isinstance(n_per_tslice, int)
            c.set_rand_fsel(self, rs, total_site, n_per_tslice)
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

    def set_uniform(self, total_site, val=0):
        """
        default (val = 0) select every sites
        val = -1 deselection everything
        """
        c.set_uniform_fsel(self, total_site, val)

    def set_rand(self, rs, total_site, n_per_tslice):
        assert isinstance(rs, RngState)
        assert isinstance(n_per_tslice, int)
        c.set_rand_fsel(self, rs, total_site, n_per_tslice)
        self.update()
        self.update(n_per_tslice)

    def add_psel(self, psel, rank_psel=1024 * 1024 * 1024 * 1024 * 1024):
        """
        Add psel points to the selection, with the rank specified as rank_psel.
        If the point is already selected with lower rank, the rank is unchanged.
        """
        c.add_psel_fsel(self, psel, rank_psel)
        self.update()

    def update(self, n_per_tslice=-1):
        """
        if n_per_tslice < 0: only update various indices
        if n_per_tslice >= 0: only update parameters (n_per_tslice and prob)
        """
        c.update_fsel(self, n_per_tslice)

    def select_rank_range(self, rank_start=0, rank_stop=-1):
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

    def select_t_range(self, rank_start=0, rank_stop=-1):
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

    def save(self, path):
        return c.save_fsel(self, path)

    def load(self, path, n_per_tslice):
        return c.load_fsel(self, path, n_per_tslice)

    def geo(self):
        geo = Geometry((0, 0, 0, 0))
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

    def idx_from_coordinate(self, xg):
        return c.get_idx_from_coordinate_fsel(self, xg)

    def coordinate_from_idx(self, idx):
        return c.get_coordinate_from_idx_fsel(self, idx)

### -------------------------------------------------------------------
