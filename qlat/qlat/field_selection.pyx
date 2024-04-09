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

    def __init__(self, xg_arr=None, Geometry geo=None, cc.Bool distributed=False):
        """
        PointsSelection()
        PointsSelection(None, geo)
        PointsSelection(n_points, geo)
        PointsSelection(xg, geo)
        PointsSelection(xg_arr, geo)
        PointsSelection(xg_list, geo)
        PointsSelection(xg_arr, geo, distributed)
        """
        self.set_xg_arr(xg_arr, geo, distributed)

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

    def set_rand(self, Coordinate total_site not None, cc.Long n_points, RngState rs not None):
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

    def distributed(self):
        return self.xx.distributed

    def set_distributed(self, const cc.Bool distributed):
        """
        set the distributed flag
        """
        self.xx.distributed = distributed

    def n_points(self):
        return self.xx.size()

    def set_n_points(self, const cc.Long n_points):
        """
        will clear contents
        """
        self.xx.init(n_points)

    def xg_arr(self):
        """
        return xg for all selected points
        shape = (psel.n_points(), 4,)
        """
        return np.asarray(self, dtype=np.int32)

    @q.timer
    def set_xg_arr(self, xg_arr=None, Geometry geo=None, cc.Bool distributed=False):
        """
        psel.set_xg_arr()
        psel.set_xg_arr(None, geo)
        psel.set_xg_arr(n_points, geo)
        psel.set_xg_arr(xg, geo)
        psel.set_xg_arr(xg_arr, geo)
        psel.set_xg_arr(xg_list, geo)
        psel.set_xg_arr(xg_arr, geo, distributed)
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
            self.xx.init(n_points)
        elif isinstance(xg_arr, Coordinate):
            xg = xg_arr
            n_points = 1
            self.xx.init(n_points)
            self.xx[0] = xg.xx
        elif isinstance(xg_arr, np.ndarray):
            n_points = len(xg_arr)
            self.xx.init(n_points)
            assert xg_arr.shape == (n_points, 4,)
            np.asarray(self, dtype=np.int32)[:] = xg_arr
        elif isinstance(xg_arr, list):
            n_points = len(xg_arr)
            self.xx.init(n_points)
            for i in range(n_points):
                self.xx[i] = Coordinate(xg_arr[i]).xx
        else:
            raise Exception(f"PointsSelection.set_xg_arr({xg_arr},{geo})")
        self.xx.distributed = distributed

    def set_geo(self, Geometry geo):
        self.geo = geo

    def __setitem__(self, idx, val):
        """
        Implemented in terms of ``np.asarray``
        """
        np.asarray(self)[idx] = val

    def __getitem__(self, idx):
        """
        Implemented in terms of ``np.asarray``
        """
        return np.asarray(self)[idx]

    def __iter__(self):
        cdef cc.Long idx
        cdef cc.Long n_points = self.n_points()
        for idx in range(n_points):
            yield self.coordinate_from_idx(idx)

    def __len__(self):
        return self.n_points()

    def coordinate_from_idx(self, cc.Long idx):
        cdef Coordinate xg = Coordinate()
        xg.xx = self.xx[idx]
        return xg

    def intersect(self, FieldSelection fsel):
        """
        return new psel
        """
        cdef PointsSelection psel_new = self.copy()
        psel_new.xx = cc.intersect(fsel.xx, self.xx)
        return psel_new

    def crc32(self):
        return 0

    def __getstate__(self):
        """
        Only work when single node (or if all nodes has the same data).
        """
        geo = self.geo
        xg_arr = self.xg_arr()
        distributed = self.distributed()
        return [ xg_arr, geo, distributed, ]

    def __setstate__(self, state):
        """
        Only work when single node (or if all nodes has the same data).
        """
        if self.view_count > 0:
            raise ValueError("can't load while being viewed")
        self.__init__()
        [ xg_arr, geo, distributed, ] = state
        self.set_xg_arr(xg_arr, geo)
        if distributed:
            self.xx.distributed = True

### -------------------------------------------------------------------

cdef class FieldSelection:

    def __cinit__(self):
        self.cdata = <cc.Long>&(self.xx)
        self.view_count = 0

    def __init__(self, Geometry geo=None, cc.Long val=-1):
        if geo is not None:
            self.set_uniform(geo, val);

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        """
        Get buffer view of fsel.f_rank field as 1-D array of np.int64
        The values are rank (rank >= 0 means selected, rank == -1 means not selected)
        Need to call fsel.update() after modifying the f_rank via this buffer view.
        """
        cdef int ndim = 1
        cdef Buffer buf = Buffer(self, ndim)
        buf.format = 'q'
        buf.itemsize = sizeof(cc.Int64t)
        buf.buf = <char*>(self.xx.f_rank.field.data())
        cdef int multiplicity = self.xx.f_rank.get_geo().multiplicity
        assert multiplicity == 1
        buf.set_dim_size(0, self.xx.f_rank.field.size())
        buf.update_strides_from_shape()
        buf.set_buffer(buffer, flags)
        self.view_count += 1

    def release_buffer(self, Buffer buf):
        assert buf.obj is self
        self.view_count -= 1

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

    def update(self):
        """
        update various indices based on f_rank
        """
        cc.update_field_selection(self.xx)

    def set_empty(self, Geometry geo not None):
        """
        set an empty fsel with geo (all rank=-1)
        """
        if self.view_count > 0:
            raise Exception("FieldSelection: self.view_count>0")
        self.set_uniform(geo, -1)

    def set_uniform(self, Geometry geo not None, cc.Long val=0):
        """
        default (val = 0) select every sites
        val = -1 deselection everything
        """
        if self.view_count > 0:
            raise Exception("FieldSelection: self.view_count>0")
        self.xx.init()
        cc.mk_field_selection(self.xx.f_rank, geo.xx, val)
        self.update()

    def set_rand(self, Coordinate total_site not None, cc.Long n_per_tslice, RngState rs not None):
        if self.view_count > 0:
            raise Exception("FieldSelection: self.view_count>0")
        cc.mk_field_selection(self.xx.f_rank, total_site.xx, n_per_tslice, rs.xx)
        self.update()

    def set_rand_psel(self, Coordinate total_site not None, cc.Long n_per_tslice, RngState rs not None,
                      PointsSelection psel=None):
        if self.view_count > 0:
            raise Exception("FieldSelection: self.view_count>0")
        self.set_rand(total_site, n_per_tslice, rs)
        if psel is not None:
            self.add_psel(psel)
        self.update()

    def add_psel(self, PointsSelection psel, cc.Long rank_psel=1024 * 1024 * 1024 * 1024 * 1024):
        """
        Add psel points to the selection, with the rank specified as rank_psel.
        If the point is already selected with lower rank, the rank is unchanged.
        """
        cc.add_field_selection(self.xx.f_rank, psel.xx, rank_psel)
        self.update()

    def add_fsel(self, FieldSelection fsel):
        """
        Add fsel points to the selection, with the rank specified in fsel.
        If the point is already selected with lower rank, the rank is unchanged.
        """
        cc.add_field_selection(self.xx.f_rank, fsel.xx)
        self.update()

    def intersect_with(self, FieldSelection fsel):
        """
        Modify the `self`.
        More efficient if `self` is smaller than `fsel`.
        """
        cc.intersect_with(self.xx, fsel.xx)

    def intersect(self, FieldSelection fsel):
        """
        Do NOT change the `self`, but return a new one
        More efficient if `self` is smaller than `fsel`
        """
        cdef FieldSelection fsel_new = self.copy()
        fsel_new.intersect_with(fsel)
        return fsel_new

    def is_containing_psel(self, PointsSelection psel):
        cdef cc.Bool x = cc.is_containing(self.xx, psel.xx)
        return x

    def is_containing_fsel(self, FieldSelection fsel_small):
        cdef cc.Bool x = cc.is_containing(self.xx, fsel_small.xx)
        return x

    def is_containing(self, sel_small):
        if isinstance(sel_small, PointsSelection):
            return self.is_containing_psel(sel_small)
        elif isinstance(sel_small, FieldSelection):
            return self.is_containing_fsel(sel_small)
        else:
            raise Exception("'sel_small' not PointsSelection or FieldSelection")

    def to_psel(self):
        cdef PointsSelection psel = PointsSelection(None, self.geo())
        cc.assign_direct(psel.xx, cc.psel_from_fsel(self.xx))
        return psel

    def to_psel_local(self):
        cdef PointsSelection psel = PointsSelection(None, self.geo())
        cc.assign_direct(psel.xx, cc.psel_from_fsel_local(self.xx))
        return psel

    def save(self, const cc.std_string& path):
        cdef cc.Long total_bytes = cc.write_field_selection(self.xx, path)
        return total_bytes

    def load(self, const cc.std_string& path):
        if self.view_count > 0:
            raise Exception("FieldSelection: self.view_count>0")
        cdef cc.Long total_bytes = cc.read_field_selection(self.xx, path)
        return total_bytes

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

    def __setitem__(self, idx, val):
        """
        Implemented in terms of ``np.asarray``
        """
        np.asarray(self)[idx] = val

    def __getitem__(self, idx):
        """
        Implemented in terms of ``np.asarray``
        """
        return np.asarray(self)[idx]

    def __iter__(self):
        """
        iterate over all local selected coordinate as xg
        """
        cdef cc.Long idx
        cdef cc.Long n_elems = self.n_elems()
        for idx in range(n_elems):
            yield self.coordinate_from_idx(idx)

    def __len__(self):
        return self.n_elems()

    def idx_from_coordinate(self, Coordinate xg not None):
        cdef cc.Coordinate xl_xx = self.xx.f_local_idx.get_geo().coordinate_l_from_g(xg.xx)
        cdef cc.Long idx = self.xx.f_local_idx.get_elem(xl_xx)
        return idx

    def coordinate_from_idx(self, cc.Long idx):
        cdef Coordinate xg = Coordinate()
        cdef cc.Long index = self.xx.indices[idx]
        cdef cc.Coordinate xl_xx = self.xx.f_local_idx.get_geo().coordinate_from_index(index)
        cc.assign_direct(xg.xx, self.xx.f_local_idx.get_geo().coordinate_g_from_l(xl_xx))
        return xg

    def __getstate__(self):
        """
        Only work when single node.
        """
        geo = self.geo()
        fsel_arr = self[:].copy()
        return [ fsel_arr, geo, ]

    def __setstate__(self, state):
        """
        Only work when single node.
        """
        if self.view_count > 0:
            raise ValueError("can't load while being viewed")
        self.__init__()
        cdef Geometry geo
        [ fsel_arr, geo, ] = state
        self.set_empty(geo)
        self[:] = fsel_arr
        self.update()

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
