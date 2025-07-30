# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc
from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

from .geometry cimport *
from .field_types cimport *
from .field_base cimport (
        SelectedPointsBase,
        SelectedFieldBase,
        FieldBase,
        )
from .selected_points_types cimport (
        SelectedPointsChar,
        )

import cqlat as c
import numpy as np

class q:
    from qlat_utils import (
            timer,
            get_id_node,
            get_num_node,
            hash_sha256,
            mk_cache,
            get_fname,
            displayln_info,
            )
    from .mpi_utils import (
            get_comm,
            )

### -------------------------------------------------------------------

cdef class SelectedShufflePlan:

    def __init__(self, *args):
        """
        SelectedShufflePlan()
        SelectedShufflePlan("r_from_l", psel_src, geo, rs)
        SelectedShufflePlan("r_from_l", psel_src_list, geo_src_list, rs)
        SelectedShufflePlan("t_slice_from_l", psel_src_list, geo_src_list)
        SelectedShufflePlan("dist_t_slice_from_l", psel_src, geo, num_field)
        """
        self.psel_src_list = None
        self.psel_dst_list = None
        self.fsel_src_list = None
        self.fsel_dst_list = None
        self.geo_src_list = None
        self.geo_dst_list = None
        self.xx.init()
        if len(args) == 0:
            return
        elif len(args) == 1:
            return
        elif args[0] == "r_from_l" and isinstance(args[1], PointsSelection):
            self.init_from_psel_r_from_l(*args[1:])
        elif args[0] == "r_from_l" and isinstance(args[1], list) and len(args[1]) > 0 and isinstance(args[1][0], PointsSelection):
            self.init_from_psel_list_r_from_l(*args[1:])
        elif args[0] == "t_slice_from_l" and isinstance(args[1], list) and len(args[1]) > 0 and isinstance(args[1][0], PointsSelection):
            self.init_from_psel_list_t_slice_from_l(*args[1:])
        elif args[0] == "dist_t_slice_from_l" and isinstance(args[1], PointsSelection):
            self.init_from_psel_list_dist_t_slice_from_l(*args[1:])
        else:
            raise Exception(f"SelectedShufflePlan.__init__ {args}")

    @q.timer
    def init_from_psel_r_from_l(self, PointsSelection psel_src, Geometry geo_src, RngState rs):
        """
        shuffle to PointsDistType::Random ("r") from PointsDistType::Local ("l").
        """
        cdef PointsSelection psel_dst = PointsSelection()
        self.xx.init()
        cc.set_selected_shuffle_plan_r_from_l(self.xx, psel_src.xx, geo_src.xx, rs.xx)
        psel_dst.xx.points_dist_type = self.xx.points_dist_type_recv
        cc.shuffle_points_selection(psel_dst.xx, psel_src.xx, self.xx)
        self.psel_src_list = [ psel_src, ]
        self.psel_dst_list = [ psel_dst, ]

    @q.timer
    def init_from_psel_list_r_from_l(self, list psel_src_list, list geo_src_list, RngState rs):
        """
        shuffle to PointsDistType::Random ("r") from PointsDistType::Local ("l").
        """
        assert len(geo_src_list) == len(psel_src_list)
        cdef cc.std_vector[cc.PointsSelection] psel_src_vec
        psel_src_vec.resize(len(psel_src_list))
        cdef PointsSelection psel
        for i in range(psel_src_vec.size()):
            psel = psel_src_list[i]
            cc.qswap(psel.xx, psel_src_vec[i])
        cdef cc.std_vector[cc.Geometry] geo_src_vec
        geo_src_vec.resize(len(geo_src_list))
        cdef Geometry geo
        for i in range(geo_src_vec.size()):
            geo = geo_src_list[i]
            geo_src_vec[i] = geo.xx
        self.xx.init()
        cc.set_selected_shuffle_plan_r_from_l(self.xx, psel_src_vec, geo_src_vec, rs.xx)
        cdef cc.std_vector[cc.PointsSelection] psel_dst_vec
        cc.shuffle_points_selection(psel_dst_vec, psel_src_vec, self.xx)
        psel_dst_list = [ PointsSelection() for i in range(psel_dst_vec.size()) ]
        for i in range(psel_src_vec.size()):
            psel = psel_src_list[i]
            cc.qswap(psel.xx, psel_src_vec[i])
        for i in range(psel_dst_vec.size()):
            psel = psel_dst_list[i]
            cc.qswap(psel.xx, psel_dst_vec[i])
        self.psel_src_list = psel_src_list
        self.psel_dst_list = psel_dst_list

    @q.timer
    def init_from_psel_list_t_slice_from_l(self, list psel_src_list, list geo_src_list):
        """
        shuffle to PointsDistType::Local ("l") from PointsDistType::Local ("l").
        """
        assert len(geo_src_list) == len(psel_src_list)
        cdef cc.std_vector[cc.PointsSelection] psel_src_vec
        cdef cc.std_vector[cc.PointsSelection] psel_dst_vec
        cdef PointsSelection psel
        psel_src_vec.resize(len(psel_src_list))
        for i in range(psel_src_vec.size()):
            psel = psel_src_list[i]
            cc.qswap(psel.xx, psel_src_vec[i])
        cdef cc.std_vector[cc.Geometry] geo_src_vec
        geo_src_vec.resize(len(geo_src_list))
        cdef Geometry geo
        for i in range(geo_src_vec.size()):
            geo = geo_src_list[i]
            geo_src_vec[i] = geo.xx
        self.xx.init()
        cc.set_selected_shuffle_plan_t_slice_from_l(self.xx, psel_src_vec, geo_src_vec)
        cc.shuffle_points_selection(psel_dst_vec, psel_src_vec, self.xx)
        psel_dst_list = [ PointsSelection() for i in range(psel_dst_vec.size()) ]
        for i in range(psel_src_vec.size()):
            psel = psel_src_list[i]
            cc.qswap(psel.xx, psel_src_vec[i])
        for i in range(psel_dst_vec.size()):
            psel = psel_dst_list[i]
            cc.qswap(psel.xx, psel_dst_vec[i])
        self.psel_src_list = psel_src_list
        self.psel_dst_list = psel_dst_list

    @q.timer
    def init_from_psel_list_dist_t_slice_from_l(self, PointsSelection psel_src, Geometry geo, int num_field):
        """
        shuffle to PointsDistType::Local ("l") from PointsDistType::Local ("l").
        match `init_from_psel_list_t_slice_from_l`.
        Example use case:
            Perform spatial smearing for propagators:
            (1) shuffle propagators with "t_slice_from_l".
            (2) shuffle gauge field with "dist_t_slice_from_l".
        """
        cdef cc.std_vector[cc.PointsSelection] psel_src_vec
        cdef cc.std_vector[cc.PointsSelection] psel_dst_vec
        cdef PointsSelection psel
        psel_src_vec.resize(1)
        cc.qswap(psel_src.xx, psel_src_vec[0])
        self.xx.init()
        cc.set_selected_shuffle_plan_dist_t_slice_from_l(self.xx, psel_src_vec[0], geo.xx, num_field)
        cc.shuffle_points_selection(psel_dst_vec, psel_src_vec, self.xx)
        psel_dst_list = [ PointsSelection() for i in range(psel_dst_vec.size()) ]
        cc.qswap(psel_src.xx, psel_src_vec[0])
        for i in range(psel_dst_vec.size()):
            psel = psel_dst_list[i]
            cc.qswap(psel.xx, psel_dst_vec[i])
        self.psel_src_list = [ psel_src, ]
        self.psel_dst_list = psel_dst_list

    @property
    def points_dist_type_send(self):
        return cc.show(self.xx.points_dist_type_send)

    @property
    def points_dist_type_recv(self):
        return cc.show(self.xx.points_dist_type_recv)

    @property
    def num_selected_points_send(self):
        return self.xx.num_selected_points_send

    @property
    def num_selected_points_recv(self):
        return self.xx.num_selected_points_recv

    @property
    def psel_send_list(self):
        return self.psel_src_list

    @property
    def psel_recv_list(self):
        return self.psel_dst_list

    @property
    def geo_send_list(self):
        if self.geo_src_list is not None:
            return self.geo_src_list
        self.geo_src_list = []
        assert self.xx.size_node_send.size() == self.xx.coor_node_send.size()
        cdef cc.Int num_geo = self.xx.size_node_send.size()
        cdef Coordinate total_site = Coordinate()
        total_site.xx = self.xx.total_site
        cdef Coordinate size_node = Coordinate()
        cdef Coordinate coor_node = Coordinate()
        cdef Coordinate node_site = Coordinate()
        cdef Geometry geo
        if total_site == Coordinate():
            return self.geo_src_list
        cdef list geo_list = []
        for i in range(num_geo):
            size_node.xx = self.xx.size_node_send[i]
            coor_node.xx = self.xx.coor_node_send[i]
            node_site.xx = total_site.xx / size_node.xx
            assert node_site * size_node == total_site
            geo = Geometry(coor_node, size_node, node_site)
            geo_list.append(geo)
        self.geo_src_list = geo_list
        return self.geo_src_list

    @property
    def geo_recv_list(self):
        if self.geo_dst_list is not None:
            return self.geo_dst_list
        self.geo_dst_list = []
        assert self.xx.size_node_recv.size() == self.xx.coor_node_recv.size()
        cdef cc.Int num_geo = self.xx.size_node_recv.size()
        cdef Coordinate total_site = Coordinate()
        total_site.xx = self.xx.total_site
        cdef Coordinate size_node = Coordinate()
        cdef Coordinate coor_node = Coordinate()
        cdef Coordinate node_site = Coordinate()
        cdef Geometry geo
        if total_site == Coordinate():
            return self.geo_dst_list
        cdef list geo_list = []
        for i in range(num_geo):
            size_node.xx = self.xx.size_node_recv[i]
            coor_node.xx = self.xx.coor_node_recv[i]
            node_site.xx = total_site.xx / size_node.xx
            assert node_site * size_node == total_site
            geo = Geometry(coor_node, size_node, node_site)
            geo_list.append(geo)
        self.geo_dst_list = geo_list
        return self.geo_dst_list

    @property
    def fsel_send_list(self):
        if self.fsel_src_list is None:
            return self.fsel_src_list
        self.fsel_src_list = []
        if self.points_dist_type_send not in [ "l", "f", ]:
            return self.fsel_src_list
        if len(self.psel_send_list) != len(self.geo_send_list):
            return self.fsel_src_list
        cdef FieldSelection fsel
        cdef list fsel_list = []
        for psel, geo in zip(self.psel_send_list, self.geo_send_list):
            fsel = FieldSelection(psel, geo)
            fsel_list.append(fsel)
        self.fsel_src_list = fsel_list
        return self.fsel_src_list

    @property
    def fsel_recv_list(self):
        if self.fsel_dst_list is None:
            return self.fsel_dst_list
        self.fsel_dst_list = []
        if self.points_dist_type_recv not in [ "l", "f", ]:
            return self.fsel_dst_list
        if len(self.psel_recv_list) != len(self.geo_recv_list):
            return self.fsel_dst_list
        cdef FieldSelection fsel
        cdef list fsel_list = []
        for psel, geo in zip(self.psel_recv_list, self.geo_recv_list):
            fsel = FieldSelection(psel, geo)
            fsel_list.append(fsel)
        self.fsel_dst_list = fsel_list
        return self.fsel_dst_list

    @q.timer
    def shuffle(self, SelectedPointsChar sp_src, *, bint is_reverse=False):
        """
        shuffle `sp_src` with this plan
        return `sp_dst`
        The type for all `sp` is `SelectedPointsChar`
        """
        assert 1 == len(self.psel_send_list)
        assert 1 == len(self.psel_recv_list)
        assert isinstance(sp_src, SelectedPointsChar)
        psel_src = self.psel_send_list[0]
        psel_dst = self.psel_recv_list[0]
        psel = sp_src.psel
        cdef SelectedPointsChar sp_dst
        if is_reverse:
            assert (psel is None) or (psel == psel_dst)
            sp_dst = SelectedPointsChar(psel_src)
            cc.shuffle_selected_points_back(sp_dst.xx, sp_src.xx, self.xx)
        else:
            assert (psel is None) or (psel == psel_src)
            sp_dst = SelectedPointsChar(psel_dst)
            cc.shuffle_selected_points(sp_dst.xx, sp_src.xx, self.xx)
        return sp_dst

    @q.timer
    def shuffle_list(self, list sp_src_list, *, bint is_reverse=False):
        """
        shuffle `sp_src_list` with this plan
        return `sp_dst_list`
        The type for all `sp` is `SelectedPointsChar`
        """
        if is_reverse:
            psel_src_list = self.psel_recv_list
            psel_dst_list = self.psel_send_list
        else:
            psel_src_list = self.psel_send_list
            psel_dst_list = self.psel_recv_list
        assert len(sp_src_list) == len(psel_src_list)
        for sp_src, psel_src in zip(sp_src_list, psel_src_list):
            assert isinstance(sp_src, SelectedPointsChar)
            psel = sp_src.psel
            assert (psel is None) or (psel == psel_src)
        cdef cc.std_vector[cc.SelectedPoints[cc.Char]] sp_src_vec
        cdef cc.std_vector[cc.SelectedPoints[cc.Char]] sp_dst_vec
        cdef SelectedPointsChar sp
        sp_src_vec.resize(len(psel_src_list))
        for i in range(sp_src_vec.size()):
            sp = sp_src_list[i]
            cc.qswap(sp.xx, sp_src_vec[i])
        if is_reverse:
            cc.shuffle_selected_points_back(sp_dst_vec, sp_src_vec, self.xx)
        else:
            cc.shuffle_selected_points(sp_dst_vec, sp_src_vec, self.xx)
        assert <cc.Long>sp_dst_vec.size() == len(psel_dst_list)
        sp_dst_list = [ SelectedPointsChar(psel) for psel in psel_dst_list ]
        for i in range(sp_src_vec.size()):
            sp = sp_src_list[i]
            cc.qswap(sp.xx, sp_src_vec[i])
        for i in range(sp_dst_vec.size()):
            sp = sp_dst_list[i]
            assert sp.psel is psel_dst_list[i]
            cc.qswap(sp.xx, sp_dst_vec[i])
        return sp_dst_list

    @q.timer
    def shuffle_sp(self, cls, object src, *, bint is_reverse=False):
        """
        shuffle `src` with this plan
        return `sp_dst`
        The type for the returned object is `cls`
        """
        assert len(self.psel_send_list) == 1
        assert len(self.psel_recv_list) == 1
        cdef SelectedPointsChar spc_src
        cdef SelectedPointsChar spc_dst
        cdef Geometry geo
        cdef list geo_send_list
        cdef list geo_recv_list
        if is_reverse:
            geo_send_list = self.geo_recv_list
            geo_recv_list = self.geo_send_list
        else:
            geo_send_list = self.geo_send_list
            geo_recv_list = self.geo_recv_list
        cdef Geometry geo_send = None
        cdef Geometry geo_recv = None
        if geo_send_list is not None:
            if len(geo_send_list) == 1:
                geo_send = geo_send_list[0]
        if geo_recv_list is not None:
            if len(geo_recv_list) == 1:
                geo_recv = geo_recv_list[0]
        spc_src = SelectedPointsChar()
        if isinstance(src, SelectedPointsBase):
            src.swap_cast(spc_src)
            spc_dst = self.shuffle(spc_src, is_reverse=is_reverse)
            src.swap_cast(spc_src)
        elif isinstance(src, SelectedFieldBase):
            geo = Geometry()
            src.swap_sp_cast(spc_src, geo)
            if geo_send is not None:
                assert geo == geo_send
            spc_dst = self.shuffle(spc_src, is_reverse=is_reverse)
            src.swap_sp_cast(spc_src, geo)
        elif isinstance(src, FieldBase):
            geo = Geometry()
            src.swap_sp_cast(spc_src, geo)
            if geo_send is not None:
                assert geo == geo_send
            spc_dst = self.shuffle(spc_src, is_reverse=is_reverse)
            src.swap_sp_cast(spc_src, geo)
        else:
            assert False
        dst = cls()
        if geo_recv is None:
            geo_recv = Geometry()
        geo = geo_recv.copy()
        if isinstance(src, SelectedPointsBase):
            dst.swap_cast(spc_dst)
            if is_reverse:
                psel_list = self.psel_send_list
            else:
                psel_list = self.psel_recv_list
            assert len(psel_list) == 1
            dst.psel = psel_list[0]
        elif isinstance(src, SelectedFieldBase):
            dst.swap_sp_cast(spc_dst, geo)
            if is_reverse:
                fsel_list = self.fsel_send_list
            else:
                fsel_list = self.fsel_recv_list
            if len(fsel_list) == 1:
                dst.fsel = fsel_list[0]
        elif isinstance(src, FieldBase):
            dst.swap_sp_cast(spc_dst, geo)
        else:
            assert False
        return dst

    @q.timer
    def shuffle_sp_list(self, cls, list src_list, *, bint is_reverse=False):
        """
        shuffle `src_list` with this plan
        return `sp_dst_list`
        The type for the returned list of objects are `cls`
        """
        cdef list geo_send_list
        cdef list geo_recv_list
        cdef list psel_send_list
        cdef list psel_recv_list
        cdef int num_send
        cdef int num_recv
        if is_reverse:
            geo_send_list = self.geo_recv_list
            geo_recv_list = self.geo_send_list
            psel_send_list = self.psel_recv_list
            psel_recv_list = self.psel_send_list
            num_send = self.num_selected_points_recv
            num_recv = self.num_selected_points_send
        else:
            geo_send_list = self.geo_send_list
            geo_recv_list = self.geo_recv_list
            psel_send_list = self.psel_send_list
            psel_recv_list = self.psel_recv_list
            num_send = self.num_selected_points_send
            num_recv = self.num_selected_points_recv
        assert len(psel_send_list) == num_send
        assert len(psel_recv_list) == num_recv
        assert len(src_list) == num_send
        if len(geo_send_list) != 0:
            assert len(geo_send_list) == num_send
        else:
            geo_send_list = [ None for i in range(num_send) ]
        if len(geo_recv_list) != 0:
            assert len(geo_recv_list) == num_recv
        else:
            geo_recv_list = [ None for i in range(num_recv) ]
        spc_src_list = []
        geo_src_list = []
        for idx, (src, geo_send,) in enumerate(zip(src_list, geo_send_list)):
            spc_src = SelectedPointsChar()
            geo_src = Geometry()
            if isinstance(src, SelectedPointsBase):
                src.swap_cast(spc_src)
                if src.psel is not None:
                    assert src.psel == psel_send_list[idx]
            elif isinstance(src, SelectedFieldBase):
                src.swap_sp_cast(spc_src, geo_src)
                assert geo_send is not None
                assert geo_src == geo_send
                if is_reverse:
                    fsel_list = self.fsel_recv_list
                else:
                    fsel_list = self.fsel_send_list
                assert len(fsel_list) == num_send
                if src.fsel is not None:
                    if src.fsel is not fsel_list[idx]:
                        assert PointsSelection(src.fsel) == self.psel_send_list[idx]
            elif isinstance(src, FieldBase):
                src.swap_sp_cast(spc_src, geo_src)
                assert geo_send is not None
                assert geo_src == geo_send
                assert spc_src.n_points == geo_send.local_volume_expanded
            else:
                assert False
            spc_src_list.append(spc_src)
            geo_src_list.append(geo_src)
        assert len(spc_src_list) == num_send
        spc_dst_list = self.shuffle_list(spc_src_list, is_reverse=is_reverse)
        assert len(spc_dst_list) == num_recv
        for src, spc_src, geo_src in zip(src_list, spc_src_list, geo_src_list):
            if isinstance(src, SelectedPointsBase):
                src.swap_cast(spc_src)
            elif isinstance(src, SelectedFieldBase):
                src.swap_sp_cast(spc_src, geo_src)
            elif isinstance(src, FieldBase):
                src.swap_sp_cast(spc_src, geo_src)
            else:
                assert False
        dst_list = []
        cdef Geometry geo
        for idx, (spc_dst, geo_recv,) in enumerate(zip(spc_dst_list, geo_recv_list)):
            dst = cls()
            if isinstance(src, SelectedPointsBase):
                dst.swap_cast(spc_dst)
                assert dst.psel is psel_recv_list[idx]
            elif isinstance(src, SelectedFieldBase):
                assert geo_recv is not None
                geo = geo_recv.copy()
                dst.swap_sp_cast(spc_dst, geo)
                if is_reverse:
                    fsel_list = self.fsel_send_list
                else:
                    fsel_list = self.fsel_recv_list
                assert len(fsel_list) == num_recv
                dst.fsel = fsel_list[idx]
            elif isinstance(src, FieldBase):
                assert geo_recv is not None
                geo = geo_recv.copy()
                assert geo.local_volume_expanded == spc_dst.n_points
                dst.swap_sp_cast(spc_dst, geo)
            else:
                assert False
            dst_list.append(dst)
        assert len(dst_list) == num_recv
        return dst_list

###

cdef class PointsSelection:

    def __cinit__(self):
        self.cdata = <cc.Long>&(self.xx)
        self.view_count = 0

    def __init__(self, *args):
        """
        PointsSelection()
        PointsSelection(total_site)
        PointsSelection(total_site, xg_arr)
        PointsSelection(total_site, xg_list)
        PointsSelection(total_site, xg)
        PointsSelection(total_site, n_points)
        PointsSelection(total_site, xg_arr, points_dist_type)
        PointsSelection(geo) # full selection (type "f")
        PointsSelection(fsel) # (type "l")
        PointsSelection(fsel, ssp)
        PointsSelection(psel) # A copy of `psel`.
        PointsSelection(psel, ssp)
        PointsSelection(psel, ssp, is_reverse)
        #
        points_dist_type in [ "g", "f", "l", "r", "o", ]
        """
        cdef cc.Int len_args = len(args)
        if len_args == 0:
            self.xx.init()
        elif isinstance(args[0], Coordinate):
            total_site = args[0]
            self.init_from_total_site(*args)
        elif isinstance(args[0], Geometry):
            self.init_from_geo(*args)
        elif isinstance(args[0], FieldSelection):
            self.init_from_fsel(*args)
        elif isinstance(args[0], PointsSelection):
            self.init_from_psel(*args)
        else:
            raise Exception(f"PointsSelection::__init__: {args}")

    def init_from_total_site(self, Coordinate total_site, object xg_arr=None, str points_dist_type=None):
        """
        xg_arr can be n_points, xg, xg_arr, xg_list.
        points_dist_type in [ "g", "f", "l", "r", "o", ]
        """
        self.xx.init()
        self.xx.init(total_site.xx, 0)
        if xg_arr is None:
            return
        self.xg_arr = xg_arr
        if points_dist_type is None:
            return
        self.points_dist_type = points_dist_type

    def init_from_geo(self, Geometry geo):
        """
        points_dist_type in [ "f", ]
        """
        cc.set_psel_full(self.xx, geo.xx)

    def init_from_psel(self, PointsSelection psel, SelectedShufflePlan ssp=None, bint is_reverse=False):
        """
        Shuffle according to `ssp`.
        """
        fname = q.get_fname()
        self.xx.init()
        if ssp is None:
            self @= psel
        elif (len(ssp.psel_send_list) == 1 and ssp.psel_send_list[0] is psel) and (not is_reverse):
            assert len(ssp.psel_recv_list) == 1
            self @= ssp.psel_recv_list[0]
        elif (len(ssp.psel_recv_list) == 1 and ssp.psel_recv_list[0] is psel) and is_reverse:
            assert len(ssp.psel_send_list) == 1
            self @= ssp.psel_send_list[0]
        else:
            q.displayln_info(f"WARNING: {fname}: psel is not ssp.psel_src")
            assert cc.show(psel.xx.points_dist_type) == cc.show(ssp.xx.points_dist_type_send)
            self.xx.points_dist_type = ssp.xx.points_dist_type_recv
            if is_reverse:
                cc.shuffle_points_selection_back(self.xx, psel.xx, ssp.xx)
            else:
                cc.shuffle_points_selection(self.xx, psel.xx, ssp.xx)

    def init_from_fsel(self, FieldSelection fsel, SelectedShufflePlan ssp=None):
        """
        Shuffle according to `ssp`.
        self.points_dist_type == "l" for PointsDistType::Local (if ssp is None)
        """
        self.xx.init()
        if ssp is None:
            cc.set_psel_from_fsel(self.xx, fsel.xx)
        else:
            self.xx.points_dist_type = ssp.xx.points_dist_type_recv
            cc.shuffle_field_selection(self.xx, fsel.xx, ssp.xx)

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

    def swap(self, PointsSelection other):
        cc.qswap(self.xx, other.xx)

    def swap_sp_cast(self, SelectedPointsChar other, Coordinate total_site):
        assert other.psel is None
        cc.qswap_cast(self.xx, other.xx, total_site.xx)

    @property
    def points_dist_type(self):
        """
        points_dist_type in [ "g", "f", "l", "r", "o", ]
        Meaning: Global, Full, Local, Random, Other
        """
        return cc.show(self.xx.points_dist_type)

    @points_dist_type.setter
    def points_dist_type(self, str value):
        """
        set the points_dist_type flag
        value in [ "g", "f", "l", "r", "o", ]
        """
        self.xx.points_dist_type = cc.read_points_dist_type(value)

    @property
    def total_site(self):
        cdef Coordinate value = Coordinate()
        value.xx = self.xx.total_site
        return value

    @total_site.setter
    def total_site(self, Coordinate value):
        self.xx.total_site = value.xx

    @property
    def geo(self):
        cdef Geometry geo = Geometry(self.total_site)
        return geo

    @property
    def n_points(self):
        return self.xx.size()

    @property
    def xg_arr(self):
        """
        return xg for all selected points
        shape = (psel.n_points, 4,)
        """
        return np.asarray(self, dtype=np.int32)

    @xg_arr.setter
    def xg_arr(self, object xg_arr):
        """
        xg_arr can be n_points, xg, xg_arr, xg_list.
        """
        if self.view_count > 0:
            raise ValueError("can't re-init while being viewed")
        cdef cc.Coordinate total_site = self.xx.total_site
        cdef cc.PointsDistType points_dist_type = self.xx.points_dist_type
        cdef cc.Long n_points
        cdef cc.Long i
        cdef Coordinate xg
        if isinstance(xg_arr, int):
            n_points = xg_arr
            self.xx.init(total_site, n_points, points_dist_type)
        elif isinstance(xg_arr, Coordinate):
            xg = xg_arr
            n_points = 1
            self.xx.init(total_site, n_points, points_dist_type)
            self.xx[0] = xg.xx
        elif isinstance(xg_arr, np.ndarray):
            n_points = len(xg_arr)
            assert xg_arr.shape == (n_points, 4,)
            self.xx.init(total_site, n_points, points_dist_type)
            np.asarray(self, dtype=np.int32)[:] = xg_arr
        elif isinstance(xg_arr, list):
            n_points = len(xg_arr)
            self.xx.init(total_site, n_points, points_dist_type)
            for i in range(n_points):
                self.xx[i] = Coordinate(xg_arr[i]).xx
        else:
            raise Exception(f"PointsSelection.xg_arr = {xg_arr}")

    def set_rand(self, Coordinate total_site not None, cc.Long n_points, RngState rs not None):
        if self.view_count > 0:
            raise ValueError("can't re-init while being viewed")
        cc.assign_direct(self.xx, cc.mk_random_points_selection(total_site.xx, n_points, rs.xx))

    def to_lat_data(self):
        cdef LatDataInt ld = LatDataInt()
        cc.lat_data_from_points_selection(ld.xx, self.xx)
        return ld

    def from_lat_data(self, LatDataInt ld not None):
        cc.points_selection_from_lat_data(self.xx, ld.xx)

    @q.timer
    def bcast(self, cc.Int root=0):
        cdef LatDataInt ld
        if cc.get_num_node() != 1:
            if cc.get_id_node() == root:
                ld = self.to_lat_data()
            else:
                ld = LatDataInt()
            ld.bcast(root)
            if cc.get_id_node() != root:
                self.from_lat_data(ld)
        return self

    def save_str(self):
        """
        only return str at node 0
        """
        cdef LatDataInt ld
        if q.get_id_node() == 0:
            ld = self.to_lat_data()
            return ld.save_str()
        else:
            return bytes()

    def load_str(self, bytes content):
        """
        only need str at node 0
        """
        cdef LatDataInt ld = LatDataInt()
        if q.get_id_node() == 0:
            ld.load_str(content)
        ld.bcast()
        self.from_lat_data(ld)

    def save(self, const cc.std_string& path, *, is_sync_node=True):
        if is_sync_node:
            cc.save_points_selection_info(self.xx, path)
        else:
            cc.save_points_selection(self.xx, path)

    def load(self, const cc.std_string& path, Geometry geo=None, *, is_sync_node=True):
        if self.view_count > 0:
            raise ValueError("can't re-init while being viewed")
        if is_sync_node:
            cc.assign_direct(self.xx, cc.load_points_selection_info(path))
        else:
            cc.assign_direct(self.xx, cc.load_points_selection(path))
        cdef Coordinate total_site
        if path.endswith(".lati"):
            if geo is not None:
                total_site = geo.total_site
                assert self.xx.total_site == total_site.xx
        else:
            assert self.xx.total_site == cc.Coordinate()
            assert geo is not None
            total_site = geo.total_site
            self.xx.total_site = total_site.xx

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
        cdef cc.Long n_points = self.n_points
        for idx in range(n_points):
            yield self.coordinate_from_idx(idx)

    def __len__(self):
        return self.xx.size()

    def __eq__(self, PointsSelection other):
        if self is other:
            return True
        if self.xx == other.xx:
            return True
        return False

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

    def is_containing_psel(self, PointsSelection psel_small):
        cdef cc.Bool x = cc.is_containing(self.xx, psel_small.xx)
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
            raise Exception("PointsSelection: 'sel_small' not PointsSelection or FieldSelection sel_small={sel_small}")

    def __repr__(self):
        """
        only show information of this node.
        """
        return f"PointsSelection({self.points_dist_type}, {self.total_site}, {self.xg_arr.tolist()})"

    def hash_sha256(self):
        """
        Return hash of psel.
        Always return the same hash from all the node.
        If points_dist_type == "g", then return the hash for each node, and verify they are all the same.
        Otherwise, return the hash of the gathered information from all the node.
        """
        v = q.hash_sha256((
            "PointsSelection:",
            self.xg_arr,
            self.total_site,
            self.points_dist_type,
            ))
        v_list = q.get_comm().allgather(v)
        if self.points_dist_type == "g":
            for v1 in v_list:
                assert v == v1
            return v
        return q.hash_sha256(v_list)

    def __getstate__(self):
        """
        Only work if run with single node (or if all nodes has the same data).
        """
        xg_arr = self.xg_arr
        total_site = self.total_site
        points_dist_type = self.points_dist_type
        return [ xg_arr, total_site, points_dist_type, ]

    def __setstate__(self, state):
        """
        Only work if run with single node (or if all nodes has the same data).
        """
        if self.view_count > 0:
            raise ValueError("can't load while being viewed")
        [ xg_arr, total_site, points_dist_type, ] = state
        self.__init__(total_site, xg_arr, points_dist_type)

### -------------------------------------------------------------------

cdef class FieldSelection:

    def __cinit__(self):
        self.cdata = <cc.Long>&(self.xx)
        self.view_count = 0

    def __init__(self, *args):
        """
        FieldSelection()
        FieldSelection(geo, 0) # selecting all points
        FieldSelection(geo, -1) # no points being selected
        FieldSelection(geo) # same as `FieldSelection(geo, -1)`, no points being selected
        FieldSelection(psel) # require `psel.points_dist_type in [ "l", "f", "g", ]`
        FieldSelection(psel, geo) # require `psel.points_dist_type in [ "l", "f", "g", ]`
        """
        cdef cc.Int len_args = len(args)
        if len_args == 0:
            self.xx.init()
            return
        elif isinstance(args[0], Geometry):
            self.init_from_geo(*args)
        elif isinstance(args[0], PointsSelection):
            self.init_from_psel(*args)
        else:
            raise Exception("FieldSelection.__init__: {args}")

    def init_from_geo(self, Geometry geo, cc.Long rank=-1):
        """
        By default `rank=-1` means empty selection.
        Non-negative `rank` means full selection.
        """
        self.set_uniform(geo, rank)

    def init_from_psel(self, PointsSelection psel, Geometry geo=None, cc.Long rank_psel=1024 * 1024 * 1024 * 1024 * 1024):
        """
        psel.points_dist_type [ "g", "f", "l", ]
        """
        if geo is None:
            geo = psel.geo
        cc.set_fsel_from_psel(self.xx, psel.xx, geo.xx, rank_psel)

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
        cdef int multiplicity = self.xx.f_rank.multiplicity
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

    def is_containing_psel(self, PointsSelection psel_small):
        cdef cc.Bool x = cc.is_containing(self.xx, psel_small.xx)
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
            raise Exception("PointsSelection: 'sel_small' not PointsSelection or FieldSelection sel_small={sel_small}")

    def to_psel(self):
        cdef PointsSelection psel = PointsSelection(self.total_site)
        cc.assign_direct(psel.xx, cc.psel_from_fsel(self.xx))
        return psel

    def to_psel_local(self):
        return PointsSelection(self)

    def save(self, const cc.std_string& path):
        cdef cc.Long total_bytes = cc.write_field_selection(self.xx, path)
        return total_bytes

    def load(self, const cc.std_string& path):
        if self.view_count > 0:
            raise Exception("FieldSelection: self.view_count>0")
        cdef cc.Long total_bytes = cc.read_field_selection(self.xx, path)
        return total_bytes

    @property
    def geo(self):
        cdef Geometry geo = Geometry()
        geo.xx = self.xx.get_geo()
        return geo

    @property
    def total_site(self):
        cdef Coordinate total_site = Coordinate()
        cc.assign_direct(total_site.xx, self.xx.get_geo().total_site())
        return total_site

    @property
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
        cdef cc.Long n_elems = self.n_elems
        for idx in range(n_elems):
            yield self.coordinate_from_idx(idx)

    def __len__(self):
        return self.n_elems

    def idx_from_coordinate(self, Coordinate xg not None):
        cdef cc.Coordinate xl_xx = self.xx.get_geo().coordinate_l_from_g(xg.xx)
        cdef cc.Long idx = self.xx.f_local_idx.get_elem(xl_xx)
        return idx

    def coordinate_from_idx(self, cc.Long idx):
        cdef Coordinate xg = Coordinate()
        cdef cc.Long index = self.xx.indices[idx]
        cdef cc.Coordinate xl_xx = self.xx.get_geo().coordinate_from_index(index)
        cc.assign_direct(xg.xx, self.xx.get_geo().coordinate_g_from_l(xl_xx))
        return xg

    def __getstate__(self):
        """
        Only work when single node.
        """
        geo = self.geo
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

cache_points_selection = q.mk_cache("points_selection")

@q.timer
def mk_xg_field(Geometry geo):
    cdef FieldInt f = FieldInt()
    cc.set_xg_field(f.xx, geo.xx)
    return f

def get_psel_single(Coordinate total_site, Coordinate xg=None):
    """
    [ xg, ]
    """
    cdef PointsSelection psel
    if xg is None:
        xg = Coordinate([ -1, -1, -1, -1, ])
    param_tuple = (total_site[0], total_site[1], total_site[2], total_site[3], xg[0], xg[1], xg[2], xg[3],)
    if param_tuple not in cache_points_selection:
        psel = PointsSelection(total_site, xg,)
        cache_points_selection[param_tuple] = psel
    return cache_points_selection[param_tuple]

def get_psel_tslice(Coordinate total_site, *, int t_dir=3):
    """
    if t_dir = 3, then [ [-1,-1,-1,0,], [-1,-1,-1,1,], ..., [-1,-1,-1,total_site[3]-1], ]
    if t_dir = 2, then [ [-1,-1,0,-1,], [-1,-1,1,-1,], ..., [-1,-1,total_site[2]-1,-1], ]
    """
    cdef PointsSelection psel
    assert 0 <= t_dir and t_dir < 4
    param_tuple = (total_site[0], total_site[1], total_site[2], total_site[3], t_dir,)
    if param_tuple not in cache_points_selection:
        psel = PointsSelection(total_site)
        cc.assign_direct(psel.xx, cc.mk_tslice_points_selection(total_site.xx, t_dir))
        cache_points_selection[param_tuple] = psel
    return cache_points_selection[param_tuple]

def is_matching_fsel(FieldSelection fsel1, FieldSelection fsel2):
    return cc.is_matching_fsel(fsel1.xx, fsel2.xx)

### -------------------------------------------------------------------
