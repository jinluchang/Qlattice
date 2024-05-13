# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc
from .geometry cimport Geometry
from .field_types cimport (
        FieldColorMatrix,
        )

from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

import cqlat as c
import qlat_utils as q
import numpy as np

from .field_utils import (
        field_expanded,
        refresh_expanded_1,
        )

cdef class GaugeField(FieldColorMatrix):

    def __init__(self, Geometry geo=None, multiplicity=4):
        assert multiplicity == 4
        super().__init__(geo, 4)

    cdef cc.Handle[cc.GaugeField] xxx(self):
        assert self.xx.multiplicity == 4
        return cc.Handle[cc.GaugeField](<cc.GaugeField&>self.xx)

    @q.timer
    def save(self, path):
        """
        Save with the standard NERSC format
        """
        return c.save_gauge_field(self, path)

    @q.timer
    def load(self, path):
        """
        Load with the standard NERSC format
        """
        return c.load_gauge_field(self, path)

    def set_rand(self, RngState rng, cc.RealD sigma=0.5, cc.Int n_step=1):
        set_g_rand_color_matrix_field(self, rng, sigma, n_step)

    def unitarize(self):
        c.unitarize_color_matrix_field(self)

    def plaq(self):
        return gf_avg_plaq(self)

    def link_trace(self):
        return gf_avg_link_trace(self)

    def twist_boundary_at_boundary(self, cc.RealD lmom=-0.5, cc.Int mu=3):
        """
        modify in place
        """
        gf_twist_boundary_at_boundary(self, lmom, mu)

    def show_info(self):
        gf_show_info(self)

###

cdef class GaugeTransform(FieldColorMatrix):

    def __init__(self, Geometry geo=None, multiplicity=1):
        assert multiplicity == 1
        super().__init__(geo, 1)

    cdef cc.Handle[cc.GaugeTransform] xxx(self):
        assert self.xx.multiplicity == 1
        return cc.Handle[cc.GaugeTransform](<cc.GaugeTransform&>self.xx)

    @q.timer
    def save(self, path):
        """
        Save as double precision with the generic Field format ``save_double``
        """
        return self.save_double(path)

    @q.timer
    def load(self, path):
        """
        Load as double precision with the generic Field format ``load_double``
        """
        return self.load_double(path)

    @q.timer
    def save_cps(self, path):
        """
        Save with the format used in CPS
        """
        return c.save_gauge_transform_cps(self, path)

    @q.timer
    def load_cps(self, path):
        """
        Load with the format used in CPS
        """
        return c.load_gauge_transform_cps(self, path)

    def set_rand(self, RngState rng, cc.RealD sigma=0.5, cc.Int n_step=1):
        set_g_rand_color_matrix_field(self, rng, sigma, n_step)

    def unitarize(self):
        c.unitarize_color_matrix_field(self)

    def __mul__(self, other):
        """
        other can be GaugeTransform, GaugeField, Prop, SelProp, PselProp, list
        """
        from qlat.propagator import Prop, SelProp, PselProp
        if isinstance(other, GaugeTransform):
            gt = GaugeTransform()
            c.apply_gt_gt(gt, self, other)
            return gt
        elif isinstance(other, GaugeField):
            gf = GaugeField()
            c.apply_gt_gf(gf, self, other)
            return gf
        elif isinstance(other, Prop):
            prop = Prop()
            c.apply_gt_prop(prop, self, other)
            return prop
        elif isinstance(other, SelProp):
            prop = SelProp(other.fsel)
            c.apply_gt_sprop(prop, self, other)
            return prop
        elif isinstance(other, PselProp):
            prop = PselProp(other.psel)
            c.apply_gt_psprop(prop, self, other)
            return prop
        elif isinstance(other, list):
            return [ self * p for p in other ]
        else:
            raise Exception("GaugeTransform.__mul__")

    def inv(self):
        gt = GaugeTransform()
        c.gt_invert(gt, self)
        return gt

###

@q.timer_verbose
def gf_show_info(GaugeField gf):
    assert gf is not None
    q.displayln_info(f"gf_show_info: plaq = {gf.plaq():.16F} ; link_trace = {gf.link_trace():.16F}.")

def gf_avg_plaq(GaugeField gf):
    assert gf is not None
    return cc.gf_avg_plaq(gf.xxx().val())

def gf_avg_spatial_plaq(GaugeField gf):
    assert gf is not None
    return cc.gf_avg_spatial_plaq(gf.xxx().val())

def gf_avg_link_trace(GaugeField gf):
    assert gf is not None
    return cc.gf_avg_link_trace(gf.xxx().val())

def gf_wilson_line_no_comm(wlf, m, gf_ext, path, path_n=None):
    """
    wlf = FieldColorMatrix(geo)
    will only modify the m'th component
    e.g. path = [ mu, mu, nu, -mu-1, -mu-1, ]
    e.g. path = [ mu, nu, -mu-1, ], path_n = [ 2, 1, 2, ]
    """
    if path_n is None:
        c.gf_wilson_line_no_comm(wlf, m, gf_ext, path)
    else:
        c.gf_wilson_line_no_comm(wlf, m, gf_ext, path, path_n)

def gf_wilson_lines_no_comm(gf_ext, path_list):
    """
    path_list = [ path_spec, ... ]
    e.g. path_spec = [ mu, mu, nu, -mu-1, -mu-1, ]
    e.g. path_spec = ([ mu, nu, -mu-1, ], [ 2, 1, 2, ],)
    return wlf
    """
    multiplicity = len(path_list)
    geo = q.geo_resize(gf_ext.geo)
    wlf = FieldColorMatrix(geo, multiplicity)
    for m, p in enumerate(path_list):
        if isinstance(p, tuple) and len(p) == 2:
            path, path_n = p
            gf_wilson_line_no_comm(wlf, m, gf_ext, path, path_n)
        else:
            path = p
            gf_wilson_line_no_comm(wlf, m, gf_ext, path)
    return wlf

def gf_avg_wilson_loop_normalized_tr(gf, l, t):
    assert isinstance(gf, GaugeField)
    assert isinstance(l, int)
    assert isinstance(t, int)
    return c.gf_avg_wilson_loop_normalized_tr(gf, l, t)

def set_g_rand_color_matrix_field(fc, rng, sigma, n_steps=1):
    assert isinstance(fc, FieldColorMatrix)
    assert isinstance(rng, RngState)
    return c.set_g_rand_color_matrix_field(fc, rng, sigma, n_steps)

def gf_twist_boundary_at_boundary(GaugeField gf, double lmom=-0.5, int mu=3):
    """
    modify gf in place
    """
    c.gf_twist_boundary_at_boundary(gf, lmom, mu)

def mk_left_expanded_gauge_field(GaugeField gf):
    cdef GaugeField gf1 = field_expanded(gf, 1, 0)
    refresh_expanded_1(gf1)
    return gf1

