# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc
from .geometry cimport Geometry
from .field_types cimport (
    FieldRealD,
    FieldColorMatrix,
    FieldComplexD,
)
from .propagator cimport (
    Prop,
    SelProp,
    PselProp,
    FermionField4d,
    SpinProp,
)
from .field_selection cimport (
    FieldSelection,
    PointsSelection,
)

from .field_utils import (
    field_expanded,
    refresh_expanded_1,
)

from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

import cqlat as c
import qlat_utils as q
import numpy as np

cdef class GaugeField(FieldColorMatrix):

    def __init__(self, Geometry geo=None, multiplicity=4):
        assert multiplicity == 4
        super().__init__(geo, 4)

    cdef cc.Handle[cc.GaugeField] xxx(self):
        assert self.xx.multiplicity == 4 or self.xx.multiplicity == 0
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
        assert self.xx.multiplicity == 1 or self.xx.multiplicity == 0
        return cc.Handle[cc.GaugeTransform](<cc.GaugeTransform&>self.xx)

    @q.timer
    def save(self, path):
        """
        Save as RealD precision with the generic Field format ``save_double``
        """
        return self.save_double(path)

    @q.timer
    def load(self, path):
        """
        Load as RealD precision with the generic Field format ``load_double``
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
        cdef GaugeTransform gt, gt1
        cdef GaugeField gf, gf1
        cdef Prop prop, prop1
        cdef SelProp s_prop, s_prop1
        cdef PselProp ps_prop, ps_prop1
        cdef FermionField4d ff, ff1
        cdef FieldSelection fsel
        cdef PointsSelection psel
        if isinstance(other, GaugeTransform):
            gt1 = other
            gt = GaugeTransform()
            cc.gt_apply_gauge_transformation(gt.xxx().val(), gt1.xxx().val(), self.xxx().val())
            return gt
        elif isinstance(other, GaugeField):
            gf1 = other
            gf = GaugeField()
            cc.gf_apply_gauge_transformation(gf.xxx().val(), gf1.xxx().val(), self.xxx().val(), False)
            return gf
        elif isinstance(other, Prop):
            prop1 = other
            prop = Prop()
            cc.prop_apply_gauge_transformation(prop.xxx().val(), prop1.xxx().val(), self.xxx().val())
            return prop
        elif isinstance(other, SelProp):
            s_prop1 = other
            fsel = other.fsel
            s_prop = SelProp(fsel)
            cc.prop_apply_gauge_transformation(s_prop.xx, s_prop1.xx, self.xxx().val(), fsel.xx)
            return s_prop
        elif isinstance(other, PselProp):
            ps_prop1 = other
            psel = other.psel
            ps_prop = PselProp(psel)
            cc.prop_apply_gauge_transformation(ps_prop.xx, ps_prop1.xx, self.xxx().val(), psel.xx)
            return ps_prop
        elif isinstance(other, FermionField4d):
            ff1 = other
            ff = FermionField4d()
            cc.ff_apply_gauge_transformation(ff.xxx().val(), ff1.xxx().val(), self.xxx().val())
            return ff
        elif isinstance(other, list):
            return [ self * p for p in other ]
        else:
            raise Exception("GaugeTransform.__mul__")

    def inv(self):
        cdef GaugeTransform gt_inv = GaugeTransform()
        cc.gt_invert(gt_inv.xxx().val(), self.xxx().val())
        return gt_inv

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

def gf_plaq_field(GaugeField gf):
    assert gf is not None
    cdef FieldRealD f_plaq = FieldRealD()
    cc.gf_plaq_field(f_plaq.xx, gf.xxx().val())
    return f_plaq

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

def gf_twist_boundary_at_boundary(GaugeField gf, cc.RealD lmom=-0.5, int mu=3):
    """
    modify gf in place
    """
    c.gf_twist_boundary_at_boundary(gf, lmom, mu)

def mk_left_expanded_field(gf):
    """
    Return left expanded field.
    Similar to ``set_left_expanded_gauge_field`` in C++
    """
    gf1 = field_expanded(gf, 1, 0)
    refresh_expanded_1(gf1)
    return gf1

def gf_reduce_half(GaugeField gf):
    """
    return hgf with half of the size in all directions
    Do not modify `gf`.
    Can perform some shift before call this function, e.g.
    `gf_reduce_half(gf.shift())`
    """
    cdef GaugeField hgf = GaugeField()
    cc.gf_reduce_half(hgf.xxx().val(), gf.xxx().val())
    return hgf

###

@q.timer
def invert_dwf_qed(
        FieldComplexD f_in4d, FieldComplexD gf1,
        cc.RealD mass, cc.RealD m5, cc.Int ls,
        *,
        t_wick_phase_factor_arr=None,
        cc.Bool is_dagger=False,
        cc.RealD stop_rsd=1e-8, cc.Long max_num_iter=50000,
        ):
    """
    properly project to 4d fermion field
    if is_dagger is false (default), then M out = in
    if is_dagger is true, then M^dag out = in
    #
    gf1 = q.mk_left_expanded_field(gf)
    """
    cdef FieldComplexD f_out4d = FieldComplexD()
    cdef cc.vector[cc.ComplexD] t_wick_phase_factor_vec = cc.vector[cc.ComplexD]()
    cdef cc.Int t_size
    cdef cc.Int i
    if t_wick_phase_factor_arr is not None:
        t_size = len(t_wick_phase_factor_arr)
        t_wick_phase_factor_vec.resize(t_size)
        for i in range(t_size):
            t_wick_phase_factor_vec[i] = cc.ccpy_d(t_wick_phase_factor_arr[i])
    cc.invert_dwf_qed(
        f_out4d.xx, f_in4d.xx, gf1.xx,
        mass, m5, ls, t_wick_phase_factor_vec, is_dagger, stop_rsd, max_num_iter)
    return f_out4d

@q.timer
def cg_with_m_dwf_qed(
        FieldComplexD f_in5d, FieldComplexD gf1,
        cc.RealD mass, cc.RealD m5, cc.Int ls,
        *,
        t_wick_phase_factor_arr=None,
        cc.Bool is_dagger=False,
        cc.RealD stop_rsd=1e-8, cc.Long max_num_iter=50000,
        ):
    """
    if is_dagger is false, then M^dag M out = in
    if is_dagger is true, then M M^dag out = in
    #
    gf1 = q.mk_left_expanded_field(gf)
    """
    cdef FieldComplexD f_out5d = FieldComplexD()
    cdef cc.vector[cc.ComplexD] t_wick_phase_factor_vec = cc.vector[cc.ComplexD]()
    cdef cc.Int t_size
    cdef cc.Int i
    if t_wick_phase_factor_arr is not None:
        t_size = len(t_wick_phase_factor_arr)
        t_wick_phase_factor_vec.resize(t_size)
        for i in range(t_size):
            t_wick_phase_factor_vec[i] = cc.ccpy_d(t_wick_phase_factor_arr[i])
    cc.cg_with_m_dwf_qed(
        f_out5d.xx, f_in5d.xx, gf1.xx,
        mass, m5, ls, t_wick_phase_factor_vec, is_dagger, stop_rsd, max_num_iter)
    return f_out5d

@q.timer
def multiply_m_dwf_qed(
        FieldComplexD f_in5d, FieldComplexD gf1,
        cc.RealD mass, cc.RealD m5, cc.Int ls,
        *,
        t_wick_phase_factor_arr=None,
        cc.Bool is_dagger=False,
        ):
    """
    gf1 = q.mk_left_expanded_field(gf)
    """
    cdef FieldComplexD f_out5d = FieldComplexD()
    cdef cc.vector[cc.ComplexD] t_wick_phase_factor_vec = cc.vector[cc.ComplexD]()
    cdef cc.Int t_size
    cdef cc.Int i
    if t_wick_phase_factor_arr is not None:
        t_size = len(t_wick_phase_factor_arr)
        t_wick_phase_factor_vec.resize(t_size)
        for i in range(t_size):
            t_wick_phase_factor_vec[i] = cc.ccpy_d(t_wick_phase_factor_arr[i])
    cc.multiply_m_dwf_qed(f_out5d.xx, f_in5d.xx, gf1.xx, mass, m5, ls, t_wick_phase_factor_vec, is_dagger)
    return f_out5d

###