# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc
from .geometry cimport Geometry
from .field_types cimport (
        FieldComplexD,
        )
from .field_selection cimport (
        FieldSelection,
        PointsSelection,
        )

from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT

import cqlat as c
import qlat_utils as q
import numpy as np

from .field_utils import mk_fft

cdef class Prop(FieldWilsonMatrix):

    def __init__(self, Geometry geo=None, multiplicity=1):
        assert multiplicity == 1
        super().__init__(geo, 1)

    cdef cc.Handle[cc.Prop] xxx(self):
        assert self.xx.multiplicity == 1
        return cc.Handle[cc.Prop](<cc.Prop&>self.xx)

    def get_elem_wm(self, cc.Long index, int m=0):
        cdef WilsonMatrix wm = WilsonMatrix()
        np.asarray(wm)[:] = self[index, m]
        return wm

    def __getstate__(self):
        """
        Only work when single node (or if all nodes has the same data).
        """
        return super().__getstate__()

    def __setstate__(self, state):
        """
        Only work when single node (or if all nodes has the same data).
        """
        super().__setstate__(state)

###

cdef class SelProp(SelectedFieldWilsonMatrix):

    def __init__(self, FieldSelection fsel=None, multiplicity=1):
        assert multiplicity == 1
        super().__init__(fsel, 1)

    cdef cc.Handle[cc.SelProp] xxx(self):
        assert self.xx.multiplicity == 1
        return cc.Handle[cc.SelProp](<cc.SelProp&>self.xx)

    def get_elem_wm(self, cc.Long idx, int m=0):
        cdef WilsonMatrix wm = WilsonMatrix()
        np.asarray(wm)[:] = self[idx, m]
        return wm

    def __getstate__(self):
        """
        Only work when single node (or if all nodes has the same data).
        """
        return super().__getstate__()

    def __setstate__(self, state):
        """
        Only work when single node (or if all nodes has the same data).
        """
        super().__setstate__(state)

###

cdef class PselProp(SelectedPointsWilsonMatrix):

    def __init__(self, *args):
        cdef cc.Int len_args = len(args)
        if len_args == 0:
            super().__init__()
        elif isinstance(args[0], PointsSelection):
            if len(args) == 1:
                psel, = args
            else:
                psel, multiplicity, = args
                assert multiplicity == 1
            super().__init__(psel, 1)
        elif isinstance(args[0], SelProp):
            super().__init__(*args)

    cdef cc.Handle[cc.PselProp] xxx(self):
        assert self.xx.multiplicity == 1
        return cc.Handle[cc.PselProp](<cc.PselProp&>self.xx)

    def get_elem_wm(self, cc.Long idx, int m=0):
        cdef WilsonMatrix wm = WilsonMatrix()
        np.asarray(wm)[:] = self[idx, m]
        return wm

    def __getstate__(self):
        """
        Only work when single node (or if all nodes has the same data).
        """
        return super().__getstate__()

    def __setstate__(self, state):
        """
        Only work when single node (or if all nodes has the same data).
        """
        super().__setstate__(state)

###

def set_point_src(Prop prop_src not None, Geometry geo not None, Coordinate xg not None, cc.PyComplexD value=1.0):
    cc.set_point_src(prop_src.xxx().p[0], geo.xx, xg.xx, cc.ccpy_d(value))

def set_wall_src(Prop prop_src not None, Geometry geo not None, int tslice, CoordinateD lmom=None):
    if lmom is None:
        lmom = CoordinateD()
    cc.set_wall_src(prop_src.xxx().p[0], geo.xx, tslice, lmom.xx)

def mk_point_src(Geometry geo not None, Coordinate xg not None, cc.PyComplexD value=1.0):
    cdef Prop prop_src = Prop(geo)
    set_point_src(prop_src, geo, xg, value)
    return prop_src

def mk_wall_src(Geometry geo not None, int tslice, CoordinateD lmom=None):
    cdef Prop prop_src = Prop(geo)
    set_wall_src(prop_src, geo, tslice, lmom)
    return prop_src

@q.timer
def mk_rand_u1_src(sel, rs):
    """
    return (prop_src, fu1,) where prop_src = Prop() and fu1 = FieldComplex
    fu1 stores the random u1 numbers (fu1.multiplicity == 1)
    sel can be psel or fsel
    """
    prop_src = Prop()
    fu1 = FieldComplexD()
    if isinstance(sel, FieldSelection):
        fsel = sel
        c.set_rand_u1_src_fsel(prop_src, fu1, fsel, rs)
    elif isinstance(sel, PointsSelection):
        psel = sel
        geo = psel.geo
        assert isinstance(geo, Geometry)
        c.set_rand_u1_src_psel(prop_src, fu1, psel, geo, rs)
    else:
        raise Exception(f"mk_rand_u1_src {type(sel)}")
    return (prop_src, fu1,)

@q.timer
def get_rand_u1_sol(Prop prop_sol, FieldComplexD fu1, sel):
    assert isinstance(prop_sol, Prop)
    assert isinstance(fu1, FieldComplexD)
    if isinstance(sel, FieldSelection):
        fsel = sel
        s_prop = SelProp(fsel)
        c.set_rand_u1_sol_fsel(s_prop, prop_sol, fu1, fsel)
        return s_prop
    elif isinstance(sel, PointsSelection):
        psel = sel
        sp_prop = PselProp(psel)
        c.set_rand_u1_sol_psel(sp_prop, prop_sol, fu1, psel)
        return sp_prop
    else:
        raise Exception(f"get_rand_u1_sol {type(sel)}")

@q.timer_verbose
def mk_rand_u1_prop(inv, sel, rs):
    """
    interface function
    return s_prop
    sel can be psel or fsel
    """
    prop_src, fu1 = mk_rand_u1_src(sel, rs)
    prop_sol = inv * prop_src
    return get_rand_u1_sol(prop_sol, fu1, sel)

@q.timer
def free_invert(prop_src, mass, m5=1.0, momtwist=None):
    assert isinstance(prop_src, Prop)
    if momtwist is None:
        momtwist = [ 0.0, 0.0, 0.0, 0.0, ]
    prop_sol = Prop()
    c.free_invert_prop(prop_sol, prop_src, mass, m5, momtwist)
    return prop_sol

def convert_mspincolor_from_wm(prop_wm):
    prop_msc = prop_wm.copy(False)
    if isinstance(prop_wm, Prop):
        c.convert_mspincolor_from_wm_prop(prop_msc, prop_wm)
    elif isinstance(prop_wm, SelProp):
        c.convert_mspincolor_from_wm_s_prop(prop_msc, prop_wm)
    elif isinstance(prop_wm, PselProp):
        c.convert_mspincolor_from_wm_sp_prop(prop_msc, prop_wm)
    else:
        raise Exception("prop type match failed")
    return prop_msc

def convert_wm_from_mspincolor(prop_msc):
    prop_wm = prop_msc.copy(False)
    if isinstance(prop_msc, Prop):
        c.convert_wm_from_mspincolor_prop(prop_wm, prop_msc)
    elif isinstance(prop_msc, SelProp):
        c.convert_wm_from_mspincolor_s_prop(prop_wm, prop_msc)
    elif isinstance(prop_msc, PselProp):
        c.convert_wm_from_mspincolor_sp_prop(prop_wm, prop_msc)
    else:
        raise Exception("prop type match failed")
    return prop_wm

@q.timer
def flip_tpbc_with_tslice(prop, tslice_flip_tpbc):
    if isinstance(prop, SelProp):
        c.flip_tpbc_with_tslice_s_prop(prop, tslice_flip_tpbc)
    elif isinstance(prop, PselProp):
        c.flip_tpbc_with_tslice_sp_prop(prop, tslice_flip_tpbc)
    else:
        print(type(prop))
        assert False

@q.timer
def free_scalar_invert_mom_cfield(FieldComplexD f, mass):
    assert isinstance(f, FieldComplexD)
    c.free_scalar_invert_mom_cfield(f, mass)

@q.timer
def free_scalar_invert_cfield(src, mass, *, mode_fft=1):
    fft_f = mk_fft(is_forward=True, is_normalizing=True, mode_fft=mode_fft)
    fft_b = mk_fft(is_forward=False, is_normalizing=True, mode_fft=mode_fft)
    f = fft_f * src
    free_scalar_invert_mom_cfield(f, mass)
    sol = fft_b * f
    return sol

cdef class FermionField4d(FieldWilsonVector):

    def __init__(self, Geometry geo=None, multiplicity=1):
        assert multiplicity == 1
        super().__init__(geo, 1)

    cdef cc.Handle[cc.FermionField4d] xxx(self):
        assert self.xx.multiplicity == 1
        return cc.Handle[cc.FermionField4d](<cc.FermionField4d&>self.xx)

###
