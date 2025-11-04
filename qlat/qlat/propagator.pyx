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
        assert self.xx.multiplicity == 1 or self.xx.multiplicity == 0
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

    def glb_sum_tslice(self, *, cc.Int t_dir=3):
        cdef SelectedPointsWilsonMatrix sp = super().glb_sum_tslice(t_dir=t_dir)
        cdef PselProp sp_prop = PselProp(sp.psel)
        sp_prop @= sp
        return sp_prop

###

cdef class SelProp(SelectedFieldWilsonMatrix):

    def __init__(self, FieldSelection fsel=None, multiplicity=1):
        assert multiplicity == 1
        super().__init__(fsel, 1)

    cdef cc.Handle[cc.SelProp] xxx(self):
        assert self.xx.multiplicity == 1 or self.xx.multiplicity == 0
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
        assert self.xx.multiplicity == 1 or self.xx.multiplicity == 0
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

cdef class SpinProp(FieldSpinMatrix):

    def __init__(self, Geometry geo=None, multiplicity=1):
        assert multiplicity == 1
        super().__init__(geo, 1)

    cdef cc.Handle[cc.SpinProp] xxx(self):
        assert self.xx.multiplicity == 1 or self.xx.multiplicity == 0
        return cc.Handle[cc.SpinProp](<cc.SpinProp&>self.xx)

    def get_elem_sm(self, cc.Long index, int m=0):
        cdef SpinMatrix sm = SpinMatrix()
        np.asarray(sm)[:] = self[index, m]
        return sm

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

    def glb_sum_tslice(self, *, cc.Int t_dir=3):
        cdef SelectedPointsSpinMatrix sp = super().glb_sum_tslice(t_dir=t_dir)
        return sp

###

def set_point_src(Prop prop_src not None, Geometry geo not None, Coordinate xg not None, cc.PyComplexD value=1.0):
    cc.set_point_src(prop_src.xxx().val(), geo.xx, xg.xx, cc.ccpy_d(value))

def set_wall_src(Prop prop_src not None, Geometry geo not None, int tslice, CoordinateD lmom=None):
    if lmom is None:
        lmom = CoordinateD()
    cc.set_wall_src(prop_src.xxx().val(), geo.xx, tslice, lmom.xx)

def set_rand_vol_u1(
        FieldComplexD fu1 not None,
        Geometry geo not None,
        RngState rs not None,
        ):
    cc.set_rand_vol_u1(fu1.xx, geo.xx, rs.xx)

def set_rand_vol_u1_src(
        Prop prop_src not None,
        FieldComplexD fu1 not None,
        ):
    """
    prop_src ~ fu1
    """
    cc.set_rand_vol_u1_src(prop_src.xxx().val(), fu1.xx)

@q.timer
def mk_point_src(Geometry geo not None, Coordinate xg not None, cc.PyComplexD value=1.0):
    cdef Prop prop_src = Prop(geo)
    set_point_src(prop_src, geo, xg, value)
    return prop_src

@q.timer
def mk_wall_src(Geometry geo not None, int tslice, CoordinateD lmom=None):
    cdef Prop prop_src = Prop(geo)
    set_wall_src(prop_src, geo, tslice, lmom)
    return prop_src

@q.timer
def mk_rand_vol_u1(
        Geometry geo not None,
        RngState rs not None,
        ):
    """
    return prop_src, fu1
    prop_src ~ fu1
    """
    cdef FieldComplexD fu1 = FieldComplexD(geo, 1)
    set_rand_vol_u1(fu1, geo, rs)
    return fu1

@q.timer
def mk_rand_vol_u1_src(
        FieldComplexD fu1 not None,
        ):
    """
    return prop_src
    prop_src ~ fu1
    """
    cdef Prop prop_src = Prop(fu1.geo)
    set_rand_vol_u1_src(prop_src, fu1)
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
def free_invert(prop_src, cc.RealD mass, cc.RealD m5=1.0, CoordinateD momtwist=None):
    cdef Prop qcd_prop_src
    cdef Prop qcd_prop_sol
    cdef SpinProp spin_prop_src
    cdef SpinProp spin_prop_sol
    if momtwist is None:
        momtwist = CoordinateD([ 0.0, 0.0, 0.0, 0.0, ])
    if isinstance(prop_src, Prop):
        qcd_prop_src = prop_src
        qcd_prop_sol = Prop()
        cc.free_invert(qcd_prop_sol.xxx().val(), qcd_prop_src.xxx().val(), mass, m5, momtwist.xx)
        return qcd_prop_sol
    elif isinstance(prop_src, SpinProp):
        spin_prop_src = prop_src
        spin_prop_sol = SpinProp()
        cc.free_invert(spin_prop_sol.xxx().val(), spin_prop_src.xxx().val(), mass, m5, momtwist.xx)
        return spin_prop_sol
    else:
        assert False

@q.timer
def invert_qed(
        SpinProp sp_src, FieldComplexD gf1,
        cc.RealD mass, cc.RealD m5, cc.Int ls,
        cc.Bool is_dagger=False,
        cc.RealD stop_rsd=1e-8, cc.Long max_num_iter=50000,
        ):
    """
    gf1 = q.mk_left_expanded_field(gf)
    """
    cdef SpinProp sp_sol = SpinProp()
    cc.invert_qed(
        sp_sol.xxx().val(), sp_src.xxx().val(), gf1.xx,
        mass, m5, ls, is_dagger, stop_rsd, max_num_iter)
    return sp_sol

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
def mk_ff_list_from_prop(Prop prop):
    """
    return ff_list
    isinstance(ff_list, list)
    len(ff_list) == 12
    """
    cdef cc.Int num_field = 12
    cdef list ff_list = []
    cdef FermionField4d ff
    cdef cc.std_vector[cc.FermionField4d] ff_vec
    cc.set_ff_vec_from_prop(ff_vec, prop.xxx().val())
    assert <cc.Int>ff_vec.size() == num_field
    for i in range(num_field):
        ff = FermionField4d()
        cc.qswap(ff_vec[i], ff.xx)
        ff_list.append(ff)
    return ff_list

@q.timer
def mk_prop_from_ff_list(list ff_list):
    """
    return prop
    isinstance(prop, Prop)
    """
    cdef cc.Int num_field = 12
    assert len(ff_list) == num_field
    cdef Prop prop = Prop()
    cdef FermionField4d ff
    cdef cc.std_vector[cc.FermionField4d] ff_vec
    ff_vec.resize(num_field)
    for i in range(num_field):
        ff = ff_list[i]
        cc.qswap(ff_vec[i], ff.xx)
    cc.set_prop_from_ff_vec(prop.xxx().val(), ff_vec)
    for i in range(num_field):
        ff = ff_list[i]
        cc.qswap(ff_vec[i], ff.xx)
    return prop

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
        assert self.xx.multiplicity == 1 or self.xx.multiplicity == 0
        return cc.Handle[cc.FermionField4d](<cc.FermionField4d&>self.xx)

###
