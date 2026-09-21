# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.propagator``
===========================\n
Lattice propagator types and utilities for domain-wall fermion QCD
simulations.  Defines the core propagator containers (`Prop`, `SelProp`,
`PselProp`, `SpinProp`, `FermionField4d`) together with functions for
constructing point, wall, and random-U1 sources, performing free and
conjugate-gradient inversions, converting between Wilson-matrix and
spin-color layouts, and splitting/reassembling propagators into
individual fermion fields.\n
Documentation: ``docs/qlat/qlat_propagator.md``\n
.. note:: Update the documentation when updating this source file.
"""

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

class q:
    from qlat_utils import (
        timer,
    )

import numpy as np

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
        else:
            assert False

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
        cc.Int multiplicity,
        RngState rs not None,
        ):
    cc.set_rand_vol_u1(fu1.xx, geo.xx, multiplicity, rs.xx)

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
        cc.Int multiplicity,
        RngState rs not None,
        ):
    """
    return prop_src, fu1
    prop_src ~ fu1
    """
    cdef FieldComplexD fu1 = FieldComplexD(geo, multiplicity)
    set_rand_vol_u1(fu1, geo, multiplicity, rs)
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
    cdef Prop prop_src = Prop()
    cdef FieldComplexD fu1 = FieldComplexD()
    cdef FieldSelection fsel
    cdef PointsSelection psel
    cdef Geometry geo
    cdef RngState rs_c = <RngState>rs
    if isinstance(sel, FieldSelection):
        fsel = <FieldSelection>sel
        cc.py_set_rand_u1_src_fsel(
            prop_src.xxx().val(), fu1.xx, fsel.xx, rs_c.xx)
    elif isinstance(sel, PointsSelection):
        psel = <PointsSelection>sel
        geo = psel.geo
        assert isinstance(geo, Geometry)
        cc.py_set_rand_u1_src_psel(
            prop_src.xxx().val(), fu1.xx, psel.xx, geo.xx, rs_c.xx)
    else:
        raise Exception(f"mk_rand_u1_src {type(sel)}")
    return (prop_src, fu1,)

@q.timer
def get_rand_u1_sol(Prop prop_sol, FieldComplexD fu1, sel):
    assert isinstance(prop_sol, Prop)
    assert isinstance(fu1, FieldComplexD)
    cdef SelProp s_prop
    cdef PselProp sp_prop
    cdef FieldSelection fsel
    cdef PointsSelection psel
    if isinstance(sel, FieldSelection):
        fsel = <FieldSelection>sel
        s_prop = SelProp(fsel)
        cc.py_set_rand_u1_sol_fsel(
            s_prop.xxx().val(), prop_sol.xxx().val(), fu1.xx, fsel.xx)
        return s_prop
    elif isinstance(sel, PointsSelection):
        psel = <PointsSelection>sel
        sp_prop = PselProp(psel)
        cc.py_set_rand_u1_sol_psel(
            sp_prop.xxx().val(), prop_sol.xxx().val(), fu1.xx, psel.xx)
        return sp_prop
    else:
        raise Exception(f"get_rand_u1_sol {type(sel)}")

@q.timer
def free_invert(prop_src, cc.RealD mass, cc.RealD m5=1.0, CoordinateD momtwist=None):
    """
    Compute the free (gauge-field-independent) inverse of a propagator source.\n
    The source is transformed to momentum space, multiplied by the analytic
    free DWF propagator, and transformed back.  Supports both `Prop` and
    `SpinProp` input; returns the same type.
    """
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
def free_mom_invert(prop_src, cc.RealD mass, cc.RealD m5=1.0, CoordinateD momtwist=None):
    """
    Apply the free DWF inverse in momentum space.\n
    Unlike `free_invert`, no Fourier transform is performed: `prop_src` is
    assumed to already be in momentum space (e.g. the output of a normalizing
    forward FFT).  Supports both `Prop` and `SpinProp` input; returns the same
    type.  All parameters of the C++ `free_mom_invert` kernel are exposed.
    """
    cdef Prop qcd_prop_src
    cdef Prop qcd_prop_sol
    cdef SpinProp spin_prop_src
    cdef SpinProp spin_prop_sol
    if momtwist is None:
        momtwist = CoordinateD([ 0.0, 0.0, 0.0, 0.0, ])
    if isinstance(prop_src, Prop):
        qcd_prop_src = prop_src
        qcd_prop_sol = Prop()
        cc.free_mom_invert(qcd_prop_sol.xxx().val(), qcd_prop_src.xxx().val(), mass, m5, momtwist.xx)
        return qcd_prop_sol
    elif isinstance(prop_src, SpinProp):
        spin_prop_src = prop_src
        spin_prop_sol = SpinProp()
        cc.free_mom_invert(spin_prop_sol.xxx().val(), spin_prop_src.xxx().val(), mass, m5, momtwist.xx)
        return spin_prop_sol
    else:
        assert False

@q.timer
def invert_qed(
        SpinProp sp_src, FieldComplexD gf1,
        cc.RealD mass, cc.RealD m5, cc.Int ls,
        *,
        t_wick_phase_factor_arr=None,
        cc.Bool is_dagger=False,
        cc.RealD stop_rsd=1e-8, cc.Long max_num_iter=50000,
        ):
    """
    gf1 = q.mk_left_expanded_field(gf)
    """
    cdef SpinProp sp_sol = SpinProp()
    cdef cc.vector[cc.ComplexD] t_wick_phase_factor_vec = cc.vector[cc.ComplexD]()
    cdef cc.Int t_size
    cdef cc.Int i
    if t_wick_phase_factor_arr is not None:
        t_size = len(t_wick_phase_factor_arr)
        t_wick_phase_factor_vec.resize(t_size)
        for i in range(t_size):
            t_wick_phase_factor_vec[i] = cc.ccpy_d(t_wick_phase_factor_arr[i])
    cc.invert_qed(
        sp_sol.xxx().val(), sp_src.xxx().val(), gf1.xx,
        mass, m5, ls, t_wick_phase_factor_vec,
        is_dagger, stop_rsd, max_num_iter)
    return sp_sol

def convert_mspincolor_from_wm(prop_wm):
    cdef Prop prop_msc_prop
    cdef SelProp prop_msc_s_prop
    cdef PselProp prop_msc_sp_prop
    prop_msc = prop_wm.copy(False)
    if isinstance(prop_wm, Prop):
        prop_msc_prop = prop_msc
        cc.convert_mspincolor_from_wm(
            prop_msc_prop.xxx().val(), (<Prop>prop_wm).xxx().val())
    elif isinstance(prop_wm, SelProp):
        prop_msc_s_prop = prop_msc
        cc.convert_mspincolor_from_wm(
            prop_msc_s_prop.xxx().val(), (<SelProp>prop_wm).xxx().val())
    elif isinstance(prop_wm, PselProp):
        prop_msc_sp_prop = prop_msc
        cc.convert_mspincolor_from_wm(
            prop_msc_sp_prop.xxx().val(), (<PselProp>prop_wm).xxx().val())
    else:
        raise Exception("prop type match failed")
    return prop_msc

def convert_wm_from_mspincolor(prop_msc):
    cdef Prop prop_wm_prop
    cdef SelProp prop_wm_s_prop
    cdef PselProp prop_wm_sp_prop
    prop_wm = prop_msc.copy(False)
    if isinstance(prop_msc, Prop):
        prop_wm_prop = prop_wm
        cc.convert_wm_from_mspincolor(
            prop_wm_prop.xxx().val(), (<Prop>prop_msc).xxx().val())
    elif isinstance(prop_msc, SelProp):
        prop_wm_s_prop = prop_wm
        cc.convert_wm_from_mspincolor(
            prop_wm_s_prop.xxx().val(), (<SelProp>prop_msc).xxx().val())
    elif isinstance(prop_msc, PselProp):
        prop_wm_sp_prop = prop_wm
        cc.convert_wm_from_mspincolor(
            prop_wm_sp_prop.xxx().val(), (<PselProp>prop_msc).xxx().val())
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
    cdef PselProp ps
    cdef Geometry geo
    if isinstance(prop, SelProp):
        cc.flip_tpbc_with_tslice((<SelProp>prop).xxx().val(),
                                 (<SelProp>prop).fsel.xx,
                                 tslice_flip_tpbc)
    elif isinstance(prop, PselProp):
        ps = <PselProp>prop
        geo = ps.psel.geo
        cc.flip_tpbc_with_tslice(ps.xxx().val(), ps.psel.xx,
                                 tslice_flip_tpbc, geo.total_site[3])
    else:
        print(type(prop))
        assert False

### -------------------------------------------------------------------
### cqlat-compatible entry points

def set_rand_u1_src_psel(Prop prop, FieldComplexD fu1, PointsSelection psel,
                         Geometry geo, RngState rs):
    cc.py_set_rand_u1_src_psel(prop.xxx().val(), fu1.xx, psel.xx, geo.xx,
                               rs.xx)

def set_rand_u1_sol_psel(PselProp sp_prop, Prop prop, FieldComplexD fu1,
                         PointsSelection psel):
    cc.py_set_rand_u1_sol_psel(sp_prop.xxx().val(), prop.xxx().val(), fu1.xx,
                               psel.xx)

def set_rand_u1_src_fsel(Prop prop, FieldComplexD fu1, FieldSelection fsel,
                         RngState rs):
    cc.py_set_rand_u1_src_fsel(prop.xxx().val(), fu1.xx, fsel.xx, rs.xx)

def set_rand_u1_sol_fsel(SelProp sf_prop, Prop prop, FieldComplexD fu1,
                         FieldSelection fsel):
    cc.py_set_rand_u1_sol_fsel(sf_prop.xxx().val(), prop.xxx().val(), fu1.xx,
                               fsel.xx)

def flip_tpbc_with_tslice_sp_prop(PselProp sp_prop,
                                  cc.Int tslice_flip_tpbc):
    cdef Geometry geo = sp_prop.psel.geo
    cc.flip_tpbc_with_tslice(sp_prop.xxx().val(), sp_prop.psel.xx,
                             tslice_flip_tpbc, geo.total_site[3])

def flip_tpbc_with_tslice_s_prop(SelProp s_prop, cc.Int tslice_flip_tpbc):
    cc.flip_tpbc_with_tslice(s_prop.xxx().val(), s_prop.fsel.xx,
                             tslice_flip_tpbc)

@q.timer
def free_scalar_mom_invert(FieldComplexD f, mass, CoordinateD momtwist=None):
    """
    Apply the free scalar inverse in momentum space, in-place.\n
    `f` is assumed to already be in momentum space (e.g. the output of a
    normalizing forward FFT).  `momtwist` is exposed and passed through to the
    C++ kernel; when None it defaults to zero.
    """
    if momtwist is None:
        momtwist = CoordinateD([ 0.0, 0.0, 0.0, 0.0, ])
    cc.free_scalar_mom_invert(f.xx, mass, momtwist.xx)

@q.timer
def free_scalar_deriv_mom(FieldComplexD f, deriv, CoordinateD momtwist=None):
    """
    Apply a lattice derivative in momentum space, in-place.\n
    `f` is assumed to already be in momentum space (e.g. the output of a
    normalizing forward FFT).  `deriv` gives the derivative order for each
    direction `x, y, z, t` and the field is multiplied by\n
        prod_mu ( 2 i sin(k_mu / 2) )^{deriv[mu]}\n
    with `k_mu = 2 pi ( smod(n_mu, L_mu) + momtwist_mu ) / L_mu`.\n
    This is the bare derivative factor: it contains no mass and no `1 / D(k)`.
    Compose it with `free_scalar_mom_invert` to differentiate the free scalar
    inverse; the two factors commute.  Note `d_mu^2 = -4 sin^2(k_mu / 2)`, i.e.
    minus the mu term of the `D(k)` used by `free_scalar_mom_invert`.\n
    At the self-conjugate momentum `k_mu = pi` the two branches of
    `2 i sin(k_mu / 2)` differ by a sign.  For an odd `deriv[mu]` that sign is
    ambiguous, so that mode is dropped (`d_mu = 0`); for an even `deriv[mu]`
    the sign squares out and the mode is kept.\n
    `deriv=None` is equivalent to `[0, 0, 0, 0]` and leaves `f` unchanged.
    """
    cdef cc.vector[cc.Int] deriv_vec = cc.vector[cc.Int]()
    cdef cc.Int i
    cdef cc.Int n
    if momtwist is None:
        momtwist = CoordinateD([ 0.0, 0.0, 0.0, 0.0, ])
    if deriv is None:
        deriv = [ 0, 0, 0, 0 ]
    else:
        deriv = list(deriv)
    if len(deriv) != 4:
        raise Exception(f"free_scalar_deriv_mom: deriv={deriv} must have length 4")
    deriv_vec.resize(4)
    for i in range(4):
        n = int(deriv[i])
        if n < 0:
            raise Exception(f"free_scalar_deriv_mom: deriv={deriv} must be non-negative")
        deriv_vec[i] = n
    cc.free_scalar_deriv_mom(f.xx, deriv_vec, momtwist.xx)

cdef class FermionField4d(FieldWilsonVector):

    def __init__(self, Geometry geo=None, multiplicity=1):
        assert multiplicity == 1
        super().__init__(geo, 1)

    cdef cc.Handle[cc.FermionField4d] xxx(self):
        assert self.xx.multiplicity == 1 or self.xx.multiplicity == 0
        return cc.Handle[cc.FermionField4d](<cc.FermionField4d&>self.xx)

###
