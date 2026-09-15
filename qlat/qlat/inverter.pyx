# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.inverter``
=========================\n
Fermion-matrix inverter framework for domain-wall fermions.\n
Provides a hierarchy of ``Inverter`` classes that apply the inverse of the
Dirac operator to propagator sources. Concrete implementations include
analytic free-field inversion, CG-based domain-wall inversion via C/C++
backends, and a gauge-transform wrapper that inverts in a different gauge.\n
Documentation: ``docs/qlat/qlat_inverter.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils.all cimport *
from qlat_utils import mk_cache
from . cimport everything as cc
from .qcd cimport GaugeField
from .qcd cimport GaugeTransform
from .propagator cimport Prop
from .propagator cimport FermionField4d
from .propagator import free_invert

from cpython.long cimport PyLong_FromVoidPtr
from cpython.long cimport PyLong_AsVoidPtr

cache_inv = mk_cache("inv")

cdef inline cc.InverterDomainWall* get_inverter_domain_wall_ptr(
        object inv) except? NULL:
    return <cc.InverterDomainWall*>PyLong_AsVoidPtr(inv.cdata)

def free_inverter_domain_wall(inv):
    cdef cc.InverterDomainWall* pinv = get_inverter_domain_wall_ptr(inv)
    del pinv

### -------------------------------------------------------------------
### cqlat-compatible entry points

def get_stop_rsd_inverter_domain_wall(inv):
    return get_inverter_domain_wall_ptr(inv).stop_rsd()

def set_stop_rsd_inverter_domain_wall(inv, stop_rsd):
    cc.py_set_stop_rsd_inverter_domain_wall(
        get_inverter_domain_wall_ptr(inv)[0], stop_rsd)

def get_max_num_iter_inverter_domain_wall(inv):
    return get_inverter_domain_wall_ptr(inv).max_num_iter()

def set_max_num_iter_inverter_domain_wall(inv, max_num_iter):
    cc.py_set_max_num_iter_inverter_domain_wall(
        get_inverter_domain_wall_ptr(inv)[0], max_num_iter)

def get_max_mixed_precision_cycle_inverter_domain_wall(inv):
    return get_inverter_domain_wall_ptr(inv).max_mixed_precision_cycle()

def set_max_mixed_precision_cycle_inverter_domain_wall(
        inv, max_mixed_precision_cycle):
    cc.py_set_max_mixed_precision_cycle_inverter_domain_wall(
        get_inverter_domain_wall_ptr(inv)[0], max_mixed_precision_cycle)

cdef inline cc.FermionAction* get_fermion_action_ptr(object fa) except? NULL:
    return <cc.FermionAction*>PyLong_AsVoidPtr(fa.cdata)

def mk_inverter_domain_wall(GaugeField gf, fa):
    cdef cc.InverterDomainWall* pinv = new cc.InverterDomainWall()
    cdef cc.FermionAction* pfa = get_fermion_action_ptr(fa)
    cc.setup_inverter(pinv[0], gf.xxx().val(), pfa[0])
    return PyLong_FromVoidPtr(<void*>pinv)

def invert_inverter_domain_wall(Prop prop_sol, Prop prop_src, inv):
    cdef cc.InverterDomainWall* pinv = get_inverter_domain_wall_ptr(inv)
    cc.invert_prop_dw(prop_sol.xxx().val(), prop_src.xxx().val(), pinv[0])

class Inverter:
    pass

## -----

class InverterDwfFreeField(Inverter):
    """
    self.mass
    self.m5
    self.momtwist
    self.timer
    """

    def __init__(self, *, mass, m5=1.0, momtwist=None, qtimer=TimerNone()):
        if momtwist is None:
            momtwist = CoordinateD([ 0.0, 0.0, 0.0, 0.0, ])
        self.mass = mass
        self.m5 = m5
        self.momtwist = momtwist
        self.timer = qtimer
        assert isinstance(self.mass, float)
        assert isinstance(self.m5, float)
        assert isinstance(self.momtwist, CoordinateD)
        assert isinstance(self.timer, (Timer, TimerNone,))

    def __mul__(self, prop_src):
        """
        prop_src: prop or [ prop, ... ]
        """
        if isinstance(prop_src, Prop):
            self.timer.start()
            prop_sol = free_invert(prop_src, self.mass, self.m5, self.momtwist)
            self.timer.stop()
            return prop_sol
        elif isinstance(prop_src, list):
            return [self * p for p in prop_src]
        else:
            raise Exception("InverterDwfFreeField")

## -----

class InverterDomainWall(Inverter):
    """
    self.cdata
    self.timer
    """

    def __init__(self, *, gf, fa, qtimer=TimerNone()):
        self.cdata = mk_inverter_domain_wall(gf, fa)
        self.timer = qtimer
        assert isinstance(self.timer, (Timer, TimerNone,))

    def __del__(self):
        assert isinstance(self.cdata, int)
        free_inverter_domain_wall(self)

    def __mul__(self, prop_src):
        """
        prop_src: prop or [ prop, ... ]
        """
        if isinstance(prop_src, Prop):
            self.timer.start()
            prop_sol = Prop()
            invert_inverter_domain_wall(prop_sol, prop_src, self)
            self.timer.stop()
            return prop_sol
        elif isinstance(prop_src, list):
            return [self * p for p in prop_src]
        else:
            raise Exception("InverterDomainWall")

    def stop_rsd(self):
        return get_inverter_domain_wall_ptr(self).stop_rsd()

    def set_stop_rsd(self, stop_rsd):
        cc.py_set_stop_rsd_inverter_domain_wall(
            get_inverter_domain_wall_ptr(self)[0], stop_rsd)

    def max_num_iter(self):
        return get_inverter_domain_wall_ptr(self).max_num_iter()

    def set_max_num_iter(self, max_num_iter):
        cc.py_set_max_num_iter_inverter_domain_wall(
            get_inverter_domain_wall_ptr(self)[0], max_num_iter)

    def max_mixed_precision_cycle(self):
        return get_inverter_domain_wall_ptr(self).max_mixed_precision_cycle()

    def set_max_mixed_precision_cycle(self, max_mixed_precision_cycle):
        cc.py_set_max_mixed_precision_cycle_inverter_domain_wall(
            get_inverter_domain_wall_ptr(self)[0], max_mixed_precision_cycle)

## -----

class InverterGaugeTransform(Inverter):
    """
    self.inverter
    self.gt
    self.gt_inv
    self.timer
    """

    def __init__(
        self,
        *,
        inverter,
        gt,
        qtimer=TimerNone(),
    ):
        self.inverter = inverter
        self.gt = gt
        self.timer = qtimer
        assert isinstance(self.inverter, Inverter)
        assert isinstance(self.gt, GaugeTransform)
        assert isinstance(self.timer, (Timer, TimerNone,))
        self.gt_inv = self.gt.inv()

    def __mul__(self, prop_src):
        assert isinstance(prop_src, (Prop, FermionField4d, list,))
        self.timer.start()
        src = self.gt_inv * prop_src
        sol = self.inverter * src
        prop_sol = self.gt * sol
        self.timer.stop()
        return prop_sol

## -----

class EigSystem:
    pass

## -----
