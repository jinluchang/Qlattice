# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.inverter``
=========================\n
Fermion-matrix inverter framework for domain-wall fermions.
``InverterDomainWall`` is a ``cdef class`` that owns its C++ object by value.
The other inverter classes (``InverterDwfFreeField``,
``InverterGaugeTransform``, ``EigSystem``) are pure Python and live in
``qlat.inverter_utils``.\n
Documentation: ``docs/qlat/qlat_inverter.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils.all cimport *
from qlat_utils import mk_cache
from . cimport everything as cc
from .qcd cimport GaugeField
from .propagator cimport Prop
from .fermion_action cimport FermionAction

cache_inv = mk_cache("inv")

cdef class Inverter:
    pass

## -----

cdef class InverterDomainWall(Inverter):
    """
    self.xx
    self.timer
    """

    cdef cc.InverterDomainWall xx
    cdef public object timer

    def __cinit__(self):
        self.xx = cc.InverterDomainWall()

    def __init__(self, *, GaugeField gf, FermionAction fa, qtimer=TimerNone()):
        cc.setup_inverter(self.xx, gf.xxx().val(), fa.xx)
        self.timer = qtimer
        assert isinstance(self.timer, (Timer, TimerNone,))

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
        return self.xx.stop_rsd()

    def set_stop_rsd(self, stop_rsd):
        cc.py_set_stop_rsd_inverter_domain_wall(self.xx, stop_rsd)

    def max_num_iter(self):
        return self.xx.max_num_iter()

    def set_max_num_iter(self, max_num_iter):
        cc.py_set_max_num_iter_inverter_domain_wall(self.xx, max_num_iter)

    def max_mixed_precision_cycle(self):
        return self.xx.max_mixed_precision_cycle()

    def set_max_mixed_precision_cycle(self, max_mixed_precision_cycle):
        cc.py_set_max_mixed_precision_cycle_inverter_domain_wall(
            self.xx, max_mixed_precision_cycle)

### -------------------------------------------------------------------
### cqlat-compatible entry points

def invert_inverter_domain_wall(
        Prop prop_sol, Prop prop_src, InverterDomainWall inv):
    cc.invert_prop_dw(prop_sol.xxx().val(), prop_src.xxx().val(), inv.xx)

def get_stop_rsd_inverter_domain_wall(InverterDomainWall inv):
    return inv.stop_rsd()

def set_stop_rsd_inverter_domain_wall(InverterDomainWall inv, stop_rsd):
    inv.set_stop_rsd(stop_rsd)

def get_max_num_iter_inverter_domain_wall(InverterDomainWall inv):
    return inv.max_num_iter()

def set_max_num_iter_inverter_domain_wall(InverterDomainWall inv, max_num_iter):
    inv.set_max_num_iter(max_num_iter)

def get_max_mixed_precision_cycle_inverter_domain_wall(InverterDomainWall inv):
    return inv.max_mixed_precision_cycle()

def set_max_mixed_precision_cycle_inverter_domain_wall(
        InverterDomainWall inv, max_mixed_precision_cycle):
    inv.set_max_mixed_precision_cycle(max_mixed_precision_cycle)
