# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.fthmc``
=====================\n
Flow-time Hamiltonian Monte Carlo (ftHMC) utilities: gradient flow
evolution, flowed Hamiltonian evaluation, and flowed force computation.
``FlowInfo`` is a ``cdef class`` that owns its C++ object by value.\n
Documentation: ``docs/qlat/qlat_fthmc.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils.all cimport *
from . cimport everything as cc
from .qcd cimport GaugeField
from .hmc cimport GaugeMomentum
from .gauge_action cimport GaugeAction

cdef class FlowInfo:

    cdef cc.FlowInfo xx

    def __cinit__(self):
        self.xx = cc.FlowInfo()

    def add_flow(self, eo, mu, epsilon, flow_size=1):
        self.xx.v.push_back(cc.FlowStepInfo(eo, mu, epsilon, flow_size))

    def add_rand_order_flow(self, rng, epsilon, *args):
        if len(args) == 0:
            add_rand_order_flow_flow_info(self, rng, epsilon)
        elif len(args) == 1:
            epsilon2 = args[0]
            add_rand_order_flow2_flow_info(self, rng, epsilon, epsilon2)
        else:
            raise Exception("add_rand_order_flow")

    def show(self):
        return cc.show(self.xx)

### -------------------------------------------------------------------

def add_flow_flow_info(FlowInfo fi, eo, mu, epsilon, flow_size=1):
    fi.add_flow(eo, mu, epsilon, flow_size)

def show_flow_info(FlowInfo fi):
    return fi.show()

def add_rand_order_flow_flow_info(FlowInfo fi, RngState rs, cc.RealD epsilon):
    cdef cc.FlowInfo step = cc.mk_flow_info_step(rs.xx, epsilon)
    cdef cc.Long i
    for i in range(step.v.size()):
        fi.xx.v.push_back(step.v[i])

def add_rand_order_flow2_flow_info(
        FlowInfo fi, RngState rs, cc.RealD epsilon, cc.RealD epsilon2):
    cdef cc.FlowInfo step = cc.mk_flow_info_step(rs.xx, epsilon, epsilon2)
    cdef cc.Long i
    for i in range(step.v.size()):
        fi.xx.v.push_back(step.v[i])

def gf_flow(GaugeField gf, GaugeField gf0, FlowInfo fi):
    cc.gf_flow(gf.xxx().val(), gf0.xxx().val(), fi.xx)

def gf_flow_inv(GaugeField gf, GaugeField gf1, FlowInfo fi):
    cc.gf_flow_inv(gf.xxx().val(), gf1.xxx().val(), fi.xx)

def gf_hamilton_flowed_node(GaugeField gf0, GaugeAction ga, FlowInfo fi):
    return cc.gf_hamilton_flowed_node(gf0.xxx().val(), ga.xx, fi.xx)

def set_gm_force_flowed(
        GaugeMomentum gm_force, GaugeField gf0, GaugeAction ga, FlowInfo fi):
    cc.set_gm_force_flowed(
        gm_force.xxx().val(), gf0.xxx().val(), ga.xx, fi.xx)

def set_gm_force_flowed_no_det(
        GaugeMomentum gm_force, GaugeMomentum gm_force_pre, GaugeField gf0,
        FlowInfo fi):
    cc.set_gm_force_flowed_no_det(
        gm_force.xxx().val(), gm_force_pre.xxx().val(), gf0.xxx().val(),
        fi.xx)
