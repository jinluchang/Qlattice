# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.fthmc``
=====================\n
Flow-time Hamiltonian Monte Carlo (ftHMC) utilities: gradient flow
evolution, flowed Hamiltonian evaluation, and flowed force computation.\n
Documentation: ``docs/qlat/qlat_fthmc.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils.all cimport *
from . cimport everything as cc
from .qcd cimport GaugeField
from .hmc cimport GaugeMomentum
from .gauge_action cimport GaugeAction

from cpython.long cimport PyLong_FromVoidPtr
from cpython.long cimport PyLong_AsVoidPtr

cdef inline cc.FlowInfo* get_flow_info_ptr(object fi) except? NULL:
    return <cc.FlowInfo*>PyLong_AsVoidPtr(fi.cdata)

def free_flow_info(fi):
    cdef cc.FlowInfo* pfi = get_flow_info_ptr(fi)
    del pfi

def add_flow_flow_info(fi, eo, mu, epsilon, flow_size=1):
    cdef cc.FlowInfo* pfi = get_flow_info_ptr(fi)
    pfi.v.push_back(cc.FlowStepInfo(eo, mu, epsilon, flow_size))

class FlowInfo:

    def __init__(self):
        self.cdata = mk_flow_info()

    def __del__(self):
        free_flow_info(self)

    def add_flow(self, eo, mu, epsilon, flow_size=1):
        add_flow_flow_info(self, eo, mu, epsilon, flow_size)

    def add_rand_order_flow(self, rng, epsilon, *args):
        if len(args) == 0:
            add_rand_order_flow_flow_info(self, rng, epsilon)
        elif len(args) == 1:
            epsilon2 = args[0]
            add_rand_order_flow2_flow_info(self, rng, epsilon, epsilon2)
        else:
            raise Exception("add_rand_order_flow")

    def show(self):
        return show_flow_info(self)

def mk_flow_info():
    cdef cc.FlowInfo* pfi = new cc.FlowInfo()
    return PyLong_FromVoidPtr(<void*>pfi)

def show_flow_info(fi):
    cdef cc.FlowInfo* pfi = get_flow_info_ptr(fi)
    return cc.show(pfi[0])

def add_rand_order_flow_flow_info(fi, RngState rs, cc.RealD epsilon):
    cdef cc.FlowInfo* pfi = get_flow_info_ptr(fi)
    cdef cc.FlowInfo step = cc.mk_flow_info_step(rs.xx, epsilon)
    cdef cc.Long i
    for i in range(step.v.size()):
        pfi.v.push_back(step.v[i])

def add_rand_order_flow2_flow_info(
        fi, RngState rs, cc.RealD epsilon, cc.RealD epsilon2):
    cdef cc.FlowInfo* pfi = get_flow_info_ptr(fi)
    cdef cc.FlowInfo step = cc.mk_flow_info_step(rs.xx, epsilon, epsilon2)
    cdef cc.Long i
    for i in range(step.v.size()):
        pfi.v.push_back(step.v[i])

def gf_flow(GaugeField gf, GaugeField gf0, fi):
    cdef cc.FlowInfo* pfi = get_flow_info_ptr(fi)
    assert isinstance(fi, FlowInfo)
    cc.gf_flow(gf.xxx().val(), gf0.xxx().val(), pfi[0])

def gf_flow_inv(GaugeField gf, GaugeField gf1, fi):
    cdef cc.FlowInfo* pfi = get_flow_info_ptr(fi)
    assert isinstance(fi, FlowInfo)
    cc.gf_flow_inv(gf.xxx().val(), gf1.xxx().val(), pfi[0])

def gf_hamilton_flowed_node(GaugeField gf0, GaugeAction ga, fi):
    cdef cc.FlowInfo* pfi = get_flow_info_ptr(fi)
    assert isinstance(fi, FlowInfo)
    return cc.gf_hamilton_flowed_node(gf0.xxx().val(), ga.xx, pfi[0])

def set_gm_force_flowed(
        GaugeMomentum gm_force, GaugeField gf0, GaugeAction ga, fi):
    cdef cc.FlowInfo* pfi = get_flow_info_ptr(fi)
    assert isinstance(fi, FlowInfo)
    cc.set_gm_force_flowed(
        gm_force.xxx().val(), gf0.xxx().val(), ga.xx, pfi[0])

def set_gm_force_flowed_no_det(
        GaugeMomentum gm_force, GaugeMomentum gm_force_pre, GaugeField gf0, fi):
    cdef cc.FlowInfo* pfi = get_flow_info_ptr(fi)
    assert isinstance(fi, FlowInfo)
    cc.set_gm_force_flowed_no_det(
        gm_force.xxx().val(), gm_force_pre.xxx().val(), gf0.xxx().val(),
        pfi[0])
