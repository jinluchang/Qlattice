from qlat import *

import qlat.c as c

class FlowInfo:

    def __init__(self):
        self.cdata = c.mk_flow_info()

    def __del__(self):
        c.free_flow_info(self)

    def add_flow(self, eo, mu, epsilon, flow_size = 1):
        c.add_flow_flow_info(self, eo, mu, epsilon, flow_size)

    def add_rand_order_flow(self, rng, epsilon, *args):
        if len(args) == 0:
            c.add_rand_order_flow_flow_info(self, rng, epsilon)
        elif len(args) == 1:
            epsilon2 = args[0]
            c.add_rand_order_flow2_flow_info(self, rng, epsilon, epsilon2)
        else:
            raise Exception("add_rand_order_flow")

    def show(self):
        return c.show_flow_info(self)

def gf_flow(gf, gf0, fi):
    assert isinstance(gf, GaugeField)
    assert isinstance(gf0, GaugeField)
    assert isinstance(fi, FlowInfo)
    c.gf_flow(gf, gf0, fi)

def gf_flow_inv(gf, gf1, fi):
    assert isinstance(gf, GaugeField)
    assert isinstance(gf1, GaugeField)
    assert isinstance(fi, FlowInfo)
    c.gf_flow_inv(gf, gf1, fi)

def gf_hamilton_flowed_node(gf0, ga, fi):
    assert isinstance(gf0, GaugeField)
    assert isinstance(ga, GaugeAction)
    assert isinstance(fi, FlowInfo)
    return c.gf_hamilton_flowed_node(gf0, ga, fi)

def set_gm_force_flowed(gm_force, gf0, ga, fi):
    assert isinstance(gm_force, GaugeMomentum)
    assert isinstance(gf0, GaugeField)
    assert isinstance(ga, GaugeAction)
    assert isinstance(fi, FlowInfo)
    c.set_gm_force_flowed(gm_force, gf0, ga, fi)

def set_gm_force_flowed_no_det(gm_force, gm_force_pre, gf0, fi):
    assert isinstance(gm_force, GaugeMomentum)
    assert isinstance(gm_force_pre, GaugeMomentum)
    assert isinstance(gf0, GaugeField)
    assert isinstance(fi, FlowInfo)
    c.set_gm_force_flowed_no_det(gm_force, gm_force_pre, gf0, fi)
