import cqlat as c

from qlat.field import *
from qlat.qcd import *
from qlat.rng import *
from qlat.gauge_action import *

class GaugeMomentum(Field):
    def __init__(self, geo = None):
        Field.__init__(self, "ColorMatrix", geo, 4)

def set_rand_gauge_momentum(gm, sigma, rng):
    assert type(gm) == GaugeMomentum
    assert type(sigma) == float
    assert type(rng) == RngState
    return c.set_rand_gauge_momentum(gm, sigma, rng)

def gm_hamilton_node(gm):
    assert type(gm) == GaugeMomentum
    return c.gm_hamilton_node(gm)

def gf_hamilton_node(gf, ga):
    assert type(gf) == GaugeField
    assert type(ga) == GaugeAction
    return c.gf_hamilton_node(gf, ga)

def gf_evolve(gf, gm, step_size):
    assert type(gm) == GaugeMomentum
    return c.gf_evolve(gf, gm, step_size)

def set_gm_force(gm_force, gf, ga):
    assert type(gm_force) == GaugeMomentum
    assert type(gf) == GaugeField
    assert type(ga) == GaugeAction
    return c.set_gm_force(gm_force, gf, ga)

