from qlat_utils import *
import qlat.c as c

from qlat.field import *
from qlat.qcd import *
from qlat.gauge_action import *
from qlat.elem_type import *

class GaugeMomentum(FieldColorMatrix):

    def __init__(self, geo = None):
        super().__init__(geo, 4)

    def set_rand(self, rng, sigma = 1.0):
        set_rand_gauge_momentum(self, sigma, rng)

def set_rand_gauge_momentum(gm, sigma, rng):
    assert isinstance(gm, GaugeMomentum)
    assert isinstance(sigma, float)
    assert isinstance(rng, RngState)
    return c.set_rand_gauge_momentum(gm, sigma, rng)

def gm_hamilton_node(gm):
    assert isinstance(gm, GaugeMomentum)
    return c.gm_hamilton_node(gm)

def gf_hamilton_node(gf, ga):
    assert isinstance(gf, GaugeField)
    assert isinstance(ga, GaugeAction)
    return c.gf_hamilton_node(gf, ga)

def gf_evolve(gf, gm, step_size):
    assert isinstance(gm, GaugeMomentum)
    return c.gf_evolve(gf, gm, step_size)

def set_gm_force(gm_force, gf, ga):
    assert isinstance(gm_force, GaugeMomentum)
    assert isinstance(gf, GaugeField)
    assert isinstance(ga, GaugeAction)
    return c.set_gm_force(gm_force, gf, ga)

