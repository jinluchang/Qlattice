import cqlat as c

from qlat.field import *
from qlat.qcd import *
from qlat.rng_state import *
from qlat.gauge_action import *
from qlat.scalar_action import *

class GaugeMomentum(Field):

    def __init__(self, geo = None, *, ctype = None, multiplicity = None):
        Field.__init__(self, "ColorMatrix", geo, 4)

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

def sm_hamilton_node(sm):
    assert isinstance(sm, Field)
    return c.sm_hamilton_node(sm)

def sf_hamilton_node(sf, sa):
    assert isinstance(sf, Field)
    assert isinstance(sa, ScalarAction)
    return c.sf_hamilton_node(sf, sa)
    
def sf_evolve(sf, sm, step_size):
    assert isinstance(sf, Field)
    assert isinstance(sm, Field)
    return c.sf_evolve(sf, sm, step_size)

def set_sm_force(sm_force, sf, sa):
    assert isinstance(sm_force, Field)
    assert isinstance(sf, Field)
    assert isinstance(sa, ScalarAction)
    return c.set_sm_force(sm_force, sf, sa)

