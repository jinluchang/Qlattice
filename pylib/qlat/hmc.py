import cqlat as c

from qlat.field import *
from qlat.rng import *
from qlat.gauge_action import *

def set_rand_gauge_momentum(gm, sigma, rng):
    assert type(gm) == Field and gm.ctype == "ColorMatrix"
    assert type(sigma) == float
    assert type(rng) == RngState
    return c.set_rand_gauge_momentum(gm, sigma, rng)

def gm_hamilton_node(gm):
    assert type(gm) == Field and gm.ctype == "ColorMatrix"
    return c.gm_hamilton_node(gm)

def gf_hamilton_node(gf, ga):
    assert type(gf) == Field and gf.ctype == "ColorMatrix"
    assert type(ga) == GaugeAction
    return c.gf_hamilton_node(gf, ga)

def gf_evolve(gf, gm, step_size):
    assert type(gm) == Field and gm.ctype == "ColorMatrix"
    return c.gf_evolve(gf, gm, step_size)

def set_gm_force(gm_force, gf, ga):
    assert type(gm_force) == Field and gm_force.ctype == "ColorMatrix"
    assert type(gf) == Field and gf.ctype == "ColorMatrix"
    assert type(ga) == GaugeAction
    return c.set_gm_force(gm_force, gf, ga)

