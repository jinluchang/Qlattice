import cqlat as c

from qlat.field import *
from qlat.rng_state import *

from qlat.selected_field import *
from qlat.selected_points import *

class Propagator4d(Field):

    def __init__(self, geo = None):
        Field.__init__(self, "WilsonMatrix", geo, 1)

class Prop(Propagator4d):

    pass

class SelProp(SelectedField):

    def __init__(self, fsel):
        SelectedField.__init__(self, "WilsonMatrix", fsel, 1)

class PselProp(SelectedPoints):

    def __init__(self, psel):
        SelectedPoints.__init__(self, "WilsonMatrix", psel, 1)

def set_point_src(prop_src, geo, xg, value = 1.0):
    c.set_point_src_prop(prop_src, geo, xg, value)

def set_wall_src(prop_src, geo, tslice, lmom = [0.0, 0.0, 0.0, 0.0]):
    c.set_wall_src_prop(prop_src, geo, tslice, lmom)

def mk_point_src(geo, xg, value = 1.0):
    prop_src = Prop()
    set_point_src(prop_src, geo, xg, value);
    return prop_src

def mk_wall_src(geo, tslice, lmom = [0.0, 0.0, 0.0, 0.0]):
    prop_src = Prop()
    set_wall_src(prop_src, geo, tslice, lmom);
    return prop_src

def free_invert(prop_src, mass,
        m5 = 1.0, momtwist = [0.0, 0.0, 0.0, 0.0]):
    assert isinstance(prop_src, Propagator4d)
    prop_sol = Prop()
    c.free_invert_prop(prop_sol, prop_src, mass, m5, momtwist)
    return prop_sol

def convert_mspincolor_from_wm_prop(prop_msc, prop_wm):
    assert isinstance(prop_msc, Propagator4d)
    assert isinstance(prop_wm, Propagator4d)
    return c.convert_mspincolor_from_wm_prop(prop_msc, prop_wm)

def convert_wm_from_mspincolor_prop(prop_wm, prop_msc):
    assert isinstance(prop_wm, Propagator4d)
    assert isinstance(prop_msc, Propagator4d)
    return c.convert_wm_from_mspincolor_prop(prop_wm, prop_msc)
