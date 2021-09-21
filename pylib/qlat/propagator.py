import cqlat as c

from qlat.timer import *
from qlat.field import *
from qlat.field_utils import *
from qlat.rng_state import *
from qlat.selected_field import *
from qlat.selected_points import *

class Prop(Field):

    def __init__(self, geo = None):
        Field.__init__(self, "WilsonMatrix", geo, 1)

    def copy(self, is_copying_data = True):
        f = Prop()
        if is_copying_data:
            f @= self
        return f

    def sparse(self, sel):
        if isinstance(sel, PointSelection):
            psel = sel
            sp = PselProp(psel)
            set_selected_points(sp, self)
            return sp
        elif isinstance(sel, FieldSelection):
            fsel = sel
            sf = SelProp(fsel)
            set_selected_field(sf, self)
            return sf
        else:
            raise Exception("Prop.sparse")

class SelProp(SelectedField):

    def __init__(self, fsel):
        SelectedField.__init__(self, "WilsonMatrix", fsel, 1)

    def copy(self, is_copying_data = True):
        f = SelProp(self.fsel)
        if is_copying_data:
            f @= self
        return f

    def sparse(self, sel):
        if isinstance(sel, PointSelection):
            psel = sel
            sp = PselProp(psel)
            set_selected_points(sp, self)
            return sp
        elif isinstance(sel, FieldSelection):
            fsel = sel
            sf = SelProp(fsel)
            set_selected_field(sf, self)
            return sf
        else:
            raise Exception("SelProp.sparse")

class PselProp(SelectedPoints):

    def __init__(self, psel):
        SelectedPoints.__init__(self, "WilsonMatrix", psel, 1)

    def copy(self, is_copying_data = True):
        f = PselProp(self.psel)
        if is_copying_data:
            f @= self
        return f

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

@timer
def free_invert(prop_src, mass,
        m5 = 1.0, momtwist = [0.0, 0.0, 0.0, 0.0]):
    assert isinstance(prop_src, Prop)
    prop_sol = Prop()
    c.free_invert_prop(prop_sol, prop_src, mass, m5, momtwist)
    return prop_sol

def convert_mspincolor_from_wm_prop(prop_msc, prop_wm):
    assert isinstance(prop_msc, Prop)
    assert isinstance(prop_wm, Prop)
    return c.convert_mspincolor_from_wm_prop(prop_msc, prop_wm)

def convert_wm_from_mspincolor_prop(prop_wm, prop_msc):
    assert isinstance(prop_wm, Prop)
    assert isinstance(prop_msc, Prop)
    return c.convert_wm_from_mspincolor_prop(prop_wm, prop_msc)

@timer
def free_scalar_invert_mom_cfield(f, mass):
    assert isinstance(f, Field)
    assert f.ctype == "Complex"
    c.free_scalar_invert_mom_cfield(f, mass)

@timer
def free_scalar_invert_cfield(src, mass):
    fft_f = mk_fft(True, is_normalizing = True)
    fft_b = mk_fft(False, is_normalizing = True)
    f = fft_f * src
    free_scalar_invert_mom_cfield(f, mass)
    sol = fft_b * f
    return sol

class FermionField4d(Field):

    def __init__(self, geo = None):
        Field.__init__(self, "WilsonVector", geo, 1)

    def copy(self, is_copying_data = True):
        f = FermionField4d()
        if is_copying_data:
            f @= self
        return f
