import qlat.cqlat as c

from qlat_utils import *

from qlat.field import *
from qlat.field_utils import *
from qlat.selected_field import *
from qlat.selected_points import *
from qlat.mat import *

class Prop(Field):

    def __init__(self, geo = None, *, ctype = None, multiplicity = None):
        Field.__init__(self, "WilsonMatrix", geo, 1)

    def copy(self, is_copying_data = True):
        f = Prop()
        if is_copying_data:
            f @= self
        return f

    def get_elem_wm(self, xg, m = 0):
        cdata = c.get_elem_wm_prop(self, xg, m)
        wm = WilsonMatrix(cdata = cdata, base = self)
        return wm

    def sparse(self, sel):
        # deprecated
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

###

class SelProp(SelectedField):

    def __init__(self, fsel, *, ctype = None, multiplicity = None):
        SelectedField.__init__(self, "WilsonMatrix", fsel, 1)

    def copy(self, is_copying_data = True):
        f = SelProp(self.fsel)
        if is_copying_data:
            f @= self
        return f

    def get_elem_wm(self, idx, m = 0):
        cdata = c.get_elem_wm_sprop(self, idx, m)
        wm = WilsonMatrix(cdata = cdata, base = self)
        return wm

    def sparse(self, sel):
        # deprecated
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

###

class PselProp(SelectedPoints):

    def __init__(self, psel, *, ctype = None, multiplicity = None):
        SelectedPoints.__init__(self, "WilsonMatrix", psel, 1)

    def copy(self, is_copying_data = True):
        f = PselProp(self.psel)
        if is_copying_data:
            f @= self
        return f

    def get_elem_wm(self, idx, m = 0):
        cdata = c.get_elem_wm_psprop(self, idx, m)
        wm = WilsonMatrix(cdata = cdata, base = self)
        return wm

###

def set_point_src(prop_src, geo, xg, value = 1.0):
    c.set_point_src_prop(prop_src, geo, xg, value)

def set_wall_src(prop_src, geo, tslice, lmom = None):
    if lmom is None:
        lmom = [ 0.0, 0.0, 0.0, 0.0, ]
    c.set_wall_src_prop(prop_src, geo, tslice, lmom)

def mk_point_src(geo, xg, value = 1.0):
    prop_src = Prop()
    set_point_src(prop_src, geo, xg, value)
    return prop_src

def mk_wall_src(geo, tslice, lmom = None):
    if lmom is None:
        lmom = [ 0.0, 0.0, 0.0, 0.0, ]
    prop_src = Prop()
    set_wall_src(prop_src, geo, tslice, lmom)
    return prop_src

@timer
def mk_rand_u1_src(sel, rs):
    # return (prop_src, fu1,) where prop_src = Prop() and fu1 = Field("Complex")
    # fu1 stores the random u1 numbers (fu1.multiplicity() == 1)
    # sel can be psel or fsel
    prop_src = Prop()
    fu1 = Field("Complex")
    if isinstance(sel, FieldSelection):
        fsel = sel
        c.set_rand_u1_src_fsel(prop_src, fu1, fsel, rs)
    elif isinstance(sel, PointSelection):
        psel = sel
        geo = psel.geo
        assert isinstance(geo, Geometry)
        c.set_rand_u1_src_psel(prop_src, fu1, psel, geo, rs)
    else:
        raise Exception(f"mk_rand_u1_src {type(sel)}")
    return (prop_src, fu1,)

@timer
def get_rand_u1_sol(prop_sol, fu1, sel):
    assert isinstance(prop_sol, Prop)
    assert isinstance(fu1, Field) and fu1.ctype == "Complex"
    if isinstance(sel, FieldSelection):
        fsel = sel
        s_prop = SelProp(fsel)
        c.set_rand_u1_sol_fsel(s_prop, prop_sol, fu1, fsel)
        return s_prop
    elif isinstance(sel, PointSelection):
        psel = sel
        sp_prop = PselProp(psel)
        c.set_rand_u1_sol_psel(sp_prop, prop_sol, fu1, psel)
        return sp_prop
    else:
        raise Exception(f"get_rand_u1_sol {type(sel)}")

@timer_verbose
def mk_rand_u1_prop(inv, sel, rs):
    # interface function
    # return s_prop
    # sel can be psel or fsel
    prop_src, fu1 = mk_rand_u1_src(sel, rs)
    prop_sol = inv * prop_src
    return get_rand_u1_sol(prop_sol, fu1, sel)

@timer
def free_invert(prop_src, mass, m5 = 1.0, momtwist = None):
    assert isinstance(prop_src, Prop)
    if momtwist is None:
        momtwist = [ 0.0, 0.0, 0.0, 0.0, ]
    prop_sol = Prop()
    c.free_invert_prop(prop_sol, prop_src, mass, m5, momtwist)
    return prop_sol

def convert_mspincolor_from_wm(prop_wm):
    prop_msc = prop_wm.copy(False)
    if isinstance(prop_wm, Prop):
        c.convert_mspincolor_from_wm_prop(prop_msc, prop_wm)
    elif isinstance(prop_wm, SelProp):
        c.convert_mspincolor_from_wm_s_prop(prop_msc, prop_wm)
    elif isinstance(prop_wm, PselProp):
        c.convert_mspincolor_from_wm_sp_prop(prop_msc, prop_wm)
    else:
        raise Exception("prop type match failed")
    return prop_msc

def convert_wm_from_mspincolor(prop_msc):
    prop_wm = prop_msc.copy(False)
    if isinstance(prop_msc, Prop):
        c.convert_wm_from_mspincolor_prop(prop_wm, prop_msc)
    elif isinstance(prop_msc, SelProp):
        c.convert_wm_from_mspincolor_s_prop(prop_wm, prop_msc)
    elif isinstance(prop_msc, PselProp):
        c.convert_wm_from_mspincolor_sp_prop(prop_wm, prop_msc)
    else:
        raise Exception("prop type match failed")
    return prop_wm

@timer
def flip_tpbc_with_tslice(prop, tslice_flip_tpbc):
    if isinstance(prop, SelProp):
        c.flip_tpbc_with_tslice_s_prop(prop, tslice_flip_tpbc)
    elif isinstance(prop, PselProp):
        c.flip_tpbc_with_tslice_sp_prop(prop, tslice_flip_tpbc)
    else:
        print(type(prop))
        assert False

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

    def __init__(self, geo = None, *, ctype = None, multiplicity = None):
        Field.__init__(self, "WilsonVector", geo, 1)

    def copy(self, is_copying_data = True):
        f = FermionField4d()
        if is_copying_data:
            f @= self
        return f

###
