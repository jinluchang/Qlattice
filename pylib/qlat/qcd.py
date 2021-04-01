import cqlat as c

from qlat.field import *
from qlat.rng_state import *
from qlat.timer import *

class GaugeField(Field):

    def __init__(self, geo = None):
        Field.__init__(self, "ColorMatrix", geo, 4)

    def save(self, path):
        return c.save_gauge_field(self, path)

    def load(self, path):
        return c.load_gauge_field(self, path)

    def set_rand(self, rng, sigma = 0.5, n_step = 1):
        set_g_rand_color_matrix_field(self, rng, sigma, n_step)

    def unitarize(self):
        c.unitarize_color_matrix_field(self)

    def plaq(self):
        return gf_avg_plaq(self)

    def link_trace(self):
        return gf_avg_link_trace(self)

    def show_info(self):
        gf_show_info(self)

class GaugeTransform(Field):

    def __init__(self, geo = None):
        Field.__init__(self, "ColorMatrix", geo, 1)

    def set_rand(self, rng, sigma = 0.5, n_step = 1):
        set_g_rand_color_matrix_field(self, rng, sigma, n_step)

    def unitarize(self):
        c.unitarize_color_matrix_field(self)

    def __mul__(self, other):
        from qlat.propagator import Propagator4d, Prop, SelProp, PselProp
        if isinstance(other, GaugeTransform):
            gt = GaugeTransform()
            c.apply_gt_gt(gt, self, other)
            return gt
        elif isinstance(other, GaugeField):
            gf = GaugeField()
            c.apply_gt_gf(gf, self, other)
            return gf
        elif isinstance(other, Propagator4d):
            prop = Prop()
            c.apply_gt_prop(prop, self, other)
            return prop
        elif isinstance(other, SelProp):
            prop = SelProp(other.fsel)
            c.apply_gt_sprop(prop, self, other)
            return prop
        elif isinstance(other, PselProp):
            prop = PselProp(other.psel)
            c.apply_gt_psprop(prop, self, other)
            return prop
        elif isinstance(other, list):
            return [ self * p for p in other ]
        else:
            raise Exception("GaugeTransform.__mul__")

    def inv(self):
        gt = GaugeTransform()
        c.gt_invert(gt, self)
        return gt

@timer_verbose
def gf_show_info(gf):
    assert isinstance(gf, GaugeField)
    displayln_info(f"gf_show_info: plaq = {gf.plaq():.16F} ; link_trace = {gf.link_trace():.16F}.")

def gf_avg_plaq(gf):
    assert isinstance(gf, GaugeField)
    return c.gf_avg_plaq(gf)

def gf_avg_link_trace(gf):
    assert isinstance(gf, GaugeField)
    return c.gf_avg_link_trace(gf)

def set_g_rand_color_matrix_field(fc, rng, sigma, n_steps = 1):
    assert isinstance(fc, Field) and fc.ctype == "ColorMatrix"
    assert isinstance(rng, RngState)
    return c.set_g_rand_color_matrix_field(fc, rng, sigma, n_steps)

def gf_twist_boundary_at_boundary(gf, lmom, mu):
    return c.gf_twist_boundary_at_boundary(gf, lmom, mu)
