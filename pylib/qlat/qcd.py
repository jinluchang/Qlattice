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

@timer_verbose
def gf_show_info(gf):
    assert isinstance(gf, GaugeField)
    plaq = gf.plaq()
    link_trace = gf.link_trace()
    displayln_info(f"gf_show_info: plaq = {plaq:.16F} ; link_trace = {link_trace:.16F}.")

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

def gf_twist_boundary_at_boundary(gf, mom, mu):
    return c.gf_twist_boundary_at_boundary(gf, mom, mu)
