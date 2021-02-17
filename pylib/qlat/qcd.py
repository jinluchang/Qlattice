import cqlat as c

from qlat.field import *
from qlat.rng import *

class GaugeField(Field):

    def __init__(self, geo = None):
        Field.__init__(self, "ColorMatrix", geo, 4)

    def save(self, path):
        c.save_gauge_field(self, path)

    def load(self, path):
        c.load_gauge_field(self, path)

    def set_rand(self, rng, sigma = 0.5, n_step = 1):
        set_g_rand_color_matrix_field(self, rng, sigma, n_step)

    def unitarize(self):
        c.unitarize_color_matrix_field(self)

def gf_show_info(gf):
    assert isinstance(gf, GaugeField)
    c.gf_show_info(gf)

def gf_avg_plaq(gf):
    assert isinstance(gf, GaugeField)
    return c.gf_avg_plaq(gf)

def gf_avg_link_trace(gf):
    assert isinstance(gf, GaugeField)
    return c.gf_avg_link_trace(gf)

def set_g_rand_color_matrix_field(fc, rng, sigma, n_steps = 1):
    assert isinstance(fc, Field) and fc.ctype == "ColorMatrix"
    assert isinstance(rng, RngState)
    c.set_g_rand_color_matrix_field(fc, rng, sigma, n_steps)
