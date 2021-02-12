import cqlat as c

from qlat.field import *
from qlat.rng import *

class GaugeField(Field):
    def __init__(self, geo = None):
        Field.__init__(self, "ColorMatrix", geo, 4)

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

def unitarize(x):
    if isinstance(x, Field) and x.ctype == "ColorMatrix":
        c.unitarize_color_matrix_field(x)
    else:
        raise Exception("unitarize")

