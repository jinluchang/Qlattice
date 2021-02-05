import cqlat as c

from qlat.field import *
from qlat.rng import *

def gauge_field(geo):
    return Field("ColorMatrix", geo, 4)

def gf_show_info(gf):
    assert type(gf) == Field and gf.ctype == "ColorMatrix"
    c.gf_show_info(gf)

def gf_avg_plaq(gf):
    assert type(gf) == Field and gf.ctype == "ColorMatrix"
    return c.gf_avg_plaq(gf)

def gf_avg_link_trace(gf):
    assert type(gf) == Field and gf.ctype == "ColorMatrix"
    return c.gf_avg_link_trace(gf)

def set_g_rand_color_matrix_field(fc, rng, sigma, n_steps = 1):
    assert type(fc) == Field and fc.ctype == "ColorMatrix"
    assert type(rng) == RngState
    c.set_g_rand_color_matrix_field(fc, rng, sigma, n_steps)

def unitarize(x):
    if type(x) == Field and x.ctype == "ColorMatrix":
        c.unitarize_color_matrix_field(x)
    else:
        raise Exception("unitarize")

