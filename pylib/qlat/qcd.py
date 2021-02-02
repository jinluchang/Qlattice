import cqlat as c

from qlat.field import *
from qlat.rng import *

def gf_show_info(x):
    if type(x) == Field and x.ctype == "ColorMatrix":
        c.gf_show_info(x.ctype, x.cdata)
    else:
        raise Exception("gf_show_info")

def set_g_rand_color_matrix_field(x, rng, sigma, n_steps = 1):
    if type(x) == Field and x.ctype == "ColorMatrix" and type(rng) == RngState:
        c.set_g_rand_color_matrix_field(x.ctype, x.cdata, rng.cdata, sigma, n_steps)
    else:
        raise Exception("set_g_rand_color_matrix_field")

def unitarize(x):
    if type(x) == Field and x.ctype == "ColorMatrix":
        c.unitarize_color_matrix_field(x.ctype, x.cdata)
    else:
        raise Exception("unitarize")

