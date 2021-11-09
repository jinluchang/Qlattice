from qlat import *

import qlat as q
import cqlat as c

@q.timer
def mk_pion_four_point_field(total_site, pion_mass, tag = "", r_pi = 0.0):
    geo = q.Geometry(total_site, 16)
    f = q.Field("Complex", geo)
    c.set_pion_four_point_mom_field(f, pion_mass, tag, r_pi)
    fft_b = mk_fft(False)
    f = fft_b * f
    f *= 1 / geo.total_volume()
    return f

def acc_four_point_func_em(ld, field, dtype, r_scaling_factor):
    c.acc_four_point_func_em(ld, field, dtype, r_scaling_factor)
