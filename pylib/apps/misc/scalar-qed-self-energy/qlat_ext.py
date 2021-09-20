from qlat import *

import qlat as q
import cqlat as c

@q.timer
def free_scalar_invert_mom_cfield(f, mass):
    assert isinstance(f, q.Field)
    assert f.ctype == "Complex"
    c.free_scalar_invert_mom_cfield(f, mass)

@q.timer
def free_scalar_invert_cfield(src, mass):
    fft_f = mk_fft(True, is_normalizing = True)
    fft_b = mk_fft(False, is_normalizing = True)
    f = fft_f * src
    free_scalar_invert_mom_cfield(f, mass)
    sol = fft_b * f
    return sol
