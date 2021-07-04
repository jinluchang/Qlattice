import cqlat as c

from qlat.field import *

def mk_phase_field(geo: Geometry, lmom):
    # lmom is in lattice momentum unit
    # exp(i * 2*pi/L * lmom \cdot xg )
    f = Field("Complex", geo, 1)
    c.set_phase_field(f, lmom)
    return f
