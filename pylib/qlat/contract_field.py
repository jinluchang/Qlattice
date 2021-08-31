import cqlat as c

from qlat.timer import *
from qlat.lat_io import *
from qlat.propagator import *

@timer
def contract_chvp_16(prop1, prop2):
    chvp_16 = q.Field("Complex")
    c.contract_chvp_16_field(chvp_16, prop1, prop2)
    return chvp_16
