import qlat.cqlat as c

from qlat_utils import *

from qlat.field import *
from qlat.propagator import *

@timer
def contract_chvp_16(prop1, prop2):
    chvp_16 = Field("Complex")
    c.contract_chvp_16_field(chvp_16, prop1, prop2)
    return chvp_16
