from qlat_utils import *
from .c import *
from . import c

@timer
def contract_chvp_16(prop1, prop2):
    chvp_16 = Field(ElemTypeComplexD)
    c.contract_chvp_16_field(chvp_16, prop1, prop2)
    return chvp_16
