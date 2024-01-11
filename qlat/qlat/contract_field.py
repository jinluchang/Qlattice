from qlat_utils import *
from .c import *
from . import c

@timer
def contract_chvp_16(prop1, prop2):
    """
    return chvp_16
    #
    inline void contract_chvp_16(
        FieldM<ComplexD, 16>& chvp_16,
        const Propagator4d& prop1_x_y,
        const Propagator4d& prop2_x_y)
    #
    chvp_16.get_elem(x, mu * 4 + nu) ==
    tr(g5_herm(prop2_x_y.get_elem(x)) * gammas[mu]
    * prop1_x_y.get_elem(x) * gammas[nu])
    #
    mu: polarization at sink location x
    nu: polarization at source location y
    """
    chvp_16 = Field(ElemTypeComplexD)
    c.contract_chvp_16_field(chvp_16, prop1, prop2)
    return chvp_16
