"""
Module ``qlat.contract_pion``
===============================\n
Compute the pion two-point correlation function from a propagator on a
single time slice.  Supports both dense ``Prop`` and selected
``SelProp`` inputs.\n
Documentation: ``docs/qlat/qlat_contract_pion.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils import *
from .c import *
from . import c

@timer
def contract_pion_field(prop, tslice):
    ld = LatData()
    if isinstance(prop, Prop):
        c.contract_pion_field(ld, prop, tslice)
    elif isinstance(prop, SelProp):
        c.contract_pion_sfield(ld, prop, tslice)
    else:
        raise Exception("contract_pion_field")
    return ld
