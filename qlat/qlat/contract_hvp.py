"""
Module ``qlat.contract_hvp``
=============================\n
Contract a conserved-point hadronic vacuum polarization (HVP) on a
single time slice between two propagators.  Returns the result as a
``LatData`` object.\n
Documentation: ``docs/qlat/qlat_contract_hvp.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils import *
from .c import *
from . import c

@timer
def contract_chvp3_field(prop1, prop2, tslice):
    ld = LatData()
    c.contract_chvp3_sfield(ld, prop1, prop2, tslice)
    return ld
