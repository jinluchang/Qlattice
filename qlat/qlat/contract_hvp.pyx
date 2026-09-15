# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.contract_hvp``
=============================\n
Contract a conserved-point hadronic vacuum polarization (HVP) on a
single time slice between two propagators.  Returns the result as a
``LatData`` object.\n
Documentation: ``docs/qlat/qlat_contract_hvp.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils.all cimport *
from . cimport everything as cc
from .propagator cimport SelProp

import qlat_utils as q

@q.timer
def contract_chvp3_field(prop1, prop2, tslice):
    cdef LatData ld = LatData()
    ld.xx = cc.contract_chvp3((<SelProp>prop1).xxx().val(),
                              (<SelProp>prop2).xxx().val(), tslice,
                              (<SelProp>prop1).fsel.xx)
    return ld
