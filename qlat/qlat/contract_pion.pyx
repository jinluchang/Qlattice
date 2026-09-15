# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.contract_pion``
===============================\n
Compute the pion two-point correlation function from a propagator on a
single time slice.  Supports both dense ``Prop`` and selected
``SelProp`` inputs.\n
Documentation: ``docs/qlat/qlat_contract_pion.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils.all cimport *
from . cimport everything as cc
from .propagator cimport Prop
from .propagator cimport SelProp

import qlat_utils as q

@q.timer
def contract_pion_field(prop, tslice):
    cdef LatData ld = LatData()
    if isinstance(prop, Prop):
        ld.xx = cc.contract_pion((<Prop>prop).xxx().val(), tslice)
        return ld
    elif isinstance(prop, SelProp):
        ld.xx = cc.contract_pion((<SelProp>prop).xxx().val(), tslice,
                                 (<SelProp>prop).fsel.xx)
        return ld
    else:
        raise Exception("contract_pion_field")
