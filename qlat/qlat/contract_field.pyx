# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.contract_field``
===============================\n
Field-level contraction routines for computing the conserved-point
hadronic vacuum polarization tensor with all 16 spin-color
polarization combinations at every lattice site.\n
Documentation: ``docs/qlat/qlat_contract_field.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils.all cimport *
from . cimport everything as cc
from .field_types cimport FieldComplexD
from .propagator cimport Prop

class q:
    from qlat_utils import (
        timer,
    )

@q.timer
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
    cdef FieldComplexD chvp_16 = FieldComplexD()
    cc.py_contract_chvp_16(chvp_16.xx, (<Prop>prop1).xxx().val(),
                           (<Prop>prop2).xxx().val())
    return chvp_16
