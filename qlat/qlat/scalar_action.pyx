# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.scalar_action``
==============================\n
Scalar field action for lattice HMC simulations, providing force evaluation,
field evolution, and Hamiltonian routines used by the molecular-dynamics
integrator.\n
Documentation: ``docs/qlat/qlat_scalar_action.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils.all cimport *
from . cimport everything as cc
from .field_base cimport FieldBase
from .field_types cimport FieldRealD
from .field_types cimport FieldComplexD

from cpython.long cimport PyLong_FromVoidPtr
from cpython.long cimport PyLong_AsVoidPtr

import cqlat as c

cdef inline cc.ScalarAction* get_scalar_action_ptr(object sa) except? NULL:
    return <cc.ScalarAction*>PyLong_AsVoidPtr(sa.cdata)

def mk_scalar_action(cc.RealD m_sq, cc.RealD lmbd, cc.RealD alpha):
    cdef cc.ScalarAction* psa = new cc.ScalarAction(m_sq, lmbd, alpha)
    return PyLong_FromVoidPtr(<void*>psa)

class ScalarAction:
    def __init__(self, m_sq, lmbd, alpha):
        self.cdata = mk_scalar_action(m_sq, lmbd, alpha)

    def __del__(self):
        assert isinstance(self.cdata, int)
        c.free_scalar_action(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, ScalarAction)
        c.set_scalar_action(self, v1)
        return self

    def m_sq(self):
        return c.get_m_sq_scalar_action(self)

    def lmbd(self):
        return c.get_lmbd_scalar_action(self)

    def alpha(self):
        return c.get_alpha_scalar_action(self)

    def action_node(self, FieldRealD sf):
        cdef cc.ScalarAction* psa = get_scalar_action_ptr(self)
        return psa.action_node(sf.xx)

    def hmc_estimate_mass(self, masses, field_ft, force_ft, phi0):
        assert isinstance(masses, FieldBase)
        assert isinstance(field_ft, FieldBase)
        assert isinstance(force_ft, FieldBase)
        return c.hmc_estimate_mass_scalar_action(
            self, masses, field_ft, force_ft, phi0)

    def to_mass_factor(self, sin_domega):
        assert isinstance(sin_domega, FieldBase)
        return c.to_mass_factor_scalar_action(self, sin_domega)

    def set_complex_from_double(self, cf, sf):
        assert isinstance(cf, FieldBase)
        assert isinstance(sf, FieldBase)
        return c.set_complex_from_double_scalar_action(self, cf, sf)

    def set_double_from_complex(self, sf, cf):
        assert isinstance(cf, FieldBase)
        assert isinstance(sf, FieldBase)
        return c.set_double_from_complex_scalar_action(self, sf, cf)

    def sum_sq(self, FieldRealD sf):
        cdef cc.ScalarAction* psa = get_scalar_action_ptr(self)
        return psa.sum_sq(sf.xx)

    def hmc_m_hamilton_node(self, FieldComplexD sf, FieldRealD masses):
        cdef cc.ScalarAction* psa = get_scalar_action_ptr(self)
        return psa.hmc_m_hamilton_node(sf.xx, masses.xx)

    def hmc_set_force(self, FieldRealD sm_force, FieldRealD sf):
        cdef cc.ScalarAction* psa = get_scalar_action_ptr(self)
        return psa.hmc_set_force(sm_force.xx, sf.xx)

    def hmc_field_evolve(self, FieldComplexD sf_ft, FieldComplexD sm_ft,
                         FieldRealD masses, cc.RealD step_size):
        cdef cc.ScalarAction* psa = get_scalar_action_ptr(self)
        cdef cc.RealD step_size_ = step_size
        psa.hmc_field_evolve(sf_ft.xx, sm_ft.xx, masses.xx, step_size_)

    def axial_current_node(self, FieldRealD axial_current, FieldRealD sf):
        cdef cc.ScalarAction* psa = get_scalar_action_ptr(self)
        return psa.axial_current_node(axial_current.xx, sf.xx)

    def hmc_set_rand_momentum(self, FieldComplexD sm_complex,
                              FieldRealD masses, RngState rs):
        cdef cc.ScalarAction* psa = get_scalar_action_ptr(self)
        return psa.hmc_set_rand_momentum(sm_complex.xx, masses.xx, rs.xx)

    def hmc_predict_field(self, FieldComplexD field_ft,
                          FieldComplexD momentum_ft, FieldRealD masses,
                          cc.RealD vev_sigma):
        cdef cc.ScalarAction* psa = get_scalar_action_ptr(self)
        cdef cc.RealD vev_sigma_ = vev_sigma
        return psa.hmc_predict_field(
            field_ft.xx, momentum_ft.xx, masses.xx, vev_sigma_)

    def get_polar_field(self, FieldRealD polar_field, FieldRealD field):
        cdef cc.ScalarAction* psa = get_scalar_action_ptr(self)
        return psa.get_polar_field(polar_field.xx, field.xx)
