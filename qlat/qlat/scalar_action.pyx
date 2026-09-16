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

from . import field_double as field_double

cdef inline cc.ScalarAction* get_scalar_action_ptr(object sa) except? NULL:
    return <cc.ScalarAction*>PyLong_AsVoidPtr(sa.cdata)

def mk_scalar_action(cc.RealD m_sq, cc.RealD lmbd, cc.RealD alpha):
    cdef cc.ScalarAction* psa = new cc.ScalarAction(m_sq, lmbd, alpha)
    return PyLong_FromVoidPtr(<void*>psa)

def free_scalar_action(sa):
    cdef cc.ScalarAction* psa = get_scalar_action_ptr(sa)
    del psa

def set_scalar_action(sa_new, sa):
    cdef cc.ScalarAction* p_sa_new = get_scalar_action_ptr(sa_new)
    cdef cc.ScalarAction* p_sa = get_scalar_action_ptr(sa)
    p_sa_new[0] = p_sa[0]

def get_m_sq_scalar_action(sa):
    return get_scalar_action_ptr(sa).m_sq

def get_lmbd_scalar_action(sa):
    return get_scalar_action_ptr(sa).lmbd

def get_alpha_scalar_action(sa):
    return get_scalar_action_ptr(sa).alpha

class ScalarAction:
    def __init__(self, m_sq, lmbd, alpha):
        self.cdata = mk_scalar_action(m_sq, lmbd, alpha)

    def __del__(self):
        assert isinstance(self.cdata, int)
        free_scalar_action(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, ScalarAction)
        set_scalar_action(self, v1)
        return self

    def m_sq(self):
        return get_scalar_action_ptr(self).m_sq

    def lmbd(self):
        return get_scalar_action_ptr(self).lmbd

    def alpha(self):
        return get_scalar_action_ptr(self).alpha

    def action_node(self, FieldRealD sf):
        cdef cc.ScalarAction* psa = get_scalar_action_ptr(self)
        return psa.action_node(sf.xx)

    def hmc_estimate_mass(self, masses, field_ft, force_ft, phi0):
        assert isinstance(masses, FieldBase)
        assert isinstance(field_ft, FieldBase)
        assert isinstance(force_ft, FieldBase)
        cdef cc.ScalarAction* psa = get_scalar_action_ptr(self)
        psa.hmc_estimate_mass((<FieldRealD>masses).xx,
                              (<FieldComplexD>field_ft).xx,
                              (<FieldComplexD>force_ft).xx, phi0)

    def to_mass_factor(self, sin_domega):
        assert isinstance(sin_domega, FieldBase)
        cdef cc.ScalarAction* psa = get_scalar_action_ptr(self)
        psa.to_mass_factor((<FieldRealD>sin_domega).xx)

    def set_complex_from_double(self, cf, sf):
        assert isinstance(cf, FieldBase)
        assert isinstance(sf, FieldBase)
        return field_double.set_complex_from_double(cf, sf)

    def set_double_from_complex(self, sf, cf):
        assert isinstance(cf, FieldBase)
        assert isinstance(sf, FieldBase)
        return field_double.set_double_from_complex(sf, cf)

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
