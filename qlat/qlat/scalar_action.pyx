# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.scalar_action``
==============================\n
Scalar field action for lattice HMC simulations, providing force evaluation,
field evolution, and Hamiltonian routines used by the molecular-dynamics
integrator.  ``ScalarAction`` is a ``cdef class`` that owns its C++ object by
value, like ``GaugeAction``.  The pure-Python ``ScalarAction`` entry points
``hmc_estimate_mass_scalar_action`` and ``to_mass_factor_scalar_action`` live
in ``qlat.scalar_action_utils``.\n
Documentation: ``docs/qlat/qlat_scalar_action.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils.all cimport *
from . cimport everything as cc
from .field_base cimport FieldBase
from .field_types cimport FieldRealD
from .field_types cimport FieldComplexD

from . import field_double as field_double

cdef class ScalarAction:

    cdef cc.ScalarAction xx

    def __cinit__(self):
        self.xx = cc.ScalarAction()

    def __init__(self, cc.RealD m_sq, cc.RealD lmbd, cc.RealD alpha):
        self.xx = cc.ScalarAction(m_sq, lmbd, alpha)

    def __imatmul__(self, ScalarAction v1):
        cc.assign_direct(self.xx, v1.xx)
        return self

    def m_sq(self):
        return self.xx.m_sq

    def lmbd(self):
        return self.xx.lmbd

    def alpha(self):
        return self.xx.alpha

    def action_node(self, FieldRealD sf):
        return self.xx.action_node(sf.xx)

    def hmc_estimate_mass(self, masses, field_ft, force_ft, phi0):
        assert isinstance(masses, FieldBase)
        assert isinstance(field_ft, FieldBase)
        assert isinstance(force_ft, FieldBase)
        self.xx.hmc_estimate_mass((<FieldRealD>masses).xx,
                                  (<FieldComplexD>field_ft).xx,
                                  (<FieldComplexD>force_ft).xx, phi0)

    def to_mass_factor(self, sin_domega):
        assert isinstance(sin_domega, FieldBase)
        self.xx.to_mass_factor((<FieldRealD>sin_domega).xx)

    def set_complex_from_double(self, cf, sf):
        assert isinstance(cf, FieldBase)
        assert isinstance(sf, FieldBase)
        return field_double.set_complex_from_double(cf, sf)

    def set_double_from_complex(self, sf, cf):
        assert isinstance(cf, FieldBase)
        assert isinstance(sf, FieldBase)
        return field_double.set_double_from_complex(sf, cf)

    def sum_sq(self, FieldRealD sf):
        return self.xx.sum_sq(sf.xx)

    def hmc_m_hamilton_node(self, FieldComplexD sf, FieldRealD masses):
        return self.xx.hmc_m_hamilton_node(sf.xx, masses.xx)

    def hmc_set_force(self, FieldRealD sm_force, FieldRealD sf):
        return self.xx.hmc_set_force(sm_force.xx, sf.xx)

    def hmc_field_evolve(self, FieldComplexD sf_ft, FieldComplexD sm_ft,
                         FieldRealD masses, cc.RealD step_size):
        cdef cc.RealD step_size_ = step_size
        self.xx.hmc_field_evolve(sf_ft.xx, sm_ft.xx, masses.xx, step_size_)

    def axial_current_node(self, FieldRealD axial_current, FieldRealD sf):
        return self.xx.axial_current_node(axial_current.xx, sf.xx)

    def hmc_set_rand_momentum(self, FieldComplexD sm_complex,
                              FieldRealD masses, RngState rs):
        return self.xx.hmc_set_rand_momentum(sm_complex.xx, masses.xx, rs.xx)

    def hmc_predict_field(self, FieldComplexD field_ft,
                          FieldComplexD momentum_ft, FieldRealD masses,
                          cc.RealD vev_sigma):
        cdef cc.RealD vev_sigma_ = vev_sigma
        return self.xx.hmc_predict_field(
            field_ft.xx, momentum_ft.xx, masses.xx, vev_sigma_)

    def get_polar_field(self, FieldRealD polar_field, FieldRealD field):
        return self.xx.get_polar_field(polar_field.xx, field.xx)

### -------------------------------------------------------------------
### cqlat-compatible entry points

def set_scalar_action(ScalarAction sa_new, ScalarAction sa):
    sa_new.xx = sa.xx

def get_m_sq_scalar_action(ScalarAction sa):
    return sa.xx.m_sq

def get_lmbd_scalar_action(ScalarAction sa):
    return sa.xx.lmbd

def get_alpha_scalar_action(ScalarAction sa):
    return sa.xx.alpha
