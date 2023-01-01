import qlat.c as c

from qlat.field import *
from qlat.field_utils import *

class ScalarAction:

    def __init__(self, m_sq, lmbd, alpha):
        self.cdata = c.mk_scalar_action(m_sq, lmbd, alpha)

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

    def lmbd(self):
        return c.get_alpha_scalar_action(self)

    def action_node(self, sf):
        assert isinstance(sf, FieldBase)
        return c.action_node_scalar_action(self, sf)

    def hmc_estimate_mass(self, masses, field_ft, force_ft, phi0):
        assert isinstance(masses, FieldBase)
        assert isinstance(field_ft, FieldBase)
        assert isinstance(force_ft, FieldBase)
        return c.hmc_estimate_mass_scalar_action(self, masses, field_ft, force_ft, phi0)

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

    def sum_sq(self, sf):
        assert isinstance(sf, FieldBase)
        return c.sum_sq_scalar_action(self, sf)

    def hmc_m_hamilton_node(self, sf, masses):
        assert isinstance(sf, FieldBase)
        assert isinstance(masses, FieldBase)
        return c.hmc_m_hamilton_node_scalar_action(self, sf, masses)

    def hmc_set_force(self, sm_force, sf):
        assert isinstance(sm_force, FieldBase)
        assert isinstance(sf, FieldBase)
        return c.hmc_set_force_scalar_action(self, sm_force, sf)

    def hmc_field_evolve(self, sf_ft, sm_ft, masses, step_size):
        assert isinstance(sf_ft, FieldBase)
        assert isinstance(sm_ft, FieldBase)
        assert isinstance(masses, FieldBase)
        c.hmc_field_evolve_scalar_action(self, sf_ft, sm_ft, masses, step_size)

    def axial_current_node(self, axial_current, sf):
        assert isinstance(sf, FieldBase)
        assert isinstance(axial_current, FieldBase)
        return c.axial_current_node_scalar_action(self, axial_current, sf)

    def hmc_set_rand_momentum(self, sm_complex, masses, rs):
        assert isinstance(sm_complex, FieldBase)
        assert isinstance(masses, FieldBase)
        assert isinstance(rs, RngState)
        return c.hmc_set_rand_momentum_scalar_action(self, sm_complex, masses, rs)

    def hmc_predict_field(self, field_ft, momentum_ft, masses, vev_sigma):
        assert isinstance(field_ft, FieldBase)
        assert isinstance(momentum_ft, FieldBase)
        assert isinstance(masses, FieldBase)
        return c.hmc_predict_field_scalar_action(self, field_ft, momentum_ft, masses, vev_sigma)

    def get_polar_field(self, polar_field, field):
        assert isinstance(polar_field, FieldBase)
        assert isinstance(field, FieldBase)
        return c.get_polar_field_scalar_action(self, polar_field, field)
