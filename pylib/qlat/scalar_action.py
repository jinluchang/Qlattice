import cqlat as c

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
    
    def action_node(self, sm):
        assert isinstance(sm, Field)
        return c.action_node_scalar_action(self, sm)
    
    def set_complex_from_double(self, cf, sf):
        assert isinstance(cf, Field)
        assert isinstance(sf, Field)
        return c.set_complex_from_double_scalar_action(self, cf, sf)
    
    def set_double_from_complex(self, sf, cf):
        assert isinstance(cf, Field)
        assert isinstance(sf, Field)
        return c.set_double_from_complex_scalar_action(self, sf, cf)
    
    def sum_sq(self, sf):
        assert isinstance(sf, Field)
        return c.sum_sq_scalar_action(self, sf)
    
    def hmc_m_hamilton_node(self, sf):
        assert isinstance(sf, Field)
        return c.hmc_m_hamilton_node_scalar_action(self, sf)
    
    def hmc_set_force(self, sm_force, sf):
        assert isinstance(sm_force, Field)
        assert isinstance(sf, Field)
        return c.hmc_set_force_scalar_action(self, sm_force, sf)
    
    def hmc_field_evolve(self, sf, sm, sf_complex, sm_complex, step_size, V):
        assert isinstance(sf, Field)
        assert isinstance(sm, Field)
        self.set_complex_from_double(sf_complex, sf)
        self.set_complex_from_double(sm_complex, sm)
        fft = mk_fft(True)
        sf_complex=fft*sf_complex
        sf_complex*=1/V**0.5
        sm_complex=fft*sm_complex
        sm_complex*=1/V**0.5
        c.hmc_field_evolve_scalar_action(self, sf_complex, sm_complex, step_size)
        ifft = mk_fft(False)
        sf_complex=ifft*sf_complex
        sf_complex*=1/V**0.5
        self.set_double_from_complex(sf, sf_complex)
    
    def axial_current_node(self, axial_current, sf):
        assert isinstance(sf, Field)
        assert isinstance(axial_current, Field)
        return c.axial_current_node_scalar_action(self, axial_current, sf)
