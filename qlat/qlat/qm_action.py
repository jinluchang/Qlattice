from qlat_utils import *
from .c import *
from . import c

class QMAction:

    def __init__(self, lmbd, v0, alpha, m_particle, dt):
        self.cdata = c.mk_qm_action(lmbd, v0, alpha, m_particle, dt)

    def __del__(self):
        assert isinstance(self.cdata, int)
        c.free_qm_action(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, QMAction)
        c.set_qm_action(self, v1)
        return self

    def lmbd(self):
        return c.get_lmbd_qm_action(self)

    def v0(self):
        return c.get_v0_qm_action(self)

    def alpha(self):
        return c.get_alpha_qm_action(self)

    def m_particle(self):
        return c.get_m_particle_qm_action(self)

    def dt(self):
        return c.get_dt_qm_action(self)

    def action_node(self, f):
        assert isinstance(f, FieldBase)
        return c.action_node_qm_action(self, f)

    def hmc_m_hamilton_node(self, m):
        assert isinstance(m, FieldBase)
        return c.hmc_m_hamilton_node_qm_action(self, m)
    
    def sum_sq(self, f):
        assert isinstance(f, FieldBase)
        return c.sum_sq_qm_action(self, f)

    def hmc_set_force(self, force, f):
        assert isinstance(force, FieldBase)
        assert isinstance(f, FieldBase)
        return c.hmc_set_force_qm_action(self, force, f)

    def hmc_field_evolve(self, f, m, step_size):
        assert isinstance(f, FieldBase)
        assert isinstance(m, FieldBase)
        c.hmc_field_evolve_qm_action(self, f, m, step_size)

    def hmc_set_rand_momentum(self, m, rs):
        assert isinstance(m, FieldBase)
        assert isinstance(rs, RngState)
        return c.hmc_set_rand_momentum_qm_action(self, m, rs)
