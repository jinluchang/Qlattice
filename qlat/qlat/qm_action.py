from qlat_utils import *
from .c import *
from . import c

class QMAction:

    def __init__(self, alpha, beta, FV_offset, TV_offset, barrier_strength, L, M, epsilon, t_FV_out, t_FV_mid, dt, measure_offset_L, measure_offset_M):
        self.cdata = c.mk_qm_action(alpha, beta, FV_offset, TV_offset, barrier_strength, L, M, epsilon, t_FV_out, t_FV_mid, dt, measure_offset_L, measure_offset_M)

    def __del__(self):
        assert isinstance(self.cdata, int)
        c.free_qm_action(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, QMAction)
        c.set_qm_action(self, v1)
        return self

    def alpha(self):
        return c.get_alpha_qm_action(self)

    def beta(self):
        return c.get_beta_qm_action(self)

    def barrier_strength(self):
        return c.get_barrier_strength_qm_action(self)

    def M(self):
        return c.get_M_qm_action(self)

    def L(self):
        return c.get_L_qm_action(self)

    def t_FV_out(self):
        return c.get_t_FV_out_qm_action(self)

    def t_FV_mid(self):
        return c.get_t_FV_mid_qm_action(self)
    
    def t_FV(self):
        return 2*self.t_FV_out() + self.t_FV_mid()

    def dt(self):
        return c.get_dt_qm_action(self)

    def V(self, x, t):
        return c.V_qm_action(self, x, t)

    def dV(self, x, t):
        return c.dV_qm_action(self, x, t)

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
