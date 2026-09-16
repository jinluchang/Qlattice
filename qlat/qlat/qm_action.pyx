# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.qm_action``
============================\n
Quantum-mechanical action for Hamiltonian Monte Carlo (HMC) simulations.
Defines a confining potential with configurable barrier strength and FV
(finite-volume) parameters.  ``QMAction`` is a ``cdef class`` that owns its
C++ object by value.\n
Documentation: ``docs/qlat/qlat_qm_action.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils.all cimport *
from .field_base cimport FieldBase
from .field_types cimport FieldRealD
from . cimport everything as cc

cdef class QMAction:

    cdef cc.QMAction xx

    def __cinit__(self):
        self.xx = cc.QMAction()

    def __init__(
        self,
        alpha,
        beta,
        V_FV_min,
        FV_offset,
        TV_offset,
        barrier_strength,
        L,
        M,
        epsilon,
        t_FV_out,
        t_FV_mid,
        dt,
        measure_offset_L,
        measure_offset_M,
    ):
        self.xx = cc.QMAction(
            alpha,
            beta,
            V_FV_min,
            FV_offset,
            TV_offset,
            barrier_strength,
            L,
            M,
            epsilon,
            t_FV_out,
            t_FV_mid,
            dt,
            measure_offset_L,
            measure_offset_M,
        )

    def __imatmul__(self, QMAction v1):
        self.xx = v1.xx
        return self

    def alpha(self):
        return self.xx.alpha

    def beta(self):
        return self.xx.beta

    def barrier_strength(self):
        return self.xx.barrier_strength

    def M(self):
        return self.xx.M

    def L(self):
        return self.xx.L

    def t_FV_out(self):
        return self.xx.t_FV_out

    def t_FV_mid(self):
        return self.xx.t_FV_mid

    def t_FV(self):
        return 2 * self.t_FV_out() + self.t_FV_mid()

    def dt(self):
        return self.xx.dt

    def V(self, x, t):
        return V_qm_action(self, x[0], x[1], t)

    def dV(self, x, t):
        return dV_qm_action(self, x[0], x[1], t, 0)

    def action_node(self, f):
        assert isinstance(f, FieldBase)
        return self.xx.action_node((<FieldRealD>f).xx)

    def hmc_m_hamilton_node(self, m):
        assert isinstance(m, FieldBase)
        return self.xx.hmc_m_hamilton_node((<FieldRealD>m).xx)

    def sum_sq(self, f):
        assert isinstance(f, FieldBase)
        return self.xx.sum_sq((<FieldRealD>f).xx)

    def hmc_set_force(self, force, f):
        assert isinstance(force, FieldBase)
        assert isinstance(f, FieldBase)
        self.xx.hmc_set_force((<FieldRealD>force).xx, (<FieldRealD>f).xx)

    def hmc_field_evolve(self, f, m, step_size):
        assert isinstance(f, FieldBase)
        assert isinstance(m, FieldBase)
        self.xx.hmc_field_evolve(
            (<FieldRealD>f).xx, (<FieldRealD>m).xx, step_size)

    def hmc_set_rand_momentum(self, m, RngState rs):
        assert isinstance(m, FieldBase)
        self.xx.hmc_set_rand_momentum((<FieldRealD>m).xx, rs.xx)

### -------------------------------------------------------------------
### cqlat-compatible entry points

def set_qm_action(QMAction qma_new, QMAction qma):
    qma_new.xx = qma.xx

def V_qm_action(QMAction qma, x0=0.0, x1=0.0, t=0):
    cdef cc.vector[cc.RealD] x_v = cc.vector[cc.RealD]()
    x_v.resize(2)
    x_v[0] = x0
    x_v[1] = x1
    return qma.xx.V(x_v.v, t)

def dV_qm_action(QMAction qma, x0=0.0, x1=0.0, t=0, idx=0):
    cdef cc.vector[cc.RealD] x_v = cc.vector[cc.RealD]()
    x_v.resize(2)
    x_v[0] = x0
    x_v[1] = x1
    return qma.xx.dV(x_v.v, t, idx)

def get_alpha_qm_action(QMAction qma):
    return qma.alpha()

def get_beta_qm_action(QMAction qma):
    return qma.beta()

def get_barrier_strength_qm_action(QMAction qma):
    return qma.barrier_strength()

def get_M_qm_action(QMAction qma):
    return qma.M()

def get_L_qm_action(QMAction qma):
    return qma.L()

def get_t_FV_out_qm_action(QMAction qma):
    return qma.t_FV_out()

def get_t_FV_mid_qm_action(QMAction qma):
    return qma.t_FV_mid()

def get_dt_qm_action(QMAction qma):
    return qma.dt()

def action_node_qm_action(QMAction qma, f):
    return qma.action_node(f)

def hmc_m_hamilton_node_qm_action(QMAction qma, m):
    return qma.hmc_m_hamilton_node(m)

def sum_sq_qm_action(QMAction qma, f):
    return qma.sum_sq(f)

def hmc_set_force_qm_action(QMAction qma, force, f):
    qma.hmc_set_force(force, f)

def hmc_field_evolve_qm_action(QMAction qma, f, m, step_size):
    qma.hmc_field_evolve(f, m, step_size)

def hmc_set_rand_momentum_qm_action(QMAction qma, m, RngState rs):
    qma.hmc_set_rand_momentum(m, rs)
