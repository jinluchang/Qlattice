# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

"""
Module ``qlat.qm_action``
============================\n
Quantum-mechanical action for Hamiltonian Monte Carlo (HMC) simulations,
wrapping the C-level ``QMAction`` that defines a confining potential with
configurable barrier strength and FV (finite-volume) parameters.\n
Documentation: ``docs/qlat/qlat_qm_action.md``\n
.. note:: Update the documentation when updating this source file.
"""

from qlat_utils.all cimport *
from .field_base cimport FieldBase
from .field_types cimport FieldRealD
from . cimport everything as cc

from cpython.long cimport PyLong_FromVoidPtr
from cpython.long cimport PyLong_AsVoidPtr

cdef inline cc.QMAction* get_qm_action_ptr(object qma) except? NULL:
    return <cc.QMAction*>PyLong_AsVoidPtr(qma.cdata)

def mk_qm_action(
    alpha,
    beta=0.0,
    V_FV_min=0.0,
    FV_offset=0.0,
    TV_offset=0.0,
    barrier_strength=1.0,
    L=1.0,
    M=0.0,
    epsilon=0.0,
    t_FV_out=10,
    t_FV_mid=5,
    dt=1.0,
    measure_offset_L=False,
    measure_offset_M=False,
):
    cdef cc.QMAction* pqma = new cc.QMAction(
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
    return PyLong_FromVoidPtr(<void*>pqma)

def free_qm_action(qma):
    cdef cc.QMAction* pqma = get_qm_action_ptr(qma)
    del pqma

def set_qm_action(qma_new, qma):
    cdef cc.QMAction* p_qma_new = get_qm_action_ptr(qma_new)
    cdef cc.QMAction* p_qma = get_qm_action_ptr(qma)
    p_qma_new[0] = p_qma[0]

def V_qm_action(qma, x0=0.0, x1=0.0, t=0):
    cdef cc.vector[cc.RealD] x_v = cc.vector[cc.RealD]()
    x_v.resize(2)
    x_v[0] = x0
    x_v[1] = x1
    return get_qm_action_ptr(qma).V(x_v.v, t)

def dV_qm_action(qma, x0=0.0, x1=0.0, t=0, idx=0):
    cdef cc.vector[cc.RealD] x_v = cc.vector[cc.RealD]()
    x_v.resize(2)
    x_v[0] = x0
    x_v[1] = x1
    return get_qm_action_ptr(qma).dV(x_v.v, t, idx)

### -------------------------------------------------------------------
### cqlat-compatible entry points

def get_alpha_qm_action(qma):
    return get_qm_action_ptr(qma).alpha

def get_beta_qm_action(qma):
    return get_qm_action_ptr(qma).beta

def get_barrier_strength_qm_action(qma):
    return get_qm_action_ptr(qma).barrier_strength

def get_M_qm_action(qma):
    return get_qm_action_ptr(qma).M

def get_L_qm_action(qma):
    return get_qm_action_ptr(qma).L

def get_t_FV_out_qm_action(qma):
    return get_qm_action_ptr(qma).t_FV_out

def get_t_FV_mid_qm_action(qma):
    return get_qm_action_ptr(qma).t_FV_mid

def get_dt_qm_action(qma):
    return get_qm_action_ptr(qma).dt

def action_node_qm_action(qma, f):
    return get_qm_action_ptr(qma).action_node((<FieldRealD>f).xx)

def hmc_m_hamilton_node_qm_action(qma, m):
    return get_qm_action_ptr(qma).hmc_m_hamilton_node((<FieldRealD>m).xx)

def sum_sq_qm_action(qma, f):
    return get_qm_action_ptr(qma).sum_sq((<FieldRealD>f).xx)

def hmc_set_force_qm_action(qma, force, f):
    get_qm_action_ptr(qma).hmc_set_force(
        (<FieldRealD>force).xx, (<FieldRealD>f).xx)

def hmc_field_evolve_qm_action(qma, f, m, step_size):
    get_qm_action_ptr(qma).hmc_field_evolve(
        (<FieldRealD>f).xx, (<FieldRealD>m).xx, step_size)

def hmc_set_rand_momentum_qm_action(qma, m, RngState rs):
    get_qm_action_ptr(qma).hmc_set_rand_momentum((<FieldRealD>m).xx, rs.xx)

class QMAction:
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
        self.cdata = mk_qm_action(
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

    def __del__(self):
        assert isinstance(self.cdata, int)
        free_qm_action(self)

    def __imatmul__(self, v1):
        assert isinstance(v1, QMAction)
        set_qm_action(self, v1)
        return self

    def alpha(self):
        return get_qm_action_ptr(self).alpha

    def beta(self):
        return get_qm_action_ptr(self).beta

    def barrier_strength(self):
        return get_qm_action_ptr(self).barrier_strength

    def M(self):
        return get_qm_action_ptr(self).M

    def L(self):
        return get_qm_action_ptr(self).L

    def t_FV_out(self):
        return get_qm_action_ptr(self).t_FV_out

    def t_FV_mid(self):
        return get_qm_action_ptr(self).t_FV_mid

    def t_FV(self):
        return 2 * self.t_FV_out() + self.t_FV_mid()

    def dt(self):
        return get_qm_action_ptr(self).dt

    def V(self, x, t):
        return V_qm_action(self, x[0], x[1], t)

    def dV(self, x, t):
        return dV_qm_action(self, x[0], x[1], t, 0)

    def action_node(self, f):
        assert isinstance(f, FieldBase)
        return get_qm_action_ptr(self).action_node((<FieldRealD>f).xx)

    def hmc_m_hamilton_node(self, m):
        assert isinstance(m, FieldBase)
        return get_qm_action_ptr(self).hmc_m_hamilton_node((<FieldRealD>m).xx)

    def sum_sq(self, f):
        assert isinstance(f, FieldBase)
        return get_qm_action_ptr(self).sum_sq((<FieldRealD>f).xx)

    def hmc_set_force(self, force, f):
        assert isinstance(force, FieldBase)
        assert isinstance(f, FieldBase)
        return get_qm_action_ptr(self).hmc_set_force(
            (<FieldRealD>force).xx, (<FieldRealD>f).xx)

    def hmc_field_evolve(self, f, m, step_size):
        assert isinstance(f, FieldBase)
        assert isinstance(m, FieldBase)
        get_qm_action_ptr(self).hmc_field_evolve(
            (<FieldRealD>f).xx, (<FieldRealD>m).xx, step_size)

    def hmc_set_rand_momentum(self, m, RngState rs):
        assert isinstance(m, FieldBase)
        return get_qm_action_ptr(self).hmc_set_rand_momentum(
            (<FieldRealD>m).xx, rs.xx)
