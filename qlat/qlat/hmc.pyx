# cython: binding=True, embedsignature=True, c_string_type=unicode, c_string_encoding=utf8

from qlat_utils.all cimport *
from . cimport everything as cc
from .geometry cimport Geometry
from .gauge_action cimport GaugeAction
from .field_types cimport (
        FieldColorMatrix,
        FieldRealD,
        )
from .qcd cimport GaugeField

import cqlat as c
import qlat_utils as q
import numpy as np
from .mpi import glb_sum_double

import math

cdef class GaugeMomentum(FieldColorMatrix):

    def __init__(self, Geometry geo=None):
        super().__init__(geo, 4)

    cdef cc.Handle[cc.GaugeMomentum] xxx(self):
        assert self.xx.multiplicity == 4 or self.xx.multiplicity == 0
        return cc.Handle[cc.GaugeMomentum](<cc.GaugeMomentum&>self.xx)

    def set_rand(self, RngState rng, cc.RealD sigma=1.0):
        set_rand_gauge_momentum(self, sigma, rng)

    def set_rand_fa(self, FieldRealD mf, RngState rng):
        set_rand_gauge_momentum_fa(self, mf, rng)

###

def set_rand_gauge_momentum(GaugeMomentum gm, cc.RealD sigma, RngState rng):
    return cc.set_rand_gauge_momentum(gm.xxx().val(), sigma, rng.xx)

def set_rand_gauge_momentum_fa(GaugeMomentum gm, FieldRealD mf, RngState rng):
    return cc.set_rand_gauge_momentum(gm.xxx().val(), mf.xx, rng.xx)

def gm_hamilton_node(GaugeMomentum gm):
    return cc.gm_hamilton_node(gm.xxx().val())

def gm_hamilton_node_fa(GaugeMomentum gm, FieldRealD mf):
    return cc.gm_hamilton_node(gm.xxx().val(), mf.xx)

def gf_hamilton_node(GaugeField gf, GaugeAction ga):
    return cc.gf_hamilton_node(gf.xxx().val(), ga.xx)

def gf_evolve(GaugeField gf, GaugeMomentum gm, cc.RealD step_size):
    return cc.gf_evolve(gf.xxx().val(), gm.xxx().val(), step_size)

def gf_evolve_dual(GaugeField gf, GaugeMomentum gm_dual, cc.RealD step_size):
    return cc.gf_evolve_dual(gf.xxx().val(), gm_dual.xxx().val(), step_size)

def gf_evolve_fa(GaugeField gf, GaugeMomentum gm, FieldRealD mf, cc.RealD step_size):
    return cc.gf_evolve(gf.xxx().val(), gm.xxx().val(), mf.xx, step_size)

def gf_evolve_fa_dual(GaugeField gf, GaugeMomentum gm_dual, FieldRealD mf_dual, cc.RealD step_size):
    return cc.gf_evolve_dual(gf.xxx().val(), gm_dual.xxx().val(), mf_dual.xx, step_size)

def set_gm_force(GaugeMomentum gm_force, GaugeField gf, GaugeAction ga):
    return cc.set_gm_force(gm_force.xxx().val(), gf.xxx().val(), ga.xx)

def set_gm_force_dual(GaugeMomentum gm_force_dual, GaugeField gf, GaugeMomentum gm_force):
    return cc.set_gm_force_dual(gm_force_dual.xxx().val(), gf.xxx().val(), gm_force.xxx().val())

def project_gauge_transform(GaugeMomentum gm, GaugeMomentum gm_dual, FieldRealD mf, FieldRealD mf_dual):
    """
    Project out the pure gauge transformation movement.
    """
    return cc.project_gauge_transform(gm.xxx().val(), gm_dual.xxx().val(), mf.xx, mf_dual.xx)

def dot_gauge_momentum(GaugeMomentum gm1, GaugeMomentum gm2):
    """
    return RealD dot field with multiplicity = 4
    """
    cdef FieldRealD f = FieldRealD()
    cc.dot_gauge_momentum(f.xx, gm1.xxx().val(), gm2.xxx().val())
    return f

@q.timer_verbose
def metropolis_accept(delta_h, traj, rs):
    flag_d = 0.0
    accept_prob = 0.0
    if q.get_id_node() == 0:
        if delta_h <= 0.0:
            accept_prob = 1.0
            flag_d = 1.0
        else:
            accept_prob = math.exp(-delta_h)
            rand_num = rs.u_rand_gen(1.0, 0.0)
            if rand_num <= accept_prob:
                flag_d = 1.0
    flag_d = glb_sum_double(flag_d)
    accept_prob = glb_sum_double(accept_prob)
    flag = flag_d > 0.5
    q.displayln_info(f"metropolis_accept: flag={flag:d} with accept_prob={accept_prob * 100.0:.1f}% delta_h={delta_h:.16f} traj={traj}")
    return flag, accept_prob

@q.timer
def gm_evolve_fg_pure_gauge(gm, gf_init, ga, fg_dt, dt):
    geo = gf_init.geo
    gf = GaugeField(geo)
    gf @= gf_init
    gm_force = GaugeMomentum(geo)
    set_gm_force(gm_force, gf, ga)
    gf_evolve(gf, gm_force, fg_dt)
    set_gm_force(gm_force, gf, ga)
    gm_force *= dt
    gm += gm_force

@q.timer_verbose
def run_hmc_evolve_pure_gauge(gm, gf, ga, rs, n_step, md_time=1.0):
    energy = gm_hamilton_node(gm) + gf_hamilton_node(gf, ga)
    dt = md_time / n_step
    lam = 0.5 * (1.0 - 1.0 / math.sqrt(3.0));
    theta = (2.0 - math.sqrt(3.0)) / 48.0;
    ttheta = theta * dt * dt * dt;
    gf_evolve(gf, gm, lam * dt)
    for i in range(n_step):
        gm_evolve_fg_pure_gauge(gm, gf, ga, 4.0 * ttheta / dt, 0.5 * dt);
        gf_evolve(gf, gm, (1.0 - 2.0 * lam) * dt);
        gm_evolve_fg_pure_gauge(gm, gf, ga, 4.0 * ttheta / dt, 0.5 * dt);
        if i < n_step - 1:
            gf_evolve(gf, gm, 2.0 * lam * dt);
        else:
            gf_evolve(gf, gm, lam * dt);
    gf.unitarize()
    delta_h = gm_hamilton_node(gm) + gf_hamilton_node(gf, ga) - energy;
    delta_h = glb_sum_double(delta_h)
    return delta_h

@q.timer_verbose
def run_hmc_pure_gauge(gf, ga, traj, rs, *, is_reverse_test=False, n_step=6, md_time=1.0, is_always_accept=False):
    fname = q.get_fname()
    rs = rs.split(f"{traj}")
    geo = gf.geo
    gf0 = GaugeField(geo)
    gf0 @= gf
    gm = GaugeMomentum(geo)
    gm.set_rand(rs.split("set_rand_gauge_momentum"), 1.0)
    delta_h = run_hmc_evolve_pure_gauge(gm, gf0, ga, rs, n_step, md_time)
    if is_reverse_test:
        gm_r = GaugeMomentum(geo)
        gm_r @= gm
        gf0_r = GaugeField(geo)
        gf0_r @= gf0
        delta_h_rev = run_hmc_evolve_pure_gauge(gm_r, gf0_r, ga, rs, n_step, -md_time)
        gf0_r -= gf
        q.displayln_info(f"{fname}: reversed delta_diff: {delta_h + delta_h_rev} / {delta_h}")
        gf_diff_norm = q.qnorm(gf0_r)
        gf_norm = q.qnorm(gf0)
        q.displayln_info(f"{fname}: reversed gf_diff: {gf_diff_norm} / {gf_norm}")
        assert gf_diff_norm <= 1e-12 * gf_norm
    flag, accept_prob = metropolis_accept(delta_h, traj, rs.split("metropolis_accept"))
    if flag or is_always_accept:
        q.displayln_info(f"{fname}: update gf (traj={traj})")
        gf @= gf0
    return flag, delta_h
