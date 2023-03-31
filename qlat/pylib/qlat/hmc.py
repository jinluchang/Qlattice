from qlat_utils import *
import qlat.c as c

from qlat.field import *
from qlat.qcd import *
from qlat.gauge_action import *
from qlat.elem_type import *

import math

class GaugeMomentum(FieldColorMatrix):

    def __init__(self, geo = None):
        super().__init__(geo, 4)

    def set_rand(self, rng, sigma = 1.0):
        set_rand_gauge_momentum(self, sigma, rng)

def set_rand_gauge_momentum(gm, sigma, rng):
    assert isinstance(gm, GaugeMomentum)
    assert isinstance(sigma, float)
    assert isinstance(rng, RngState)
    return c.set_rand_gauge_momentum(gm, sigma, rng)

def gm_hamilton_node(gm):
    assert isinstance(gm, GaugeMomentum)
    return c.gm_hamilton_node(gm)

def gf_hamilton_node(gf, ga):
    assert isinstance(gf, GaugeField)
    assert isinstance(ga, GaugeAction)
    return c.gf_hamilton_node(gf, ga)

def gf_evolve(gf, gm, step_size):
    assert isinstance(gm, GaugeMomentum)
    return c.gf_evolve(gf, gm, step_size)

def set_gm_force(gm_force, gf, ga):
    assert isinstance(gm_force, GaugeMomentum)
    assert isinstance(gf, GaugeField)
    assert isinstance(ga, GaugeAction)
    return c.set_gm_force(gm_force, gf, ga)

@timer_verbose
def metropolis_accept(delta_h, traj, rs):
    flag_d = 0.0
    accept_prob = 0.0
    if get_id_node() == 0:
        if delta_h <= 0.0:
            accept_prob = 1.0
            flag_d = 1.0
        else:
            accept_prob = math.exp(-delta_h)
            rand_num = rs.u_rand_gen(1.0, 0.0)
            if rand_num <= accept_prob:
                flag_d = 1.0
    flag_d = glb_sum(flag_d)
    accept_prob = glb_sum(accept_prob)
    flag = flag_d > 0.5
    displayln_info(f"metropolis_accept: flag={flag:d} with accept_prob={accept_prob * 100.0:.1f}% delta_h={delta_h:.16f} traj={traj}")
    return flag, accept_prob

@timer
def gm_evolve_fg_pure_gauge(gm, gf_init, ga, fg_dt, dt):
    geo = gf_init.geo()
    gf = GaugeField(geo)
    gf @= gf_init
    gm_force = GaugeMomentum(geo)
    set_gm_force(gm_force, gf, ga)
    gf_evolve(gf, gm_force, fg_dt)
    set_gm_force(gm_force, gf, ga)
    gm_force *= dt
    gm += gm_force

@timer_verbose
def run_hmc_evolve_pure_gauge(gm, gf, ga, rs, n_steps, md_time = 1.0):
    energy = gm_hamilton_node(gm) + gf_hamilton_node(gf, ga)
    dt = md_time / n_steps
    lam = 0.5 * (1.0 - 1.0 / math.sqrt(3.0));
    theta = (2.0 - math.sqrt(3.0)) / 48.0;
    ttheta = theta * dt * dt * dt;
    gf_evolve(gf, gm, lam * dt)
    for i in range(n_steps):
        gm_evolve_fg_pure_gauge(gm, gf, ga, 4.0 * ttheta / dt, 0.5 * dt);
        gf_evolve(gf, gm, (1.0 - 2.0 * lam) * dt);
        gm_evolve_fg_pure_gauge(gm, gf, ga, 4.0 * ttheta / dt, 0.5 * dt);
        if i < n_steps - 1:
            gf_evolve(gf, gm, 2.0 * lam * dt);
        else:
            gf_evolve(gf, gm, lam * dt);
    gf.unitarize()
    delta_h = gm_hamilton_node(gm) + gf_hamilton_node(gf, ga) - energy;
    delta_h = glb_sum(delta_h)
    return delta_h

@timer_verbose
def run_hmc_pure_gauge(gf, ga, traj, rs, *, is_reverse_test = False, n_steps = 6, md_time = 1.0, is_always_accept = False):
    fname = get_fname()
    rs = rs.split(f"{traj}")
    geo = gf.geo()
    gf0 = GaugeField(geo)
    gf0 @= gf
    gm = GaugeMomentum(geo)
    gm.set_rand(rs.split("set_rand_gauge_momentum"), 1.0)
    delta_h = run_hmc_evolve_pure_gauge(gm, gf0, ga, rs, n_steps, md_time)
    if is_reverse_test:
        gm_r = GaugeMomentum(geo)
        gm_r @= gm
        gf0_r = GaugeField(geo)
        gf0_r @= gf0
        delta_h_rev = run_hmc_evolve_pure_gauge(gm_r, gf0_r, ga, rs, n_steps, -md_time)
        gf0_r -= gf;
        displayln_info(f"{fname}: reversed delta_diff: {delta_h + delta_h_rev} / {delta_h}")
        displayln_info(f"{fname}: reversed gf_diff: {qnorm(gf0_r)} / {qnorm(gf0)}")
    flag, accept_prob = metropolis_accept(delta_h, traj, rs.split("metropolis_accept"))
    if flag or is_always_accept:
        displayln_info(f"{fname}: update gf (traj={traj})")
        gf @= gf0
