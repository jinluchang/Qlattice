"""
Utility functions for HMC that do not call C++ functions directly.\n
MD Integrator
=============\n
The molecular dynamics step uses the **force-gradient integrator** (a.k.a.
the "minimum norm" or "Omelyan" integrator).  Key parameters::\n
    lam = (1 - 1 / sqrt(3)) / 2
    theta = (2 - sqrt(3)) / 48\n
The algorithm, in terms of :math:`dt = md\_time / n\_step`::\n
    U += lam * dt * P                                       # half position update
    for i in range(n_step):
        P += dt/2 * F( evolve(U, 4*theta*dt^2 * F(U)) )    # force-gradient step
        U += (1 - 2*lam) * dt * P                           # main position step
        P += dt/2 * F( evolve(U, 4*theta*dt^2 * F(U)) )    # force-gradient step
        if i < n_step - 1:
            U += 2*lam * dt * P
        else:
            U += lam * dt * P                               # half position update (last)\n
``evolve(U, fg_dt)`` evolves the gauge field forward by ``fg_dt`` using
*the force* as the momentum (one step along the force direction).  The
force is then re-evaluated at this intermediate point, and the averaged
force is applied to :math:`P`.\n
Error scaling::\n
    single-step error   ~ O(dt^5)
    trajectory error    ~ O(dt^4)    (fixed md_time, n_step = md_time / dt)\n
The force-gradient term adds an O(dt^3) correction to the momentum update,
raising the effective order from the standard leapfrog O(dt^2) trajectory
error to O(dt^4).
"""

import math

import qlat_utils as q

from .mpi import glb_sum_double
from .qcd import GaugeField
from .hmc import (
    GaugeMomentum,
    set_gm_force,
    gf_evolve,
    gm_hamilton_node,
    gf_hamilton_node,
)

@q.timer_verbose
def metropolis_accept(delta_h, traj, rs):
    """Metropolis accept/reject step.\n
    Accept with probability min(1, exp(-delta_h)).  ``delta_h`` is the
    energy violation from the MD trajectory (already MPI-summed via
    :func:`run_hmc_evolve_pure_gauge`).  Returns ``(flag, accept_prob)``.
    """
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
    q.displayln_info(
        f"metropolis_accept: flag={flag:d} with accept_prob={accept_prob * 100.0:.1f}% delta_h={delta_h:.16f} traj={traj}"
    )
    return flag, accept_prob

@q.timer
def gm_evolve_fg_pure_gauge(gm, gf_init, ga, fg_dt, dt):
    """Force-gradient momentum update (one ``gm_update`` in the pseudo-code above).\n
    Evaluate the force at :math:`U`, evolve :math:`U` forward by ``fg_dt``
    along the force direction, re-evaluate the force at the intermediate
    point, and add ``dt * F_intermediate`` to :math:`P`.
    """
    geo = gf_init.geo
    gf = GaugeField(geo)
    gf @= gf_init
    gm_force = GaugeMomentum(geo)
    set_gm_force(gm_force, gf, ga)
    gf_evolve(gf, gm_force, fg_dt)
    set_gm_force(gm_force, gf, ga)
    gm_force *= dt
    gm += gm_force

@q.timer(is_timer_fork=True)
def run_hmc_evolve_pure_gauge(gm, gf, ga, rs, n_step, md_time=1.0):
    """Run the MD evolution (force-gradient integrator).\n
    Evolve ``gf`` and ``gm`` in-place for ``n_step`` steps of size
    ``md_time / n_step``.  Returns the energy violation ``delta_h``
    (already MPI-summed).  ``rs`` is not used (accepted for interface
    compatibility).
    """
    energy = gm_hamilton_node(gm) + gf_hamilton_node(gf, ga)
    dt = md_time / n_step
    lam = 0.5 * (1.0 - 1.0 / math.sqrt(3.0))
    theta = (2.0 - math.sqrt(3.0)) / 48.0
    ttheta = theta * dt * dt * dt
    gf_evolve(gf, gm, lam * dt)
    for i in range(n_step):
        gm_evolve_fg_pure_gauge(gm, gf, ga, 4.0 * ttheta / dt, 0.5 * dt)
        gf_evolve(gf, gm, (1.0 - 2.0 * lam) * dt)
        gm_evolve_fg_pure_gauge(gm, gf, ga, 4.0 * ttheta / dt, 0.5 * dt)
        if i < n_step - 1:
            gf_evolve(gf, gm, 2.0 * lam * dt)
        else:
            gf_evolve(gf, gm, lam * dt)
    gf.unitarize()
    delta_h = gm_hamilton_node(gm) + gf_hamilton_node(gf, ga) - energy
    delta_h = glb_sum_double(delta_h)
    return delta_h

@q.timer(is_timer_fork=True)
def run_hmc_pure_gauge(
    gf,
    ga,
    traj,
    rs,
    *,
    is_reverse_test=False,
    n_step=6,
    md_time=1.0,
    is_always_accept=False,
):
    """Run a single HMC trajectory for pure gauge theory.\n
    Generates random momentum, runs the force-gradient MD evolution,
    and applies the Metropolis accept/reject step.  On accept the input
    ``gf`` is updated in-place.\n
    Parameters
    ----------
    gf : GaugeField
        Input/output gauge field.
    ga : GaugeAction
    traj : int
        Trajectory number (used to split the RNG).
    rs : RngState
    is_reverse_test : bool
        If True, run the reverse trajectory and verify reversibility.
    n_step : int
        Number of MD steps per trajectory (default 6).
    md_time : float
        Total molecular dynamics time (default 1.0).
    is_always_accept : bool
        If True, skip the Metropolis test and always update.\n
    Returns
    -------
    (flag, delta_h)
        ``flag`` is True if the trajectory was accepted.
    """
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
        q.displayln_info(
            f"{fname}: reversed delta_diff: {delta_h + delta_h_rev} / {delta_h}"
        )
        gf_diff_norm = q.qnorm(gf0_r)
        gf_norm = q.qnorm(gf0)
        q.displayln_info(f"{fname}: reversed gf_diff: {gf_diff_norm} / {gf_norm}")
        assert gf_diff_norm <= 1e-12 * gf_norm
    flag, accept_prob = metropolis_accept(delta_h, traj, rs.split("metropolis_accept"))
    if flag or is_always_accept:
        q.displayln_info(f"{fname}: update gf (traj={traj})")
        gf @= gf0
    return flag, delta_h
