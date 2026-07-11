#!/usr/bin/env python3

"""
Field truncation and masked HMC utilities for lattice QCD.

This module provides functions to truncate fields in the time direction
and perform masked HMC evolution where only selected links are updated.

Functions:
    mk_field_truncated: generic field truncation to a time range
    mk_gf_truncated: symmetric gauge field truncation with divisor padding
    mk_gt_truncated: gauge transform truncation
    mk_gf_truncated_evolve: asymmetric truncation with unity-link padding and masked HMC
    mk_selected_points_truncated: point selection truncation
    gf_evolve_masked: evolve gauge field with mask (mask=True links are updated)
    gf_unitarize_masked: unitarize gauge field with mask
    gm_evolve_fg_pure_gauge_masked: accumulate gauge force with mask
    run_hmc_evolve_pure_gauge_masked: masked HMC evolution (lower-level)
    run_hmc_pure_gauge_masked: masked HMC evolution with Metropolis test (higher-level)
"""

import math

import numpy as np
import qlat as q


### -------------------------------------------------------------------
### Field truncation functions
### -------------------------------------------------------------------


@q.timer
def mk_field_truncated(field, *, t_start, t_end):
    """Truncate field to time range [t_start, t_end)."""
    geo = field.geo
    total_site = geo.total_site
    t_size = total_site[3]
    t_size_trunc = (t_end - t_start) % t_size
    if t_size_trunc == 0:
        t_size_trunc = t_size
    assert t_size_trunc <= t_size
    t_size_trunc_valid = t_size_trunc
    size_node_t = geo.size_node[3]
    t_size_trunc = ((t_size_trunc + size_node_t - 1) // size_node_t) * size_node_t
    total_site_trunc = q.Coordinate(
        [total_site[0], total_site[1], total_site[2], t_size_trunc]
    )
    geo_trunc = q.Geometry(total_site_trunc)
    field_trunc = type(field)(geo_trunc, field.multiplicity)
    t_offset = t_start % t_size
    node_site_t = geo.node_site[3]
    n_shift = geo.size_node[3]
    shift = q.Coordinate([0, 0, 0, node_site_t])
    node_site_trunc = geo_trunc.node_site
    coor_node = geo.coor_node
    current_field = field.copy()
    current_xg_field = q.mk_xg_field(geo)
    field_trunc_arr = np.asarray(field_trunc)
    for i_shift in range(n_shift):
        current_field_arr = np.asarray(current_field)
        xg_arr = current_xg_field[:]
        for index in range(current_field.geo.local_volume):
            xg = xg_arr[index]
            t_orig = xg[3] % t_size
            t_local = (t_orig - t_offset) % t_size
            if t_local >= t_size_trunc_valid:
                continue
            xl = [xg[i] - coor_node[i] * node_site_trunc[i] for i in range(3)]
            xl3 = t_local - coor_node[3] * node_site_trunc[3]
            if xl3 < 0 or xl3 >= node_site_trunc[3]:
                continue
            index_trunc = xl[0] + node_site_trunc[0] * (
                xl[1] + node_site_trunc[1] * (xl[2] + node_site_trunc[2] * xl3)
            )
            field_trunc_arr[index_trunc] = current_field_arr[index]
        if i_shift < n_shift - 1:
            current_field = current_field.shift(shift)
            current_xg_field = current_xg_field.shift(shift)
    return field_trunc


@q.timer
def mk_gf_truncated(gf, *, t_center, t_half, t_size_divisor=1):
    """Truncate gauge field symmetrically around t_center with divisor padding."""
    total_site = gf.geo.total_site
    t_size = total_site[3]
    t_size_trunc = 2 * t_half + 1
    r = t_size_trunc % t_size_divisor
    if r != 0:
        t_size_trunc += t_size_divisor - r
    assert t_size_trunc <= t_size
    t_start = (t_center - t_half) % t_size
    t_end = (t_start + t_size_trunc) % t_size
    gf_trunc = mk_field_truncated(gf, t_start=t_start, t_end=t_end)
    geo_trunc = gf_trunc.geo
    gf_arr = np.asarray(gf_trunc)
    xg_arr = q.mk_xg_field(geo_trunc)[:]
    eye3 = np.eye(3, dtype=gf_arr.dtype)
    mask = np.zeros((geo_trunc.local_volume, 4), dtype=bool)
    mask[xg_arr[:, 3] >= 2 * t_half, :] = True
    mask[xg_arr[:, 3] == 2 * t_half - 1, 3] = True
    gf_arr[mask] = eye3
    return gf_trunc, t_start, t_size_trunc


@q.timer
def mk_gt_truncated(gt, *, t_center, t_half, t_size_divisor=1):
    """Truncate gauge transform symmetrically around t_center."""
    total_site = gt.geo.total_site
    t_size = total_site[3]
    t_size_trunc = 2 * t_half + 1
    r = t_size_trunc % t_size_divisor
    if r != 0:
        t_size_trunc += t_size_divisor - r
    assert t_size_trunc <= t_size
    t_start = (t_center - t_half) % t_size
    t_end = (t_start + t_size_trunc) % t_size
    gt_trunc = mk_field_truncated(gt, t_start=t_start, t_end=t_end)
    return gt_trunc, t_start, t_size_trunc


@q.timer(is_verbose=True)
def mk_gf_truncated_evolve(
    gf,
    *,
    t_center,
    t_left,
    t_right,
    t_pad,
    ga,
    rs,
    num_traj,
    n_step,
    md_time,
):
    """Truncate a gauge field in the time direction, fill the padding/gap with
    identity links, then evolve those links via num_traj masked HMC trajectories.

    The result has total time extent ``t_left + t_right + 1 + t_pad`` (rounded
    up to the nearest multiple of ``size_node[3]``).  The valid region
    ``[0, t_left + t_right]`` is a copy of the original field *gf* and is
    frozen during HMC.  The padding region ``[t_left + t_right + 1, ...]`` is
    filled with identity links and then evolved — the boundary temporal link
    (``dir=3`` at slice ``t_left + t_right``) is also identity-initialised and
    evolved, creating a smooth transition between the physical and padding
    regions.

    Parameters
    ----------
    gf : GaugeField
        Source gauge field to truncate.
    t_center : int
        Reference time coordinate for centering the valid region.
    t_left : int
        Number of time slices to include before (and including) *t_center*.
        (suggested default: 3)
    t_right : int
        Number of time slices to include after (and including) *t_center*.
        (suggested default: 7)
    t_pad : int
        Number of additional padding slices beyond the valid region.
        (suggested default: 5)
    ga : GaugeAction
        Gauge action used for the HMC force computation.
        (suggested default: q.GaugeAction(6.0, 0.0))
    rs : RngState
        Random state; split internally before use.
        (suggested default: q.RngState(f"{get_job_seed(job_tag)}-{traj}-{t_center}-{t_left}-{t_right}-{t_pad}"))
    num_traj : int
        Number of HMC trajectories to run (each starts from the previous state).
        (suggested default: 64)
    n_step : int
        Number of molecular-dynamics steps per trajectory.
        (suggested default: 12)
    md_time : float
        Total molecular-dynamics time per trajectory.
        (suggested default: 4.0)

    Returns
    -------
    (gf_trunc, t_start, t_size_trunc)
        The truncated-and-evolved gauge field, the starting time slice, and the
        total time extent of the result."""
    total_site = gf.geo.total_site
    t_size = total_site[3]
    t_size_trunc_valid = t_left + t_right + 1
    t_size_trunc = t_size_trunc_valid + t_pad
    size_node_t = gf.geo.size_node[3]
    t_size_trunc = ((t_size_trunc + size_node_t - 1) // size_node_t) * size_node_t
    assert t_size_trunc <= t_size
    t_start = (t_center - t_left) % t_size
    t_end = (t_start + t_size_trunc) % t_size
    gf_trunc = mk_field_truncated(gf, t_start=t_start, t_end=t_end)
    geo_trunc = gf_trunc.geo
    gf_arr = np.asarray(gf_trunc)
    assert gf_arr.shape == (geo_trunc.local_volume, 4, 3, 3)
    xg_arr = q.mk_xg_field(geo_trunc)[:]
    eye3 = np.eye(3, dtype=gf_arr.dtype)
    mask = np.zeros((geo_trunc.local_volume, 4), dtype=bool)
    mask[xg_arr[:, 3] >= t_size_trunc_valid, :] = True
    mask[xg_arr[:, 3] == t_size_trunc_valid - 1, 3] = True
    gf_arr[mask] = eye3
    rs_hmc = rs.split("mk_gf_truncated_evolve-hmc")
    for traj in range(num_traj):
        run_hmc_pure_gauge_masked(
            gf_trunc,
            ga,
            traj,
            rs_hmc,
            mask,
            n_step=n_step,
            md_time=md_time,
        )
    return gf_trunc, t_start, t_size_trunc


@q.timer
def mk_selected_points_truncated(sp, *, idx_start, idx_end):
    """Truncate SelectedPoints to index range [idx_start, idx_end)."""
    n_keep = idx_end - idx_start
    total_site = sp.psel.total_site
    psel_sub = q.PointsSelection(total_site, n_keep)
    np.asarray(psel_sub)[:] = np.asarray(sp.psel)[idx_start:idx_end]
    sp_trunc = type(sp)(psel_sub)
    sp_trunc @= sp
    return sp_trunc


### -------------------------------------------------------------------
### Masked HMC evolution functions
### Convention: mask=True means the link IS evolved/updated;
###            mask=False means the link is frozen (restored to original
###            after each operation).
### -------------------------------------------------------------------


@q.timer
def gf_evolve_masked(gf, gm, dt, mask):
    """Evolve gf by dt along gm, restore links where mask is False."""
    gf_saved = q.GaugeField(gf.geo)
    gf_saved @= gf
    q.gf_evolve(gf, gm, dt)
    gf_arr = np.asarray(gf)
    gf_saved_arr = np.asarray(gf_saved)
    gf_arr[~mask] = gf_saved_arr[~mask]


@q.timer
def gf_unitarize_masked(gf, mask):
    """Unitarize gf, restore links where mask is False."""
    gf_saved = q.GaugeField(gf.geo)
    gf_saved @= gf
    gf.unitarize()
    gf_arr = np.asarray(gf)
    gf_saved_arr = np.asarray(gf_saved)
    gf_arr[~mask] = gf_saved_arr[~mask]


@q.timer
def gm_evolve_fg_pure_gauge_masked(gm, gf_init, ga, fg_dt, dt, mask):
    """Accumulate gauge force onto gm for links where mask is True."""
    geo = gf_init.geo
    gf = q.GaugeField(geo)
    gf @= gf_init
    gm_force = q.GaugeMomentum(geo)
    q.set_gm_force(gm_force, gf, ga)
    gf_evolve_masked(gf, gm_force, fg_dt, mask)
    q.set_gm_force(gm_force, gf, ga)
    gm_force *= dt
    gm_arr = np.asarray(gm)
    gm_force_arr = np.asarray(gm_force)
    gm_arr[mask] += gm_force_arr[mask]


@q.timer(is_timer_fork=True)
def run_hmc_evolve_pure_gauge_masked(gm, gf, ga, n_step, mask, md_time=1.0):
    """Run masked HMC evolution (lower-level, takes momentum as argument). mask=True links are evolved; mask=False links are frozen."""
    energy = q.gm_hamilton_node(gm) + q.gf_hamilton_node(gf, ga)
    dt = md_time / n_step
    lam = 0.5 * (1.0 - 1.0 / math.sqrt(3.0))
    theta = (2.0 - math.sqrt(3.0)) / 48.0
    ttheta = theta * dt * dt * dt
    gf_evolve_masked(gf, gm, lam * dt, mask)
    for i in range(n_step):
        gm_evolve_fg_pure_gauge_masked(gm, gf, ga, 4.0 * ttheta / dt, 0.5 * dt, mask)
        gf_evolve_masked(gf, gm, (1.0 - 2.0 * lam) * dt, mask)
        gm_evolve_fg_pure_gauge_masked(gm, gf, ga, 4.0 * ttheta / dt, 0.5 * dt, mask)
        if i < n_step - 1:
            gf_evolve_masked(gf, gm, 2.0 * lam * dt, mask)
        else:
            gf_evolve_masked(gf, gm, lam * dt, mask)
    gf_unitarize_masked(gf, mask)
    delta_h = q.gm_hamilton_node(gm) + q.gf_hamilton_node(gf, ga) - energy
    delta_h = q.glb_sum_double(delta_h)
    return delta_h


@q.timer(is_timer_fork=True)
def run_hmc_pure_gauge_masked(
    gf,
    ga,
    traj,
    rs,
    mask,
    *,
    n_step=6,
    md_time=1.0,
):
    """Run masked HMC evolution with Metropolis test (higher-level). mask=True links are evolved; mask=False links are frozen."""
    fname = q.get_fname()
    rs = rs.split(f"{traj}")
    geo = gf.geo
    gf_orig = q.GaugeField(geo)
    gf_orig @= gf
    gf0 = q.GaugeField(geo)
    gf0 @= gf
    gm = q.GaugeMomentum(geo)
    gm.set_rand(rs.split("set_rand_gauge_momentum"), 1.0)
    gm_arr = np.asarray(gm)
    gm_arr[~mask] = 0.0
    delta_h = run_hmc_evolve_pure_gauge_masked(gm, gf0, ga, n_step, mask, md_time)
    q.metropolis_accept(delta_h, traj, rs.split("metropolis_accept"))
    gf0_arr = np.asarray(gf0)
    gf_orig_arr = np.asarray(gf_orig)
    assert (gf0_arr[~mask] == gf_orig_arr[~mask]).all(), (
        "non-masked gf links should be unchanged"
    )
    q.displayln_info(f"{fname}: update gf (traj={traj})")
    gf @= gf0
    return delta_h