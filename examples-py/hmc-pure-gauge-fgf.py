#!/usr/bin/env python3

"""
Flexible gauge fixed HMC
"""

import sys
import math as m
import numpy as np

import qlat as q

from qlat_scripts.v1 import *

load_path_list[:] = [
        "results",
        ]

@q.timer
def gm_evolve_fg_pure_gauge(gm, gf_init, ga, fg_dt, dt):
    geo = gf_init.geo
    gf = q.GaugeField(geo)
    gf @= gf_init
    gm_force = q.GaugeMomentum(geo)
    q.set_gm_force(gm_force, gf, ga)
    q.gf_evolve(gf, gm_force, fg_dt)
    q.set_gm_force(gm_force, gf, ga)
    gm_force *= dt
    gm += gm_force

@q.timer_verbose
def run_hmc_evolve_pure_gauge(gm, gf, ga, rs, n_step, md_time=1.0):
    energy = q.gm_hamilton_node(gm) + q.gf_hamilton_node(gf, ga)
    dt = md_time / n_step
    lam = 0.5 * (1.0 - 1.0 / math.sqrt(3.0));
    theta = (2.0 - math.sqrt(3.0)) / 48.0;
    ttheta = theta * dt * dt * dt;
    q.gf_evolve(gf, gm, lam * dt)
    for i in range(n_step):
        gm_evolve_fg_pure_gauge(gm, gf, ga, 4.0 * ttheta / dt, 0.5 * dt);
        q.gf_evolve(gf, gm, (1.0 - 2.0 * lam) * dt);
        gm_evolve_fg_pure_gauge(gm, gf, ga, 4.0 * ttheta / dt, 0.5 * dt);
        if i < n_step - 1:
            q.gf_evolve(gf, gm, 2.0 * lam * dt);
        else:
            q.gf_evolve(gf, gm, lam * dt);
    gf.unitarize()
    delta_h = q.gm_hamilton_node(gm) + q.gf_hamilton_node(gf, ga) - energy;
    delta_h = q.glb_sum_double(delta_h)
    return delta_h

@q.timer_verbose
def run_hmc_pure_gauge(gf, ga, traj, rs, *, is_reverse_test=False, n_step=6, md_time=1.0, is_always_accept=False):
    fname = q.get_fname()
    rs = rs.split(f"{traj}")
    geo = gf.geo
    gf0 = q.GaugeField(geo)
    gf0 @= gf
    gm = q.GaugeMomentum(geo)
    gm.set_rand(rs.split("set_rand_gauge_momentum"), 1.0)
    delta_h = run_hmc_evolve_pure_gauge(gm, gf0, ga, rs, n_step, md_time)
    if is_reverse_test:
        gm_r = q.GaugeMomentum(geo)
        gm_r @= gm
        gf0_r = q.GaugeField(geo)
        gf0_r @= gf0
        delta_h_rev = run_hmc_evolve_pure_gauge(gm_r, gf0_r, ga, rs, n_step, -md_time)
        gf0_r -= gf
        q.displayln_info(f"{fname}: reversed delta_diff: {delta_h + delta_h_rev} / {delta_h}")
        gf_diff_norm = q.qnorm(gf0_r)
        gf_norm = q.qnorm(gf0)
        q.displayln_info(f"{fname}: reversed gf_diff: {gf_diff_norm} / {gf_norm}")
        assert gf_diff_norm <= 1e-12 * gf_norm
    flag, accept_prob = q.metropolis_accept(delta_h, traj, rs.split("metropolis_accept"))
    if flag or is_always_accept:
        q.displayln_info(f"{fname}: update gf (traj={traj})")
        gf @= gf0
    return flag, delta_h

@q.timer_verbose
def run_topo_info(job_tag, traj, gf):
    info_path = get_save_path(f"{job_tag}/topo-measure-wilson-flow/traj-{traj}")
    flow_time = 6
    flow_n_step = 80
    smear_info_list = [
            [ 1.0 / flow_n_step, flow_n_step, 0.0, "runge-kutta", ],
            ] * flow_time
    energy_derivative_info = [ 1.0 / flow_n_step, 0.0, "runge-kutta", ]
    topo_list, energy_list, = q.smear_measure_topo(
            gf.copy(),
            smear_info_list=smear_info_list,
            energy_derivative_info=energy_derivative_info,
            info_path=info_path,
            density_field_path=info_path,
            )

@q.timer_verbose
def run_hmc(job_tag):
    fname = q.get_fname()
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    max_traj = get_param(job_tag, "hmc", "max_traj")
    max_traj_always_accept = get_param(job_tag, "hmc", "max_traj_always_accept")
    max_traj_reverse_test= get_param(job_tag, "hmc", "max_traj_reverse_test")
    save_traj_interval = get_param(job_tag, "hmc", "save_traj_interval")
    is_saving_topo_info = get_param(job_tag, "hmc", "is_saving_topo_info")
    md_time = get_param(job_tag, "hmc", "md_time")
    n_step = get_param(job_tag, "hmc", "n_step")
    beta = get_param(job_tag, "hmc", "beta")
    c1 = get_param(job_tag, "hmc", "c1")
    ga = q.GaugeAction(beta, c1)
    geo = q.Geometry(total_site)
    rs = q.RngState(f"run_hmc-{get_job_seed(job_tag)}")
    gf = q.GaugeField(geo)
    traj_load = None
    if get_load_path(f"{job_tag}/configs") is not None:
        for traj in range(max_traj):
            fn = get_load_path(f"{job_tag}/configs/ckpoint_lat.{traj}")
            if fn is not None:
                traj_load = traj
    if traj_load is None:
        traj = 0
        gf.set_rand(rs.split("init"), 0.1, 2)
        gf.save(get_save_path(f"{job_tag}/configs/ckpoint_lat.{traj}"))
    else:
        traj = traj_load
        gf.load(get_load_path(f"{job_tag}/configs/ckpoint_lat.{traj}"))
    for traj in range(traj, max_traj):
        traj += 1
        is_always_accept = traj < max_traj_always_accept
        is_reverse_test = traj < max_traj_reverse_test
        flag, delta_h = run_hmc_pure_gauge(
                gf, ga, traj,
                rs.split("run_hmc_pure_gauge"),
                n_step=n_step,
                md_time=md_time,
                is_always_accept=is_always_accept,
                is_reverse_test=is_reverse_test,
                )
        plaq = gf.plaq()
        info = dict()
        info["traj"] = traj
        info["plaq"] = plaq
        info["flag"] = flag
        info["delta_h"] = delta_h
        q.qtouch_info(get_save_path(f"{job_tag}/configs/ckpoint_lat_info.{traj}.txt"), pformat(info))
        q.json_results_append(f"{fname}: {traj} plaq", plaq)
        if traj % save_traj_interval == 0:
            gf.save(get_save_path(f"{job_tag}/configs/ckpoint_lat.{traj}"))
            if is_saving_topo_info:
                run_topo_info(job_tag, traj, gf)
        q.timer_display()

# ----

job_tag = "test-4nt8"
set_param(job_tag, "total_site")((4, 4, 4, 8,))
set_param(job_tag, "hmc", "max_traj")(8)
set_param(job_tag, "hmc", "max_traj_always_accept")(4)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(1.0)
set_param(job_tag, "hmc", "n_step")(6)
set_param(job_tag, "hmc", "beta")(2.13)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "gauge_fixing", "block_site")((4, 4, 4, 4,))
set_param(job_tag, "hmc", "gauge_fixing", "stout_smear_step_size")(0.125)
set_param(job_tag, "hmc", "gauge_fixing", "num_smear_step")(4)
set_param(job_tag, "hmc", "save_traj_interval")(4)
set_param(job_tag, "hmc", "is_saving_topo_info")(False)

# ----

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 1, 2],
        [1, 1, 1, 3],
        [1, 1, 1, 4],
        [1, 1, 1, 6],
        [1, 1, 1, 8],
        [1, 2, 2, 4],
        [2, 2, 2, 4],
        [2, 2, 2, 4],
        ]

if __name__ == "__main__":

    q.begin_with_mpi(size_node_list)

    ##################### CMD options #####################

    job_tag_list = q.get_arg("--job_tag_list", default="").split(",")

    #######################################################

    job_tag_list_default = [
            "test-4nt8",
            ]

    if job_tag_list == [ "", ]:
        job_tag_list = job_tag_list_default

    for job_tag in job_tag_list:
        run_params(job_tag)
        run_hmc(job_tag)

    q.check_log_json(__file__)

    q.timer_display()

    q.end_with_mpi()

    q.displayln_info(f"CHECK: finished successfully.")

# ----
