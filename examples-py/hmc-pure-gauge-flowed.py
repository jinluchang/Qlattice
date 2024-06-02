#!/usr/bin/env python3

json_results = []

import sys
import math as m
import numpy as np

import qlat as q

from qlat_scripts.v1 import *

load_path_list[:] = [
        "results",
        ]

from qlat import metropolis_accept

@q.timer_verbose
def gm_evolve_fg(gm, gf_init, ga, fi, fg_dt, dt):
    geo = gf_init.geo
    gf = q.GaugeField(geo)
    gf @= gf_init
    gm_force = q.GaugeMomentum(geo)
    q.set_gm_force_flowed(gm_force, gf, ga, fi)
    q.gf_evolve(gf, gm_force, fg_dt)
    q.set_gm_force_flowed(gm_force, gf, ga, fi)
    gm_force *= dt
    gm += gm_force

@q.timer_verbose
def run_hmc_evolve_flowed(gm, gf, ga, fi, rs, steps, md_time=1.0):
    energy = q.gm_hamilton_node(gm) + q.gf_hamilton_flowed_node(gf, ga, fi)
    dt = md_time / steps
    lam = 0.5 * (1.0 - 1.0 / m.sqrt(3.0));
    theta = (2.0 - m.sqrt(3.0)) / 48.0;
    ttheta = theta * dt * dt * dt;
    q.gf_evolve(gf, gm, lam * dt)
    for i in range(steps):
        gm_evolve_fg(gm, gf, ga, fi, 4.0 * ttheta / dt, 0.5 * dt);
        q.gf_evolve(gf, gm, (1.0 - 2.0 * lam) * dt);
        gm_evolve_fg(gm, gf, ga, fi, 4.0 * ttheta / dt, 0.5 * dt);
        if i < steps - 1:
            q.gf_evolve(gf, gm, 2.0 * lam * dt);
        else:
            q.gf_evolve(gf, gm, lam * dt);
    q.unitarize(gf)
    delta_h = q.gm_hamilton_node(gm) + q.gf_hamilton_flowed_node(gf, ga, fi) - energy;
    delta_h = q.glb_sum(delta_h)
    return delta_h

@q.timer_verbose
def mk_flow_info(fp, rng):
    time = fp["time"]
    fi = q.FlowInfo()
    fi.add_rand_order_flow(rng, time)
    # fi.add_rand_order_flow(rng, 0.1, -0.01)
    # fi.add_rand_order_flow(rng, 0.1, 0.0)
    return fi

@q.timer_verbose
def run_hmc_pure_gauge(gf, ga, fp, traj, rs, *, is_reverse_test=False, n_step=6, md_time=1.0, is_always_accept=False):
    fname = q.get_fname()
    rs = rs.split(f"{traj}")
    fi = mk_flow_info(fp, rs.split("mk_flow_info"))
    geo = gf.geo
    gf0 = q.GaugeField(geo)
    q.gf_flow_inv(gf0, gf, fi)
    if is_reverse_test:
        gf_r = q.GaugeField(geo)
        q.gf_flow(gf_r, gf0, fi)
        gf_r -= gf
        q.displayln_info(f"gf_flow_inv gf_diff: {q.qnorm(gf_r)} / {q.qnorm(gf)}")
    gm = q.GaugeMomentum(geo)
    gm.set_rand(rs.split("set_rand_gauge_momentum"), 1.0)
    delta_h = run_hmc_evolve_flowed(gm, gf0, ga, fi, rs, n_step, md_time)
    if is_reverse_test:
        gm_r = q.GaugeMomentum(geo)
        gm_r @= gm
        gf0_r = q.GaugeField(geo)
        gf0_r @= gf0
        delta_h_rev = run_hmc_evolve_flowed(gm_r, gf0_r, ga, fi, rs, n_step, -md_time)
        gf0_r -= gf;
        q.displayln_info(f"{fname}: reversed delta_diff: {delta_h + delta_h_rev} / {delta_h}")
        q.displayln_info(f"{fname}: reversed gf_diff: {q.qnorm(gf0_r)} / {q.qnorm(gf0)}")
    flag, accept_prob = metropolis_accept(delta_h, traj, rs.split("metropolis_accept"))
    if flag or is_always_accept:
        q.displayln_info(f"{fname}: update gf (traj={traj})")
        q.gf_flow(gf, gf0, fi)

@q.timer_verbose
def run_hmc(job_tag):
    fname = q.get_fname()
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    max_traj = get_param(job_tag, "hmc", "max_traj")
    max_traj_always_accept = get_param(job_tag, "hmc", "max_traj_always_accept")
    save_traj_interval = get_param(job_tag, "hmc", "save_traj_interval")
    md_time = get_param(job_tag, "hmc", "md_time")
    n_step = get_param(job_tag, "hmc", "n_step")
    fp = get_param(job_tag, "hmc", "fp") # field transformation (flow) parameters
    beta = get_param(job_tag, "hmc", "beta")
    c1 = get_param(job_tag, "hmc", "c1")
    ga = q.GaugeAction(beta, c1)
    geo = q.Geometry(total_site)
    rs = q.RngState(f"run_hmc-{job_tag}")
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
        run_hmc_pure_gauge(gf, ga, fp, traj, rs.split("run_hmc_pure_gauge"), n_step=n_step, md_time=md_time, is_always_accept=is_always_accept)
        plaq = gf.plaq()
        json_results.append((f"{fname}: {traj} plaq", plaq,))
        if traj % save_traj_interval == 0:
            gf.save(get_save_path(f"{job_tag}/configs/ckpoint_lat.{traj}"))

job_tag = "test-4nt8"
set_param(job_tag, "total_site")((4, 4, 4, 8,))
set_param(job_tag, "hmc", "max_traj")(8)
set_param(job_tag, "hmc", "max_traj_always_accept")(4)
set_param(job_tag, "hmc", "md_time")(1.0)
set_param(job_tag, "hmc", "n_step")(6)
set_param(job_tag, "hmc", "beta")(2.13)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "fp", "time")(0.1)
set_param(job_tag, "hmc", "save_traj_interval")(2)

job_tag = "32I-3.5gev-ft"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "hmc", "max_traj")(5000)
set_param(job_tag, "hmc", "max_traj_always_accept")(100)
set_param(job_tag, "hmc", "md_time")(1.0)
set_param(job_tag, "hmc", "n_step")(32)
set_param(job_tag, "hmc", "beta")(3.05)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "fp", "time")(0.1)
set_param(job_tag, "hmc", "save_traj_interval")(10)

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

    job_tags = q.get_arg("--job_tags", default="").split(",")

    #######################################################

    job_tags_default = [
            "test-4nt8",
            ]

    if job_tags == [ "", ]:
        job_tags = job_tags_default

    for job_tag in job_tags:
        run_params(job_tag)
        run_hmc(job_tag)

    q.check_log_json(__file__, json_results)

    q.timer_display()

    q.end_with_mpi()

    q.displayln_info(f"CHECK: finished successfully.")

# ----
