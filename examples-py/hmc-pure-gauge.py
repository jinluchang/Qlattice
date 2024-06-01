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

@q.timer_verbose
def run_hmc(job_tag):
    fname = q.get_fname()
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    max_traj = get_param(job_tag, "hmc", "max_traj")
    max_traj_always_accept = get_param(job_tag, "hmc", "max_traj_always_accept")
    save_traj_interval = get_param(job_tag, "hmc", "save_traj_interval")
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
        is_always_accept = traj < max_traj_always_accept
        q.run_hmc_pure_gauge(gf, ga, traj, rs.split("run_hmc_pure_gauge"), is_always_accept=is_always_accept)
        traj += 1
        plaq = gf.plaq()
        json_results.append((f"{fname}: {traj} plaq", plaq,))
        if traj % save_traj_interval == 0:
            gf.save(get_save_path(f"{job_tag}/configs/ckpoint_lat.{traj}"))

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 1, 2],
        [1, 1, 1, 4],
        [1, 2, 2, 2],
        [2, 2, 2, 2],
        [2, 2, 2, 4]]

job_tag = "test-4nt8"
set_param(job_tag, "total_site")((4, 4, 4, 8,))
set_param(job_tag, "hmc", "max_traj")(8)
set_param(job_tag, "hmc", "max_traj_always_accept")(4)
set_param(job_tag, "hmc", "beta")(2.13)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(2)

job_tag = "32I-3.5gev"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "hmc", "max_traj")(5000)
set_param(job_tag, "hmc", "max_traj_always_accept")(100)
set_param(job_tag, "hmc", "beta")(3.05)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(10)

# ----

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

    q.qremove_all_info("results")

    for job_tag in job_tags:
        run_params(job_tag)
        run_hmc(job_tag)

    q.check_log_json(__file__, json_results)

    q.timer_display()

    q.end_with_mpi()

    q.displayln_info(f"CHECK: finished successfully.")

# ----
