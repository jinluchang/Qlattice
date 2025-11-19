#!/usr/bin/env python3

import sys
import math as m
import numpy as np

import qlat as q

from qlat_scripts.v1 import *

load_path_list[:] = [
        "results",
        ]

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
        flag, delta_h = q.run_hmc_pure_gauge(
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
set_param(job_tag, "hmc", "save_traj_interval")(4)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "16I_b2p8_md4"
set_param(job_tag, "total_site")((16, 16, 16, 32,))
set_param(job_tag, "a_inv_gev")(2.646) # 2003 lattice spacing 0309017.pdf
set_param(job_tag, "hmc", "max_traj")(5000)
set_param(job_tag, "hmc", "max_traj_always_accept")(100)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(4.0)
set_param(job_tag, "hmc", "n_step")(32 * 4)
set_param(job_tag, "hmc", "beta")(2.80)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(4)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "32I_b2p8"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "a_inv_gev")(2.646) # 2003 lattice spacing 0309017.pdf
set_param(job_tag, "hmc", "max_traj")(20000)
set_param(job_tag, "hmc", "max_traj_always_accept")(100)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(1.0)
set_param(job_tag, "hmc", "n_step")(32)
set_param(job_tag, "hmc", "beta")(2.80)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(10)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "32I_b2p8_md2"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "a_inv_gev")(2.646) # 2003 lattice spacing 0309017.pdf
set_param(job_tag, "hmc", "max_traj")(10000)
set_param(job_tag, "hmc", "max_traj_always_accept")(100)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(2.0)
set_param(job_tag, "hmc", "n_step")(32 * 2)
set_param(job_tag, "hmc", "beta")(2.80)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(5)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "32I_b2p8_md3"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "a_inv_gev")(2.646) # 2003 lattice spacing 0309017.pdf
set_param(job_tag, "hmc", "max_traj")(5000)
set_param(job_tag, "hmc", "max_traj_always_accept")(100)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(3.0)
set_param(job_tag, "hmc", "n_step")(32 * 3)
set_param(job_tag, "hmc", "beta")(2.80)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(4)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "32I_b2p8_md4"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "a_inv_gev")(2.646) # 2003 lattice spacing 0309017.pdf
set_param(job_tag, "hmc", "max_traj")(5000)
set_param(job_tag, "hmc", "max_traj_always_accept")(100)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(4.0)
set_param(job_tag, "hmc", "n_step")(32 * 4)
set_param(job_tag, "hmc", "beta")(2.80)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(3)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "32I_b2p8_md5"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "a_inv_gev")(2.646) # 2003 lattice spacing 0309017.pdf
set_param(job_tag, "hmc", "max_traj")(5000)
set_param(job_tag, "hmc", "max_traj_always_accept")(100)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(5.0)
set_param(job_tag, "hmc", "n_step")(32 * 5)
set_param(job_tag, "hmc", "beta")(2.80)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(2)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "32I_b2p8_md6"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "a_inv_gev")(2.646) # 2003 lattice spacing 0309017.pdf
set_param(job_tag, "hmc", "max_traj")(5000)
set_param(job_tag, "hmc", "max_traj_always_accept")(100)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(6.0)
set_param(job_tag, "hmc", "n_step")(32 * 6)
set_param(job_tag, "hmc", "beta")(2.80)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(2)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "32I_b2p8_md7"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "a_inv_gev")(2.646) # 2003 lattice spacing 0309017.pdf
set_param(job_tag, "hmc", "max_traj")(5000)
set_param(job_tag, "hmc", "max_traj_always_accept")(100)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(7.0)
set_param(job_tag, "hmc", "n_step")(32 * 7)
set_param(job_tag, "hmc", "beta")(2.80)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(2)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "32I_b2p8_md8"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "a_inv_gev")(2.646) # 2003 lattice spacing 0309017.pdf
set_param(job_tag, "hmc", "max_traj")(5000)
set_param(job_tag, "hmc", "max_traj_always_accept")(100)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(8.0)
set_param(job_tag, "hmc", "n_step")(32 * 8)
set_param(job_tag, "hmc", "beta")(2.80)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(1)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "32I_b2p8_md9"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "a_inv_gev")(2.646) # 2003 lattice spacing 0309017.pdf
set_param(job_tag, "hmc", "max_traj")(5000)
set_param(job_tag, "hmc", "max_traj_always_accept")(100)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(9.0)
set_param(job_tag, "hmc", "n_step")(32 * 9)
set_param(job_tag, "hmc", "beta")(2.80)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(1)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "32I_b2p8_md10"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "a_inv_gev")(2.646) # 2003 lattice spacing 0309017.pdf
set_param(job_tag, "hmc", "max_traj")(5000)
set_param(job_tag, "hmc", "max_traj_always_accept")(100)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(10.0)
set_param(job_tag, "hmc", "n_step")(32 * 10)
set_param(job_tag, "hmc", "beta")(2.80)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(1)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "32nt48I_b2p95_md5"
set_param(job_tag, "total_site")((32, 32, 32, 48,))
set_param(job_tag, "a_inv_gev")(3.5) # rough guess
set_param(job_tag, "hmc", "max_traj")(5000)
set_param(job_tag, "hmc", "max_traj_always_accept")(10)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(5.0)
set_param(job_tag, "hmc", "n_step")(32 * 5)
set_param(job_tag, "hmc", "beta")(2.95)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(2)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

job_tag = "32I_b2p95_md8"
set_param(job_tag, "total_site")((32, 32, 32, 64,))
set_param(job_tag, "a_inv_gev")(3.5) # rough guess
set_param(job_tag, "hmc", "max_traj")(5000)
set_param(job_tag, "hmc", "max_traj_always_accept")(10)
set_param(job_tag, "hmc", "max_traj_reverse_test")(2)
set_param(job_tag, "hmc", "md_time")(8.0)
set_param(job_tag, "hmc", "n_step")(32 * 8)
set_param(job_tag, "hmc", "beta")(2.95)
set_param(job_tag, "hmc", "c1")(-0.331)
set_param(job_tag, "hmc", "save_traj_interval")(2)
set_param(job_tag, "hmc", "is_saving_topo_info")(True)

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
