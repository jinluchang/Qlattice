#!/usr/bin/env python3

import sys, os
import numpy as np
import qlat as q

from qlat import (
        run_flow_scale,
        )

from qlat_scripts.v1 import (
        set_param,
        get_param,
        load_path_list,
        run_params,
        run_gf,
        check_job,
        get_save_path,
        get_load_path,
        )

load_path_list[:] = [
        "results",
        "qcddata-1",
        "qcddata-2",
        "qcddata-3",
        "qcddata-4",
        "qcddata-5",
        "qcddata-6",
        "qcddata-7",
        "qcddata-8",
        "qcddata-9",
        ]

# --------------------------------------------

@q.timer(is_timer_fork=True)
def run_job(job_tag, traj):
    step_size = get_param(job_tag, "flow_scale", "step_size")
    num_step = get_param(job_tag, "flow_scale", "num_step")
    t_dir_list = get_param(job_tag, "flow_scale", "t_dir_list")
    integrator_type = get_param(job_tag, "flow_scale", "integrator_type")
    #
    traj_gf = traj
    #
    fns_produce = [
            f"{job_tag}/gf-flow-record/traj-{traj_gf}.pickle",
            ]
    for t_dir in t_dir_list:
        if t_dir == 3:
            fns_produce += [
                    f"{job_tag}/gf-flow-record-spatial/traj-{traj_gf}.pickle",
                    ]
        else:
            fns_produce += [
                    f"{job_tag}/gf-flow-record-spatial-t_dir-{t_dir}/traj-{traj_gf}.pickle",
                    ]
    fns_need = [
            (f"{job_tag}/configs/ckpoint_lat.{traj_gf}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj_gf}",),
            ]
    #
    if q.is_test():
        fns_need = []
    #
    if not check_job(job_tag, traj, fns_produce, fns_need):
        return
    #
    q.check_stop()
    q.check_time_limit()
    #
    get_gf = run_gf(job_tag, traj)
    #
    fn_out = f"{job_tag}/gf-flow-record/traj-{traj_gf}.pickle"
    params = q.default_run_flow_scale_params | dict(
            step_size=step_size,
            num_step=num_step,
            integrator_type=integrator_type,
            is_spatial=False,
            )
    if get_load_path(fn_out) is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-gf-flow-record"):
            run_flow_scale(get_save_path(fn_out), get_gf=get_gf, params=params)
            q.release_lock()
    #
    for t_dir in t_dir_list:
        if t_dir == 3:
            fn_out = f"{job_tag}/gf-flow-record-spatial/traj-{traj_gf}.pickle"
        else:
            fn_out = f"{job_tag}/gf-flow-record-spatial-t_dir-{t_dir}/traj-{traj_gf}.pickle"
        params = q.default_run_flow_scale_params | dict(
                step_size=step_size,
                num_step=num_step,
                integrator_type=integrator_type,
                is_spatial=True,
                t_dir=t_dir,
                )
        if get_load_path(fn_out) is None:
            if q.obtain_lock(f"locks/{job_tag}-{traj}-gf-flow-record-spatial-t_dir-{t_dir}"):
                run_flow_scale(get_save_path(fn_out), get_gf=get_gf, params=params)
                q.release_lock()

# --------------------------------------------

job_tag = "test-4nt8-checker"
set_param(job_tag, "traj_list")([ 1000, 1100, ])
set_param(job_tag, "seed")("test-4nt8")
set_param(job_tag, "total_site")([ 4, 4, 4, 8, ])
set_param(job_tag, "load_config_params")(None)
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(2)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(8)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(1)
set_param(job_tag, "flow_scale")(dict(
    step_size=0.05,
    num_step=10,
    t_dir_list=[ 0, 3, ],
    integrator_type="runge-kutta",
    ))

job_tag = "16IH2"
set_param(job_tag, "traj_list")(list(range(1000, 4500, 10)))
set_param(job_tag, "load_config_params")(None)
set_param(job_tag, "flow_scale")(dict(
    step_size=0.05,
    num_step=400,
    t_dir_list=[ 3, ],
    integrator_type="runge-kutta",
    ))

job_tag = "48I"
set_param(job_tag, "traj_list")(list(range(1000, 3000, 10)))
set_param(job_tag, "load_config_params")(None)
set_param(job_tag, "flow_scale")(dict(
    step_size=0.05,
    num_step=400,
    t_dir_list=[ 3, ],
    integrator_type="runge-kutta",
    ))

job_tag = "64I"
set_param(job_tag, "traj_list")(list(range(1000, 4000, 10)))
set_param(job_tag, "load_config_params")(None)
set_param(job_tag, "flow_scale")(dict(
    step_size=0.05,
    num_step=800,
    t_dir_list=[ 3, ],
    integrator_type="runge-kutta",
    ))

job_tag = "96I"
set_param(job_tag, "traj_list")(list(range(800, 2000, 10)))
set_param(job_tag, "load_config_params")(None)
set_param(job_tag, "flow_scale")(dict(
    step_size=0.05,
    num_step=1000,
    t_dir_list=[ 3, ],
    integrator_type="runge-kutta",
    ))

job_tag = "48If"
set_param(job_tag, "traj_list")(list(range(1000, 2000, 10)))
set_param(job_tag, "load_config_params")(None)
set_param(job_tag, "flow_scale")(dict(
    step_size=0.05,
    num_step=1000,
    t_dir_list=[ 3, ],
    integrator_type="runge-kutta",
    ))

job_tag = "9"
set_param(job_tag, "traj_list")(list(range(500, 3000, 10)))
set_param(job_tag, "load_config_params")(None)
set_param(job_tag, "flow_scale")(dict(
    step_size=0.05,
    num_step=800,
    t_dir_list=[ 3, ],
    integrator_type="runge-kutta",
    ))

# --------------------------------------------

job_tag_list_default = [
        "test-4nt8-checker",
        ]
job_tag_list_str_default = ",".join(job_tag_list_default)
job_tag_list = q.get_arg("--job_tag_list", default=job_tag_list_str_default).split(",")

# ----

def gracefully_finish():
    q.displayln_info("Begin to gracefully_finish.")
    q.timer_display()
    if q.is_test():
        q.json_results_append(f"q.obtained_lock_history_list={q.obtained_lock_history_list}")
        q.check_log_json(__file__)
    q.end_with_mpi()
    q.displayln_info("CHECK: finished successfully.")
    exit()

def try_gracefully_finish():
    """
    Call `gracefully_finish` if not test and if some work is done (q.obtained_lock_history_list != [])
    """
    if (not q.is_test()) and (len(q.obtained_lock_history_list) > 0):
        gracefully_finish()

size_node_list = [
        [ 1, 1, 1, 1, ],
        [ 1, 1, 1, 2, ],
        [ 1, 1, 1, 4, ],
        ]

if __name__ == "__main__":

    q.begin_with_mpi(size_node_list)
    q.check_time_limit()

    job_tag_traj_list = []
    for job_tag in job_tag_list:
        run_params(job_tag)
        traj_list = get_param(job_tag, "traj_list")
        for traj in traj_list:
            job_tag_traj_list.append((job_tag, traj,))
    if not q.is_test():
        job_tag_traj_list = q.random_permute(job_tag_traj_list, q.RngState(f"{q.get_time()}"))
        job_tag_traj_list = q.get_comm().bcast(job_tag_traj_list)
    for job_tag, traj in job_tag_traj_list:
        q.check_time_limit()
        run_job(job_tag, traj)
        q.clean_cache()
        try_gracefully_finish()

    gracefully_finish()

# ----
