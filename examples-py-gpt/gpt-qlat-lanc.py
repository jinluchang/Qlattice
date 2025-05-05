#!/usr/bin/env python3

# Need --mpi X.X.X.X runtime option

import qlat as q
import gpt as g
import qlat_gpt as qg

from qlat_scripts.v1 import (
        get_load_path,
        set_param,
        get_param,
        check_job,
        run_gf,
        run_eig,
        run_params,
        test_eig,
        )

import pprint
import os
import subprocess

@q.timer
def run_job(job_tag, traj):
    is_test = job_tag[:5] == "test-"
    #
    traj_gf = traj
    #
    fns_produce = [
            f"{job_tag}/eig/traj-{traj_gf}/metadata.txt",
            ]
    fns_need = [
            (f"{job_tag}/configs/ckpoint_lat.{traj_gf}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj_gf}",),
            ]
    #
    if is_test:
        traj_gf = 1000
        fns_need = []
    #
    if not check_job(job_tag, traj, fns_produce, fns_need):
        return
    #
    q.check_stop()
    q.check_time_limit()
    #
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    q.json_results_append(f"geo.show() = {geo.show()}")
    #
    get_gf = run_gf(job_tag, traj_gf)
    get_gf().show_info()
    #
    get_eig = run_eig(job_tag, traj_gf, get_gf)
    test_eig(get_gf(), get_eig(), job_tag, inv_type=0)
    #
    # test repartition
    path = get_load_path(f"{job_tag}/eig/traj-{traj_gf}")
    q.check_compressed_eigen_vectors(path)
    #
    new_size_node = q.Coordinate([ 2, 2, 2, 2, ])
    path_new = path
    q.eigen_system_repartition(new_size_node, path, path_new)
    q.check_compressed_eigen_vectors(path)
    #
    get_eig = run_eig(job_tag, traj_gf, get_gf)
    test_eig(get_gf(), get_eig(), job_tag, inv_type=0)
    #
    # test zip folder
    q.sync_node()
    for i in range(32):
        if i % q.get_num_node() == q.get_id_node():
            if q.is_directory(f"{path}/{i:02}"):
                subprocess.run([ "zip", "-r", "-Z", "store", f"{i:02}.zip", f"{i:02}" ], cwd=path)
                subprocess.run([ "rm", "-rf", f"{i:02}/" ], cwd=path)
    q.sync_node()
    #
    get_eig = run_eig(job_tag, traj_gf, get_gf)
    test_eig(get_gf(), get_eig(), job_tag, inv_type=0)

# ----

set_param("test-4nt8", "trajs")([ 1000, ])
set_param("test-4nt8", "mk_sample_gauge_field", "rand_n_step")(2)
set_param("test-4nt8", "mk_sample_gauge_field", "flow_n_step")(8)
set_param("test-4nt8", "mk_sample_gauge_field", "hmc_n_traj")(1)
for inv_type in [ 0, 1, 2, ]:
    for inv_acc in [ 0, 1, 2, ]:
        set_param("test-4nt8", "fermion_params", inv_type, inv_acc, "Ls")(8)
        set_param("test-4nt8", f"cg_params-{inv_type}-{inv_acc}", "maxiter")(10)
        set_param("test-4nt8", f"cg_params-{inv_type}-{inv_acc}", "maxcycle")(1 + inv_acc)

# ----

def gracefully_finish():
    q.displayln_info("Begin to gracefully_finish.")
    q.timer_display()
    qg.end_with_gpt()
    q.displayln_info("CHECK: finished successfully.")
    exit()

if __name__ == "__main__":

    qg.begin_with_gpt()

    job_tags = q.get_arg("--job_tags", default="").split(",")

    job_tags_default = [
            "test-4nt8",
            ]

    if job_tags == [ "", ]:
        job_tags = job_tags_default

    q.check_time_limit()

    for job_tag in job_tags:
        run_params(job_tag)
        for traj in get_param(job_tag, "trajs"):
            q.check_time_limit()
            run_job(job_tag, traj)
            q.clean_cache()
            if q.obtained_lock_history_list:
                q.json_results_append(f"q.obtained_lock_history_list={q.obtained_lock_history_list}")
                if job_tag[:5] != "test-":
                    gracefully_finish()

    q.check_log_json(__file__, check_eps=1e-5)

    gracefully_finish()

# ----
