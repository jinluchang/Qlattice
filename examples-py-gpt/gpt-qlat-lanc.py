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
        is_test_job_tag,
        )

import pprint
import os
import subprocess

@q.timer
def run_job(job_tag, traj):
    traj_gf = traj
    #
    fns_produce = [
            f"{job_tag}/eig/traj-{traj_gf}/metadata.txt",
            ]
    fns_need = [
            (f"{job_tag}/configs/ckpoint_lat.{traj_gf}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj_gf}",),
            ]
    #
    if is_test_job_tag(job_tag):
        traj_gf = 1000
        fns_need = []
    #
    if not check_job(job_tag, traj, fns_produce, fns_need):
        return
    #
    q.check_stop()
    q.check_time_limit()
    #
    if is_test_job_tag(job_tag):
        total_site = q.Coordinate(get_param(job_tag, "total_site"))
        geo = q.Geometry(total_site)
        q.json_results_append(f"geo.show() = {geo.show()}")
    #
    get_gf = run_gf(job_tag, traj_gf)
    get_gf().show_info()
    #
    get_eig = run_eig(job_tag, traj_gf, get_gf)
    #
    # test repartition
    if is_test_job_tag(job_tag):
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

job_tag = "test-4nt8-checker"
set_param(job_tag, "seed")("test-4nt8")
#
set_param(job_tag, "traj_list")([ 1000, ])
#
set_param(job_tag, "total_site")([ 4, 4, 4, 8, ])
set_param(job_tag, "load_config_params", "twist_boundary_at_boundary")([ 0.0, 0.0, 0.0, -0.5, ])
#
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(2)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(8)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(1)
#
set_param(job_tag, "quark_flavor_list")([ "light", "strange", "charm-1", "charm-2", ])
set_param(job_tag, "quark_mass_list")([ 0.01, 0.04, 0.1, 0.2, ])
set_param(job_tag, "fermion_params", 0, 0)({ 'Ls': 8, 'M5': 1.8, 'b': 1.5, 'c': 0.5, 'boundary_phases': [1.0, 1.0, 1.0, 1.0], })
for inv_type, mass in enumerate(get_param(job_tag, "quark_mass_list")):
    set_param(job_tag, "fermion_params", inv_type, 0)(get_param(job_tag, "fermion_params", 0, 0).copy())
    set_param(job_tag, "fermion_params", inv_type, 0, "mass")(mass)
    for inv_acc in [ 0, 1, 2, ]:
        set_param(job_tag, "fermion_params", inv_type, inv_acc)(get_param(job_tag, "fermion_params", inv_type, 0).copy())
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxiter")(10)
        set_param(job_tag, f"cg_params-{inv_type}-{inv_acc}", "maxcycle")(1 + inv_acc)
#
set_param(job_tag, "lanc_params", 0, 0, "cheby_params")({"low": 0.10, "high": 5.5, "order": 50})
set_param(job_tag, "lanc_params", 0, 0, "irl_params")({ "Nstop": 20, "Nk": 25, "Nm": 30, "resid": 1e-8, "betastp": 0.0, "maxiter": 20, "Nminres": 1, })
set_param(job_tag, "lanc_params", 0, 0, "pit_params")({ "eps": 0.01, "maxiter": 500, "real": True })
set_param(job_tag, "lanc_params", 0, 0, "fermion_params")(get_param(job_tag, "fermion_params", 0, 0).copy())
#
set_param(job_tag, "clanc_params", 0, 0, "nbasis")(20)
set_param(job_tag, "clanc_params", 0, 0, "block")([ 2, 2, 2, 2, ])
set_param(job_tag, "clanc_params", 0, 0, "cheby_params")({ "low": 0.20, "high": 5.5, "order": 50, })
set_param(job_tag, "clanc_params", 0, 0, "save_params")({ "nsingle": 10, "mpi": [ 1, 1, 1, 4, ], })
set_param(job_tag, "clanc_params", 0, 0, "irl_params")({ "Nstop": 30, "Nk": 35, "Nm": 40, "resid": 1e-8, "betastp": 0.0, "maxiter": 20, "Nminres": 1, })
set_param(job_tag, "clanc_params", 0, 0, "smoother_params")({'eps': 1e-08, 'maxiter': 10})

# ----

job_tag_list_default = [
        "test-4nt8-checker",
        ]
job_tag_list_str_default = ",".join(job_tag_list_default)
job_tag_list = q.get_arg("--job_tag_list", default=job_tag_list_str_default).split(",")
if job_tag_list == job_tag_list_default:
    is_test = True
else:
    is_test = False

# ----

def gracefully_finish():
    q.displayln_info("Begin to gracefully_finish.")
    q.timer_display()
    if is_test:
        q.json_results_append(f"q.obtained_lock_history_list={q.obtained_lock_history_list}")
        q.check_log_json(__file__)
    qg.end_with_gpt()
    q.displayln_info("CHECK: finished successfully.")
    exit()

def try_gracefully_finish():
    """
    Call `gracefully_finish` if not test and if some work is done (q.obtained_lock_history_list != [])
    """
    if (not is_test) and (len(q.obtained_lock_history_list) > 0):
        gracefully_finish()

if __name__ == "__main__":

    qg.begin_with_gpt()
    q.check_time_limit()

    job_tag_traj_list = []
    for job_tag in job_tag_list:
        run_params(job_tag)
        traj_list = get_param(job_tag, "traj_list")
        for traj in traj_list:
            job_tag_traj_list.append((job_tag, traj,))
    if not is_test:
        job_tag_traj_list = q.random_permute(job_tag_traj_list, q.RngState(f"{q.get_time()}"))
        job_tag_traj_list = q.get_comm().bcast(job_tag_traj_list)
    for job_tag, traj in job_tag_traj_list:
        q.check_time_limit()
        run_job(job_tag, traj)
        q.clean_cache()
        try_gracefully_finish()

    gracefully_finish()

# ----
