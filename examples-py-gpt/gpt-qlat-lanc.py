#!/usr/bin/env python3

# Need --mpi X.X.X.X runtime option

import qlat as q
import gpt as g
import qlat_gpt as qg

from qlat_scripts.v1 import *

import pprint
import os

@q.timer
def test_eig(gf, eig, job_tag, inv_type):
    geo = gf.geo
    src = q.FermionField4d(geo)
    src.set_rand(q.RngState("test_eig:src.set_rand"))
    q.displayln_info(f"CHECK: src norm {src.qnorm():.10E}")
    sol_ref = ru.get_inv(gf, job_tag, inv_type, inv_acc=2, eig=eig, eps=1e-10, mpi_split=False, qtimer=False) * src
    q.displayln_info(f"CHECK: sol_ref norm {sol_ref.qnorm():.10E} with eig")
    for inv_acc in [0, 1, 2]:
        sol = ru.get_inv(gf, job_tag, inv_type, inv_acc, eig=eig, mpi_split=False, qtimer=False) * src
        sol -= sol_ref
        q.displayln_info(f"CHECK: sol diff norm {sol.qnorm():.1E} inv_acc={inv_acc} with eig")
        sol = ru.get_inv(gf, job_tag, inv_type, inv_acc, mpi_split=False, qtimer=False) * src
        sol -= sol_ref
        q.displayln_info(f"CHECK: sol diff norm {sol.qnorm():.1E} inv_acc={inv_acc} without eig")

@q.timer
def run_job(job_tag, traj):
    fns_produce = [
            f"{job_tag}/eig/traj-{traj}/metadata.txt",
            ]
    fns_need = [
            # (f"{job_tag}/configs/ckpoint_lat.{traj}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj}",),
            ]
    if not check_job(job_tag, traj, fns_produce, fns_need):
        return
    #
    traj_gf = traj
    if job_tag[:5] == "test-":
        # ADJUST ME
        traj_gf = 1000
        #
    #
    q.check_stop()
    q.check_time_limit()
    #
    total_site = q.Coordinate(get_param(job_tag, "total_site"))
    geo = q.Geometry(total_site)
    q.displayln_info("CHECK: geo.show() =", geo.show())
    #
    get_gf = run_gf(job_tag, traj_gf)
    get_gf().show_info()
    #
    get_eig = run_eig(job_tag, traj_gf, get_gf)
    test_eig(get_gf(), get_eig(), job_tag, inv_type=0)
    #
    # test repartition
    path = get_load_path(f"{job_tag}/eig/traj-{traj}")
    q.check_compressed_eigen_vectors(path)
    #
    new_size_node = q.Coordinate([ 2, 2, 2, 2, ])
    path_new = path
    q.eigen_system_repartition(new_size_node, path, path_new)
    q.check_compressed_eigen_vectors(path)
    #
    get_eig = run_eig(job_tag, traj_gf, get_gf)
    test_eig(get_gf(), get_eig(), job_tag, inv_type=0)

set_param("test-4nt8", "mk_sample_gauge_field", "rand_n_step")(2)
set_param("test-4nt8", "mk_sample_gauge_field", "flow_n_step")(8)
set_param("test-4nt8", "mk_sample_gauge_field", "hmc_n_traj")(1)
set_param("test-4nt8", "trajs")([ 1000, ])

qg.begin_with_gpt()

q.qremove_all_info("results")

job_tags = [
        "test-4nt8",
        ]

q.check_time_limit()

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(get_param(job_tag)))
    for v in get_param(job_tag).items():
        q.displayln_info(f"CHECK: {v}")
    for traj in get_param(job_tag, "trajs"):
        run_job(job_tag, traj)

q.timer_display()

qg.end_with_gpt()

q.displayln_info(f"CHECK: finished successfully.")
