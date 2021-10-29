#!/usr/bin/env python3

# Need --mpi X.X.X.X --mpi X.X.X runtime option

import qlat as q
import gpt as g
import qlat_gpt as qg
import rbc_ukqcd as ru
import rbc_ukqcd_params as rup
import pprint

import os

import jobs
from jobs import *

jobs.save_path_default = "results"

jobs.load_path_list = [
        "results",
        "../mk-gf-gt/results",
        "../mk-selected-data/results",
        "/sdcc/u/jluchang/qcdqedta/hlbl-data-with-cache",
        ]

@q.timer_verbose
def check_job(job_tag, traj):
    # return True if config is finished or unavailable
    fns_produce = [
            f"psel-prop-rand-u1/{job_tag}/traj={traj}/checkpoint ; type=1.txt",
            f"psel-prop-rand-u1/{job_tag}/traj={traj}/checkpoint ; type=2.txt",
            ]
    is_job_done = True
    for fn in fns_produce:
        if get_load_path(fn) is None:
            q.displayln_info(f"check_job: {job_tag} {traj} to do as file '{fn}' does not exist.")
            is_job_done = False
            break
    if is_job_done:
        return True
    #
    fns_need = [
            f"configs/{job_tag}/ckpoint_lat.{traj}",
            f"point-selection/{job_tag}/traj={traj}.txt",
            ]
    for fn in fns_need:
        if get_load_path(fn) is None:
            q.displayln_info(f"check_job: {job_tag} {traj} unavailable as {fn} does not exist.")
            return True
    #
    q.check_stop()
    q.check_time_limit()
    #
    return False

@q.timer_verbose
def compute_prop_rand_u1(*, job_tag, traj, gf, inv_type, path_sp, psel):
    n_rand_u1 = rup.dict_params[job_tag]["n_rand_u1"]
    inv_acc = 2
    total_site = rup.dict_params[job_tag]["total_site"]
    geo = q.Geometry(total_site, 1)
    rs = q.RngState(f"seed {job_tag} {traj}").split(f"compute_prop_rand_u1 {inv_type}")
    for idx_rand_u1 in range(n_rand_u1):
        tag = f"idx_rand_u1={idx_rand_u1} ; type={inv_type} ; accuracy={inv_acc}"
        fn_sp = os.path.join(path_sp, f"{tag}.lat")
        if get_load_path(fn_sp) is None:
            inv = ru.get_inv(gf, job_tag, inv_type, inv_acc)
            sp_prop = q.mk_rand_u1_prop(inv, geo, psel, rs.split(idx_rand_u1))
            sp_prop.save(get_save_path(fn_sp))
    q.qtouch(get_save_path(os.path.join(path_sp, f"checkpoint ; type={inv_type}.txt")))

@q.timer_verbose
def run_prop_rand_u1_charm(job_tag, traj, get_gf, get_psel):
    inv_type = 2
    if None in [ get_gf, get_psel, ]:
        return
    if get_load_path(f"psel-prop-rand-u1/{job_tag}/traj={traj}/checkpoint ; type={inv_type}.txt") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-rand-u1-charm"):
        gf = get_gf()
        psel = get_psel()
        compute_prop_rand_u1(job_tag = job_tag, traj = traj,
                gf = get_gf(), inv_type = inv_type,
                path_sp = f"psel-prop-rand-u1/{job_tag}/traj={traj}",
                psel = get_psel())
        q.release_lock()

@q.timer_verbose
def run_prop_rand_u1_strange(job_tag, traj, get_gf, get_psel):
    inv_type = 1
    if None in [ get_gf, get_psel, ]:
        return
    if get_load_path(f"psel-prop-rand-u1-strange/{job_tag}/traj={traj}/checkpoint ; type={inv_type}.txt") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-rand-u1-strange"):
        compute_prop_rand_u1(job_tag = job_tag, traj = traj,
                gf = get_gf(), inv_type = inv_type,
                path_sp = f"psel-prop-rand-u1/{job_tag}/traj={traj}",
                psel = get_psel())
        q.release_lock()

@q.timer_verbose
def run_job(job_tag, traj):
    if check_job(job_tag, traj):
        return
    #
    traj_gf = traj
    if job_tag[:5] == "test-":
        # ADJUST ME
        traj_gf = 1000
        #
    #
    get_gf = run_gf(job_tag, traj_gf)
    #
    get_psel = run_psel(job_tag, traj)
    #
    run_prop_rand_u1_charm(job_tag, traj, get_gf, get_psel)
    run_prop_rand_u1_strange(job_tag, traj, get_gf, get_psel)
    #
    q.clean_cache()
    q.timer_display()

rup.dict_params["test-4nt8"]["n_rand_u1"] = 4
rup.dict_params["test-4nt16"]["n_rand_u1"] = 4
rup.dict_params["48I"]["n_rand_u1"] = 4
rup.dict_params["64I"]["n_rand_u1"] = 4

rup.dict_params["test-4nt8"]["trajs"] = list(range(1000, 1400, 100))
rup.dict_params["test-4nt16"]["trajs"] = list(range(1000, 1400, 100))
rup.dict_params["48I"]["trajs"] = list(range(3000, 500, -5))
rup.dict_params["64I"]["trajs"] = list(range(3000, 500, -5))

# rup.dict_params["test-4nt8"]["fermion_params"][0][2]["Ls"] = 10
# rup.dict_params["test-4nt8"]["fermion_params"][1][2]["Ls"] = 10

# rup.dict_params["test-4nt16"]["fermion_params"][0][2]["Ls"] = 10
# rup.dict_params["test-4nt16"]["fermion_params"][1][2]["Ls"] = 10

qg.begin_with_gpt()

# ADJUST ME
job_tags = [
        "test-4nt8", "test-4nt16",
        # "64I",
        # "48I",
        # "24D",
        ]

q.check_time_limit()

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        run_job(job_tag, traj)

qg.end_with_gpt()
