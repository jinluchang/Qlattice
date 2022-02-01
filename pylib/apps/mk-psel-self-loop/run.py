#!/usr/bin/env python3

# Need --mpi X.X.X.X --mpi X.X.X runtime option

import qlat as q
import gpt as g
import qlat_gpt as qg
import rbc_ukqcd as ru
import rbc_ukqcd_params as rup
import pprint

import os

from jobs import *

load_path_list[:] = [
        "results",
        "../mk-gf-gt/results",
        "../mk-sel/results",
        "../mk-lanc/results",
        "../qcddata",
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/mk-gf-gt/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/mk-sel/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/mk-lanc/results"),
        os.path.join(os.getenv("HOME"), "qcddata"),
        "../mk-selected-data/results",
        "/sdcc/u/jluchang/qcdqedta/hlbl-data-with-cache",
        ]

@q.timer_verbose
def compute_prop_rand_u1(*, job_tag, traj, gf, inv_type, path_sp, psel, eig = None):
    n_rand_u1 = rup.dict_params[job_tag]["n_rand_u1"]
    inv_acc = 2
    total_site = rup.dict_params[job_tag]["total_site"]
    geo = q.Geometry(total_site, 1)
    rs = q.RngState(f"seed {job_tag} {traj}").split(f"compute_prop_rand_u1 {inv_type}")
    for idx_rand_u1 in range(n_rand_u1):
        tag = f"idx_rand_u1={idx_rand_u1} ; type={inv_type} ; accuracy={inv_acc}"
        fn_sp = os.path.join(path_sp, f"{tag}.lat")
        if get_load_path(fn_sp) is None:
            q.check_stop()
            q.check_time_limit()
            inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, eig = eig)
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
    if get_load_path(f"psel-prop-rand-u1/{job_tag}/traj={traj}/checkpoint ; type={inv_type}.txt") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-rand-u1-strange"):
        compute_prop_rand_u1(job_tag = job_tag, traj = traj,
                gf = get_gf(), inv_type = inv_type,
                path_sp = f"psel-prop-rand-u1/{job_tag}/traj={traj}",
                psel = get_psel())
        q.release_lock()

@q.timer_verbose
def run_prop_rand_u1_light(job_tag, traj, get_gf, get_eig, get_psel):
    inv_type = 0
    if None in [ get_gf, get_eig, get_psel, ]:
        return
    if get_load_path(f"psel-prop-rand-u1/{job_tag}/traj={traj}/checkpoint ; type={inv_type}.txt") is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-rand-u1-light"):
        gf = get_gf()
        eig = get_eig()
        compute_prop_rand_u1(job_tag = job_tag, traj = traj,
                gf = gf, inv_type = inv_type,
                path_sp = f"psel-prop-rand-u1/{job_tag}/traj={traj}",
                psel = get_psel(), eig = eig)
        q.release_lock()

@q.timer_verbose
def run_job(job_tag, traj):
    fns_produce = [
            f"psel-prop-rand-u1/{job_tag}/traj={traj}/checkpoint ; type=0.txt",
            f"psel-prop-rand-u1/{job_tag}/traj={traj}/checkpoint ; type=1.txt",
            f"psel-prop-rand-u1/{job_tag}/traj={traj}/checkpoint ; type=2.txt",
            ]
    fns_need = [
            (f"configs/{job_tag}/ckpoint_lat.{traj}", f"configs/{job_tag}/ckpoint_lat.IEEE64BIG.{traj}",),
            f"point-selection/{job_tag}/traj={traj}.txt",
            f"eig/{job_tag}/traj={traj}",
            f"eig/{job_tag}/traj={traj}/metadata.txt",
            f"eig/{job_tag}/traj={traj}/eigen-values.txt",
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
    get_gf = run_gf(job_tag, traj_gf)
    #
    get_psel = run_psel(job_tag, traj)
    #
    def run_with_eig():
        get_eig = run_eig(job_tag, traj_gf, get_gf)
        run_prop_rand_u1_light(job_tag, traj, get_gf, get_eig, get_psel)
    #
    run_with_eig()
    #
    run_prop_rand_u1_strange(job_tag, traj, get_gf, get_psel)
    run_prop_rand_u1_charm(job_tag, traj, get_gf, get_psel)
    #
    q.clean_cache()
    q.timer_display()

# rup.dict_params["test-4nt8"]["fermion_params"][0][2]["Ls"] = 10
# rup.dict_params["test-4nt8"]["fermion_params"][1][2]["Ls"] = 10

# rup.dict_params["test-4nt16"]["fermion_params"][0][2]["Ls"] = 10
# rup.dict_params["test-4nt16"]["fermion_params"][1][2]["Ls"] = 10

rup.dict_params["test-4nt8"]["n_rand_u1"] = 4
rup.dict_params["test-4nt16"]["n_rand_u1"] = 4

rup.dict_params["test-4nt8"]["trajs"] = list(range(1000, 1400, 100))
rup.dict_params["test-4nt16"]["trajs"] = list(range(1000, 1400, 100))

rup.dict_params["48I"]["n_rand_u1"] = 2
rup.dict_params["64I"]["n_rand_u1"] = 2

rup.dict_params["16IH2"]["n_rand_u1"] = 2
rup.dict_params["32IfineH"]["n_rand_u1"] = 2

for inv_acc in [ 0, 1, 2, ]:
    rup.dict_params["64I"]["fermion_params"][0][inv_acc]["mass"] = 0.0006203
    rup.dict_params["64I"]["fermion_params"][1][inv_acc]["mass"] = 0.02539

for inv_acc in [ 0, 1, 2, ]:
    rup.dict_params["48I"]["fermion_params"][0][inv_acc]["mass"] = 0.0006979
    rup.dict_params["48I"]["fermion_params"][1][inv_acc]["mass"] = 0.03580

rup.dict_params["32Dfine"]["trajs"] = list(range(500, 3000, 10))
rup.dict_params["16IH2"]["trajs"] = list(range(500, 10000, 50))
rup.dict_params["32IfineH"]["trajs"] = list(range(500, 10000, 50))

rup.dict_params["48I"]["trajs"] = list(range(3000, 500, -5))
rup.dict_params["64I"]["trajs"] = list(range(3000, 500, -5))

qg.begin_with_gpt()

# ADJUST ME
job_tags = [
        "test-4nt8", "test-4nt16",
        # "64I",
        # "48I",
        # "24D",
        # "16IH2",
        # "32IfineH",
        ]

q.check_time_limit()

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        run_job(job_tag, traj)

qg.end_with_gpt()
