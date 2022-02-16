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
def compute_prop_rand_u1_type_acc(*, sfw, job_tag, traj, gf, eig, fsel, idx_rand_u1, inv_type, inv_acc, finished_tags):
    # same rand source for different inv_type
    tag = f"idx_rand_u1={idx_rand_u1} ; type={inv_type} ; accuracy={inv_acc}"
    if tag in finished_tags:
        return
    q.check_stop()
    q.check_time_limit()
    inv = ru.get_inv(gf, job_tag, inv_type, inv_acc, eig = eig)
    rs = q.RngState(f"seed {job_tag} {traj}").split(f"compute_prop_rand_u1(rand_u1)").split(idx_rand_u1)
    s_prop = q.mk_rand_u1_prop(inv, fsel, rs)
    s_prop.save_float_from_double(sfw, tag)
    sfw.flush()
    return s_prop

@q.timer_verbose
def compute_prop_rand_u1(*, job_tag, traj, inv_type, gf, path_s, fsel, eig = None):
    # use fsel instead of fselc
    n_rand_u1_fsel = rup.dict_params[job_tag]["n_rand_u1_fsel"]
    total_site = rup.dict_params[job_tag]["total_site"]
    geo = q.Geometry(total_site, 1)
    finished_tags = q.properly_truncate_fields(get_save_path(path_s + ".acc"))
    sfw = q.open_fields(get_save_path(path_s + ".acc"), "a", [ 1, 1, 1, 4, ])
    def comp(idx_rand_u1, inv_acc):
        compute_prop_rand_u1_type_acc(
                sfw = sfw,
                job_tag = job_tag, traj = traj,
                gf = gf, eig = eig, fsel = fsel,
                idx_rand_u1 = idx_rand_u1,
                inv_type = inv_type, inv_acc = inv_acc,
                finished_tags = finished_tags,)
    rs = q.RngState(f"seed {job_tag} {traj}").split(f"compute_prop_rand_u1(ama)")
    prob1 = rup.dict_params[job_tag]["prob_acc_1_rand_u1"]
    prob2 = rup.dict_params[job_tag]["prob_acc_2_rand_u1"]
    for idx_rand_u1 in range(n_rand_u1_fsel):
        r = rs.split(idx_rand_u1).u_rand_gen()
        inv_acc = 0
        assert 0 <= r and r <= 1
        comp(idx_rand_u1, inv_acc)
        inv_acc = 1
        if r <= prob1:
            comp(idx_rand_u1, inv_acc)
        inv_acc = 2
        if r <= prob2:
            comp(idx_rand_u1, inv_acc)
    q.clean_cache(q.cache_inv)
    sfw.close()
    q.qrename_info(get_save_path(path_s + ".acc"), get_save_path(path_s))

@q.timer_verbose
def run_prop_rand_u1(job_tag, traj, *, inv_type, get_gf, get_fsel, get_eig = None):
    if None in [ get_gf, get_fsel, ]:
        return
    if inv_type == 0 and get_eig is None:
        return
    inv_type_names = [ "light", "strange", "charm", ]
    inv_type_name = inv_type_names[inv_type]
    path_s = f"prop-rand-u1-{inv_type_name}/{job_tag}/traj={traj}"
    if get_load_path(path_s) is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-rand-u1-{inv_type_name}"):
        gf = get_gf()
        fsel, fselc = get_fsel()
        if get_eig is None:
            eig = None
        else:
            eig = get_eig()
        compute_prop_rand_u1(
                job_tag = job_tag, traj = traj,
                inv_type = inv_type,
                gf = gf,
                path_s = path_s,
                fsel = fsel,
                eig = eig)
        q.release_lock()

@q.timer_verbose
def run_job(job_tag, traj):
    fns_produce = [
            f"prop-rand-u1-light/{job_tag}/traj={traj}",
            f"prop-rand-u1-strange/{job_tag}/traj={traj}",
            f"prop-rand-u1-charm/{job_tag}/traj={traj}",
            ]
    fns_need = [
            (f"configs/{job_tag}/ckpoint_lat.{traj}", f"configs/{job_tag}/ckpoint_lat.IEEE64BIG.{traj}",),
            f"point-selection/{job_tag}/traj={traj}.txt",
            f"field-selection/{job_tag}/traj={traj}.field",
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
    get_fsel = run_fsel(job_tag, traj, get_psel)
    #
    def run_with_eig():
        get_eig = run_eig(job_tag, traj_gf, get_gf)
        run_prop_rand_u1(job_tag, traj, inv_type = 0, get_gf = get_gf, get_fsel = get_fsel, get_eig = get_eig)
    #
    run_with_eig()
    #
    def run_with_eig_strange():
        get_eig_strange = run_eig_strange(job_tag, traj_gf, get_gf)
        run_prop_rand_u1(job_tag, traj, inv_type = 1, get_gf = get_gf, get_fsel = get_fsel, get_eig = get_eig_strange)
    #
    run_with_eig_strange()
    #
    run_prop_rand_u1(job_tag, traj, inv_type = 2, get_gf = get_gf, get_fsel = get_fsel)
    #
    q.clean_cache()
    q.timer_display()

# rup.dict_params["test-4nt8"]["fermion_params"][0][2]["Ls"] = 10
# rup.dict_params["test-4nt8"]["fermion_params"][1][2]["Ls"] = 10

# rup.dict_params["test-4nt16"]["fermion_params"][0][2]["Ls"] = 10
# rup.dict_params["test-4nt16"]["fermion_params"][1][2]["Ls"] = 10

tag = "n_rand_u1_fsel"
rup.dict_params["test-4nt8"][tag] = 16
rup.dict_params["test-4nt16"][tag] = 16
rup.dict_params["48I"][tag] = 16
rup.dict_params["64I"][tag] = 16
rup.dict_params["16IH2"][tag] = 16
rup.dict_params["32IfineH"][tag] = 64

tag = "prob_acc_1_rand_u1"
rup.dict_params["test-4nt8"][tag] = 1/4
rup.dict_params["test-4nt16"][tag] = 1/4
rup.dict_params["16IH2"][tag] = 1/16
rup.dict_params["32IfineH"][tag] = 1/32

tag = "prob_acc_2_rand_u1"
rup.dict_params["test-4nt8"][tag] = 1/16
rup.dict_params["test-4nt16"][tag] = 1/16
rup.dict_params["16IH2"][tag] = 1/64
rup.dict_params["32IfineH"][tag] = 1/128

tag = "fermion_params"
for inv_acc in [ 0, 1, 2, ]:
    rup.dict_params["64I"][tag][0][inv_acc]["mass"] = 0.0006203
    rup.dict_params["64I"][tag][1][inv_acc]["mass"] = 0.02539
for inv_acc in [ 0, 1, 2, ]:
    rup.dict_params["48I"][tag][0][inv_acc]["mass"] = 0.0006979
    rup.dict_params["48I"][tag][1][inv_acc]["mass"] = 0.03580

tag = "trajs"
rup.dict_params["test-4nt8"][tag] = list(range(1000, 1400, 100))
rup.dict_params["test-4nt16"][tag] = list(range(1000, 1400, 100))
rup.dict_params["32Dfine"][tag] = list(range(500, 3000, 10))
rup.dict_params["16IH2"][tag] = list(range(500, 10000, 50))
rup.dict_params["32IfineH"][tag] = list(range(500, 10000, 50))
rup.dict_params["48I"][tag] = list(range(3000, 500, -5))
rup.dict_params["64I"][tag] = list(range(3000, 500, -5))

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
