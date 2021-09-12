#!/usr/bin/env python3

# Need --mpi X.X.X.X --mpi X.X.X runtime option

import qlat as q
import gpt as g
import qlat_gpt as qg
import rbc_ukqcd as ru
import rbc_ukqcd_params as rup
import pprint

import os

def get_save_path(fn):
    return os.path.join("results", fn)

def get_load_path(fn):
    if fn is None:
        return None
    path_list = [
            "results",
            "/gpfs/alpine/lgt116/proj-shared/ljin",
            ]
    for path in path_list:
        p = os.path.join(path, fn)
        if q.does_file_exist_sync_node(p):
            return p
    return None

@q.timer_verbose
def check_job(job_tag, traj):
    # return True if config is finished or unavailable
    fns_produce = [
            get_load_path(f"configs/{job_tag}/ckpoint_lat.{traj}"),
            get_load_path(f"gauge-transform/{job_tag}/traj={traj}.field"),
            ]
    is_job_done = True
    for fn in fns_produce:
        if fn is None:
            q.displayln_info(f"check_job: {job_tag} {traj} to do as some file does not exist.")
            is_job_done = False
            break
    if is_job_done:
        return True
    #
    fns_need = [
            get_load_path(f"configs/{job_tag}/ckpoint_lat.{traj}"),
            ]
    if job_tag[:5] == "test-":
        fns_need = []
    for fn in fns_need:
        if fn is None:
            q.displayln_info(f"check_job: {job_tag} {traj} unavailable as {fn} does not exist.")
            return True
    #
    q.check_stop()
    q.check_time_limit()
    #
    q.qmkdir_info(f"locks")
    q.qmkdir_info(get_save_path(f""))
    #
    return False

@q.timer_verbose
def run_gf(job_tag, traj):
    path_gf = get_load_path(f"configs/{job_tag}/ckpoint_lat.{traj}")
    if path_gf is None:
        if job_tag[:5] == "test-":
            gf = ru.mk_sample_gauge_field(job_tag, f"{traj}")
            q.qmkdir_info(get_save_path(f"configs"))
            q.qmkdir_info(get_save_path(f"configs/{job_tag}"))
            path_gf = get_save_path(f"configs/{job_tag}/ckpoint_lat.{traj}")
            # gf.save(path_gf)
            qg.save_gauge_field(gf, path_gf)
        else:
            assert False
    get_gf = ru.load_config_lazy(job_tag, path_gf)
    return get_gf

@q.timer_verbose
def run_gt(job_tag, traj, get_gf):
    if None in [ get_gf, ]:
        return None
    path_gt = get_load_path(f"gauge-transform/{job_tag}/traj={traj}.field")
    if path_gt is None:
        if q.obtain_lock(f"locks/{job_tag}-{traj}-gauge_fix_coulomb"):
            gf = get_gf()
            q.qmkdir_info(get_save_path(f"gauge-transform"))
            q.qmkdir_info(get_save_path(f"gauge-transform/{job_tag}"))
            gt = qg.gauge_fix_coulomb(gf)
            gt.save_double(get_save_path(f"gauge-transform/{job_tag}/traj={traj}.field"))
            q.release_lock()
            return lambda : gt
        else:
            return None
    else:
        @q.timer_verbose
        def load_gt():
            gt = q.GaugeTransform()
            gt.load_double(path_gt)
            # ADJUST ME
            # qg.check_gauge_fix_coulomb(get_gf(), gt)
            #
            return gt
        get_gt = q.lazy_call(load_gt)
    return get_gt

@q.timer
def run_job(job_tag, traj):
    if check_job(job_tag, traj):
        return
    #
    get_gf = run_gf(job_tag, traj)
    get_gt = run_gt(job_tag, traj, get_gf)
    #
    q.clean_cache()
    q.timer_display()

rup.dict_params["test-4nt8"]["n_points"] = [
        [ 6, 2, 1, ],
        [ 3, 2, 1, ],
        ]

rup.dict_params["test-4nt16"]["n_points"] = [
        [ 32, 4, 2, ],
        [ 16, 4, 2, ],
        ]

rup.dict_params["test-4nt8"]["n_points"] = [
        [ 1, 1, 1, ],
        [ 1, 1, 1, ],
        ]

rup.dict_params["test-4nt16"]["n_points"] = [
        [ 1, 1, 1, ],
        [ 1, 1, 1, ],
        ]

rup.dict_params["48I"]["n_points"] = [
        [ 2048, 64, 16, ],
        [ 1024, 64, 16, ],
        ]

rup.dict_params["test-4nt8"]["trajs"] = list(range(1000, 1400, 100))
rup.dict_params["test-4nt16"]["trajs"] = list(range(1000, 1400, 100))
rup.dict_params["48I"]["trajs"] = list(range(500, 3000, 5))

rup.dict_params["test-4nt8"]["fermion_params"][0][2]["Ls"] = 10
rup.dict_params["test-4nt8"]["fermion_params"][1][2]["Ls"] = 10

rup.dict_params["test-4nt16"]["fermion_params"][0][2]["Ls"] = 10
rup.dict_params["test-4nt16"]["fermion_params"][1][2]["Ls"] = 10

qg.begin_with_gpt()

# ADJUST ME
job_tags = [
        "test-4nt8",
        "test-4nt16",
        # "test-8nt16",
        # "test-16nt32",
        # "test-32nt64",
        # "test-48nt96",
        # "test-64nt128",
        # "test-96nt192",
        # "test-128nt256",
        # "24D",
        # "48I",
        ]

q.check_time_limit()

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        run_job(job_tag, traj)

qg.end_with_gpt()
