#!/usr/bin/env python3

# Need --mpi X.X.X.X runtime option

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
        ]

@q.timer
def check_job(job_tag, traj):
    # return True if config is finished or unavailable
    fns_produce = [
            get_load_path(f"eig/{job_tag}/traj={traj}"),
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

@q.timer
def run_job(job_tag, traj):
    if check_job(job_tag, traj):
        return
    #
    get_gf = run_gf(job_tag, traj)
    get_eig = run_eig(job_tag, traj, get_gf)
    #
    q.clean_cache()
    q.timer_display()

rup.dict_params["test-4nt8"]["trajs"] = list(range(1000, 1400, 100))
rup.dict_params["test-4nt16"]["trajs"] = list(range(1000, 1400, 100))

qg.begin_with_gpt()

# ADJUST ME
job_tags = [
        "test-4nt8", "test-4nt16",
        # "test-8nt16",
        # "test-16nt32",
        # "test-32nt64",
        # "test-48nt96",
        # "test-64nt128",
        # "test-96nt192",
        # "test-128nt256",
        # "24D",
        ]

q.check_time_limit()

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        run_job(job_tag, traj)

qg.end_with_gpt()
