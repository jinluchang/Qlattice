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
        os.path.join(os.getenv("HOME"), "qcddata"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/mk-gf-gt/results"),
        ]

@q.timer
def run_job(job_tag, traj):
    fns_produce = [
            f"eig/{job_tag}/traj={traj}",
            ]
    fns_need = [
            (f"configs/{job_tag}/ckpoint_lat.{traj}", f"configs/{job_tag}/ckpoint_lat.IEEE64BIG.{traj}"),
            ]
    if not check_job(job_tag, traj, fns_produce, fns_need):
        return
    #
    get_gf = run_gf(job_tag, traj)
    get_eig = run_eig(job_tag, traj, get_gf)
    #
    q.clean_cache()
    q.timer_display()

rup.dict_params["test-4nt8"]["trajs"] = list(range(1000, 1400, 100))
rup.dict_params["test-4nt16"]["trajs"] = list(range(1000, 1400, 100))
rup.dict_params["32Dfine"]["trajs"] = list(range(500, 3000, 10))
rup.dict_params["16IH2"]["trajs"] = list(range(500, 10000, 100))

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
        # "16IH2",
        # "24D",
        # "32Dfine",
        # "32IfineH",
        ]

q.check_time_limit()

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        run_job(job_tag, traj)

qg.end_with_gpt()
