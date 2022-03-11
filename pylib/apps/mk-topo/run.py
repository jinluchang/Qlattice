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
        "../qcddata",
        os.path.join(os.getenv("HOME"), "qcddata"),
        "../mk-gf-gt/results",
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/mk-gf-gt/results"),
        ]

@q.timer
def run_topo(job_tag, traj, get_gf):
    if get_gf is None:
        return
    fn = f"topo-measure/{job_tag}/traj={traj}.pickle"
    if get_load_path(fn) is not None:
        return
    if q.obtain_lock(f"locks/{job_tag}-{traj}-topo"):
        gf = get_gf()
        topo_list = q.smear_measure_topo(gf)
        q.save_pickle_obj(topo_list, get_save_path(fn))
        q.release_lock()

@q.timer_verbose
def run_job(job_tag, traj):
    fns_produce = [
            f"topo-measure/{job_tag}/traj={traj}.pickle",
            ]
    fns_need = [
            (f"configs/{job_tag}/ckpoint_lat.{traj}", f"configs/{job_tag}/ckpoint_lat.IEEE64BIG.{traj}",),
            ]
    if job_tag[:5] == "test-":
        fns_need = []
    if not check_job(job_tag, traj, fns_produce, fns_need):
        return
    #
    get_gf = run_gf(job_tag, traj)
    #
    run_topo(job_tag, traj, get_gf)
    #
    q.clean_cache()
    q.timer_display()

tag = "trajs"
rup.dict_params["test-4nt8"][tag] = list(range(1000, 1400, 100))
rup.dict_params["test-4nt16"][tag] = list(range(1000, 1400, 100))
rup.dict_params["48I"][tag] = list(range(1000, 3000, 5))
rup.dict_params["24D"][tag] = list(range(1000, 10000, 10))
rup.dict_params["24DH"][tag] = list(range(200, 1000, 10))
rup.dict_params["32Dfine"][tag] = list(range(1000, 10000, 10))
rup.dict_params["16IH2"][tag] = list(range(1000, 10000, 10))
rup.dict_params["32IfineH"][tag] = list(range(1000, 10000, 10))
rup.dict_params["24IH1"][tag] = list(range(1000, 10000, 10))
rup.dict_params["24IH2"][tag] = list(range(1000, 10000, 10))
rup.dict_params["24IH3"][tag] = list(range(1000, 10000, 10))
rup.dict_params["24IH4"][tag] = list(range(1000, 10000, 10))
rup.dict_params["32IH1"][tag] = list(range(1000, 10000, 10))
rup.dict_params["32IH2"][tag] = list(range(1000, 10000, 10)) + list(range(1002, 10000, 10))
rup.dict_params["32IH3"][tag] = list(range(1000, 10000, 10))

qg.begin_with_gpt()

# ADJUST ME
job_tags = [
        "test-4nt8", "test-4nt16",
        # "24IH1",
        # "24IH2",
        # "24IH3",
        # "24IH4",
        # "32IH1",
        # "32IH2",
        # "32IH3",
        # "32IfineH",
        # "16IH2",
        # "24D",
        # "24DH",
        # "32Dfine",
        # "48I",
        ]

q.check_time_limit()

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        run_job(job_tag, traj)

qg.end_with_gpt()
