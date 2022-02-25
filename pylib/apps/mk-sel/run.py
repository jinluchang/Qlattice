#!/usr/bin/env python3

# Need --mpi X.X.X.X --mpi X.X.X runtime option

import qlat_gpt as qg
import rbc_ukqcd as ru

from jobs import *

load_path_list[:] = [
        "results",
        "../qcddata",
        os.path.join(os.getenv("HOME"), "qcddata"),
        "../mk-gf-gt/results",
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/mk-gf-gt/results"),
        ]

@q.timer_verbose
def run_job(job_tag, traj):
    fns_produce = [
            f"point-selection/{job_tag}/traj={traj}.txt",
            f"field-selection/{job_tag}/traj={traj}.field",
            ]
    fns_need = [
            (f"configs/{job_tag}/ckpoint_lat.{traj}", f"configs/{job_tag}/ckpoint_lat.IEEE64BIG.{traj}",),
            ]
    if not check_job(job_tag, traj, fns_produce, fns_need):
        return
    #
    get_psel = run_psel(job_tag, traj)
    get_fsel = run_fsel(job_tag, traj, get_psel)
    assert get_psel is not None
    assert get_fsel is not None
    #
    q.clean_cache()
    q.timer_display()

tag = "trajs"
rup.dict_params["test-4nt8"][tag] = list(range(1000, 1400, 100))
rup.dict_params["test-4nt16"][tag] = list(range(1000, 1400, 100))
rup.dict_params["48I"][tag] = list(range(3000, 500, -5))
rup.dict_params["24D"][tag] = list(range(1000, 10000, 10))
rup.dict_params["16IH2"][tag] = list(range(1000, 10000, 10))
rup.dict_params["32IfineH"][tag] = list(range(1000, 10000, 50))
rup.dict_params["24IH2"][tag] = list(range(1000, 10000, 10))
rup.dict_params["24IH1"][tag] = list(range(1000, 10000, 10))
rup.dict_params["32IH2"][tag] = list(range(1000, 10000, 10))

tag = "n_points_psel"
rup.dict_params["test-4nt8"][tag] = 6
rup.dict_params["test-4nt16"][tag] = 32
rup.dict_params["48I"][tag] = 2048
rup.dict_params["24D"][tag] = 1024
rup.dict_params["32IfineH"][tag] = 512
rup.dict_params["16IH2"][tag] = 256
rup.dict_params["24IH2"][tag] = 512
rup.dict_params["24IH1"][tag] = 512
rup.dict_params["32IH2"][tag] = 512

qg.begin_with_gpt()

# ADJUST ME
job_tags = [
        "test-4nt8", "test-4nt16",
        # "24IH1",
        # "24IH2",
        # "32IH2",
        # "16IH2",
        # "32IfineH",
        # "24D",
        # "48I",
        ]

q.check_time_limit()

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        run_job(job_tag, traj)

qg.end_with_gpt()
