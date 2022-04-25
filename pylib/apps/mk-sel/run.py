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
            f"point-selection-smear/{job_tag}/traj={traj}.txt",
            f"wall-src-info-light/{job_tag}/traj={traj}.txt",
            f"wall-src-info-strange/{job_tag}/traj={traj}.txt",
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
    get_wi = run_wi(job_tag, traj)
    assert get_wi is not None
    #
    get_psel_smear = run_psel_smear(job_tag, traj)
    assert get_psel_smear is not None
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
rup.dict_params["32IcoarseH1"][tag] = list(range(300, 2000, 10))
rup.dict_params["24IH1"][tag] = list(range(1000, 10000, 10))
rup.dict_params["24IH2"][tag] = list(range(1000, 10000, 10))
rup.dict_params["24IH3"][tag] = list(range(1000, 10000, 10))
rup.dict_params["24IH4"][tag] = list(range(1000, 10000, 10))
rup.dict_params["32IH1"][tag] = list(range(1000, 10000, 10))
rup.dict_params["32IH2"][tag] = list(range(1000, 10000, 10)) + list(range(1002, 10000, 10))
rup.dict_params["32IH3"][tag] = list(range(1000, 10000, 10))

tag = "n_points_psel"
rup.dict_params["test-4nt8"][tag] = 6
rup.dict_params["test-4nt16"][tag] = 32
rup.dict_params["48I"][tag] = 2048
rup.dict_params["24D"][tag] = 1024
rup.dict_params["24DH"][tag] = 1024
rup.dict_params["32IfineH"][tag] = 512
rup.dict_params["32IcoarseH1"][tag] = 512
rup.dict_params["16IH2"][tag] = 256
rup.dict_params["24IH3"][tag] = 512
rup.dict_params["24IH2"][tag] = 512
rup.dict_params["24IH1"][tag] = 512
rup.dict_params["32IH1"][tag] = 512
rup.dict_params["32IH2"][tag] = 512
rup.dict_params["32IH3"][tag] = 512

tag = "n_exact_wsrc"
rup.dict_params["test-4nt8"][tag] = 2
rup.dict_params["48I"][tag] = 2

tag = "prob_exact_wsrc"
rup.dict_params["test-4nt16"][tag] = 1/8
rup.dict_params["16IH2"][tag] = 1/16
rup.dict_params["24D"][tag] = 1/32
rup.dict_params["24DH"][tag] = 1/32
rup.dict_params["32IfineH"][tag] = 1/32
rup.dict_params["32IcoarseH1"][tag] = 1/32
rup.dict_params["24IH1"][tag] = 1/32
rup.dict_params["24IH2"][tag] = 1/32
rup.dict_params["24IH3"][tag] = 1/32
rup.dict_params["32IH1"][tag] = 1/32
rup.dict_params["32IH2"][tag] = 1/32
rup.dict_params["32IH3"][tag] = 1/32

tag = "n_per_tslice_smear"
rup.dict_params["test-4nt8"][tag] = 2
rup.dict_params["test-4nt16"][tag] = 2
rup.dict_params["24D"][tag] = 16
rup.dict_params["24DH"][tag] = 16
rup.dict_params["16IH2"][tag] = 8
rup.dict_params["32IfineH"][tag] = 8
rup.dict_params["32IcoarseH1"][tag] = 8
rup.dict_params["24IH1"][tag] = 8
rup.dict_params["24IH2"][tag] = 8
rup.dict_params["24IH3"][tag] = 8
rup.dict_params["32IH1"][tag] = 8
rup.dict_params["32IH2"][tag] = 8
rup.dict_params["32IH3"][tag] = 8

tag = "prob_acc_1_smear"
rup.dict_params["test-4nt8"][tag] = 1/4
rup.dict_params["test-4nt16"][tag] = 1/4
rup.dict_params["24D"][tag] = 1/32
rup.dict_params["24DH"][tag] = 1/32
rup.dict_params["16IH2"][tag] = 1/16
rup.dict_params["32IfineH"][tag] = 1/32
rup.dict_params["32IcoarseH1"][tag] = 1/32
rup.dict_params["24IH1"][tag] = 1/32
rup.dict_params["24IH2"][tag] = 1/32
rup.dict_params["32IH2"][tag] = 1/32
rup.dict_params["32IH1"][tag] = 1/32

tag = "prob_acc_2_smear"
rup.dict_params["test-4nt8"][tag] = 1/16
rup.dict_params["test-4nt16"][tag] = 1/16
rup.dict_params["24D"][tag] = 1/128
rup.dict_params["24DH"][tag] = 1/128
rup.dict_params["16IH2"][tag] = 1/64
rup.dict_params["32IfineH"][tag] = 1/128
rup.dict_params["32IcoarseH1"][tag] = 1/128
rup.dict_params["24IH1"][tag] = 1/128
rup.dict_params["24IH2"][tag] = 1/128
rup.dict_params["32IH2"][tag] = 1/128
rup.dict_params["32IH1"][tag] = 1/128

tag = "prob_acc_1_psrc"
rup.dict_params["test-4nt8"][tag] = 1/4
rup.dict_params["test-4nt16"][tag] = 1/4
rup.dict_params["24D"][tag] = 1/32
rup.dict_params["24DH"][tag] = 1/32
rup.dict_params["16IH2"][tag] = 1/16
rup.dict_params["32IfineH"][tag] = 1/32
rup.dict_params["32IcoarseH1"][tag] = 1/32
rup.dict_params["24IH1"][tag] = 1/32
rup.dict_params["24IH2"][tag] = 1/32
rup.dict_params["24IH3"][tag] = 1/32
rup.dict_params["32IH1"][tag] = 1/32
rup.dict_params["32IH2"][tag] = 1/32

tag = "prob_acc_2_psrc"
rup.dict_params["test-4nt8"][tag] = 1/16
rup.dict_params["test-4nt16"][tag] = 1/16
rup.dict_params["24D"][tag] = 1/128
rup.dict_params["24DH"][tag] = 1/128
rup.dict_params["16IH2"][tag] = 1/64
rup.dict_params["32IfineH"][tag] = 1/128
rup.dict_params["32IcoarseH1"][tag] = 1/128
rup.dict_params["24IH1"][tag] = 1/128
rup.dict_params["24IH2"][tag] = 1/128
rup.dict_params["24IH3"][tag] = 1/128
rup.dict_params["32IH1"][tag] = 1/128
rup.dict_params["32IH2"][tag] = 1/128

tag = "n_rand_u1_fsel"
rup.dict_params["test-4nt8"][tag] = 4
rup.dict_params["test-4nt16"][tag] = 4
rup.dict_params["24D"][tag] = 64
rup.dict_params["24DH"][tag] = 64
rup.dict_params["48I"][tag] = 64
rup.dict_params["64I"][tag] = 64
rup.dict_params["16IH2"][tag] = 16
rup.dict_params["32IfineH"][tag] = 64
rup.dict_params["32IcoarseH1"][tag] = 64
rup.dict_params["24IH1"][tag] = 64
rup.dict_params["24IH2"][tag] = 64
rup.dict_params["24IH3"][tag] = 64
rup.dict_params["32IH1"][tag] = 64
rup.dict_params["32IH2"][tag] = 64

tag = "prob_acc_1_rand_u1"
rup.dict_params["test-4nt8"][tag] = 1/4
rup.dict_params["test-4nt16"][tag] = 1/4
rup.dict_params["24D"][tag] = 1/32
rup.dict_params["24DH"][tag] = 1/32
rup.dict_params["16IH2"][tag] = 1/16
rup.dict_params["32IfineH"][tag] = 1/32
rup.dict_params["32IcoarseH1"][tag] = 1/32
rup.dict_params["24IH1"][tag] = 1/32
rup.dict_params["24IH2"][tag] = 1/32
rup.dict_params["24IH3"][tag] = 1/32
rup.dict_params["32IH1"][tag] = 1/32
rup.dict_params["32IH2"][tag] = 1/32

tag = "prob_acc_2_rand_u1"
rup.dict_params["test-4nt8"][tag] = 1/16
rup.dict_params["test-4nt16"][tag] = 1/16
rup.dict_params["24D"][tag] = 1/128
rup.dict_params["24DH"][tag] = 1/128
rup.dict_params["16IH2"][tag] = 1/64
rup.dict_params["32IfineH"][tag] = 1/128
rup.dict_params["32IcoarseH1"][tag] = 1/128
rup.dict_params["24IH1"][tag] = 1/128
rup.dict_params["24IH2"][tag] = 1/128
rup.dict_params["24IH3"][tag] = 1/128
rup.dict_params["32IH1"][tag] = 1/128
rup.dict_params["32IH2"][tag] = 1/128

qg.begin_with_gpt()

# ADJUST ME
job_tags = [
        "test-4nt8", "test-4nt16",
        # "32IH1",
        # "32IH2",
        # "24IH1",
        # "24IH2",
        # "24IH3",
        # "16IH2",
        # "32IfineH",
        # "32IcoarseH1",
        # "24D",
        # "24DH",
        # "48I",
        ]

q.check_time_limit()

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        run_job(job_tag, traj)

qg.end_with_gpt()
