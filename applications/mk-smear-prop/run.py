#!/usr/bin/env python3

# Need --mpi X.X.X.X --mpi X.X.X runtime option

import qlat_gpt as qg
import qlat_scripts.v1.rbc_ukqcd as ru

from qlat_scripts.v1.jobs import *
from qlat_scripts.v1.gen_data import *

import qlat_scripts.v1.params

load_path_list[:] = [
        "results",
        "qcddata",
        "qcddata-1",
        "qcddata-2",
        "qcddata-3",
        "qcddata-4",
        "qcddata-5",
        "../qcddata",
        "../qcddata-1",
        "../qcddata-2",
        "../qcddata-3",
        "../qcddata-4",
        "../qcddata-5",
        "../mk-gf-gt/results",
        "../mk-sel/results",
        "../mk-lanc/results",
        os.path.join(os.getenv("HOME"), "qcddata"),
        os.path.join(os.getenv("HOME"), "qcddata-1"),
        os.path.join(os.getenv("HOME"), "qcddata-2"),
        os.path.join(os.getenv("HOME"), "qcddata-3"),
        os.path.join(os.getenv("HOME"), "qcddata-4"),
        os.path.join(os.getenv("HOME"), "qcddata-5"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-gf-gt/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-sel/results"),
        os.path.join(os.getenv("HOME"), "Qlat-sample-data/default/mk-lanc/results"),
        ]

@q.timer_verbose
def run_job(job_tag, traj):
    fns_produce = [
            (f"{job_tag}/prop-smear-light/traj-{traj}.qar", f"{job_tag}/prop-smear-light/traj-{traj}/geon-info.txt",),
            (f"{job_tag}/psel-prop-smear-light/traj-{traj}.qar", f"{job_tag}/psel-prop-smear-light/traj-{traj}/checkpoint.txt",),
            (f"{job_tag}/prop-smear-strange/traj-{traj}.qar", f"{job_tag}/prop-smear-strange/traj-{traj}/geon-info.txt",),
            (f"{job_tag}/psel-prop-smear-strange/traj-{traj}.qar", f"{job_tag}/psel-prop-smear-strange/traj-{traj}/checkpoint.txt",),
            ]
    fns_need = [
            (f"{job_tag}/configs/ckpoint_lat.{traj}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj}",),
            f"{job_tag}/gauge-transform/traj-{traj}.field",
            f"{job_tag}/point-selection/traj-{traj}.txt",
            f"{job_tag}/field-selection/traj-{traj}.field",
            f"{job_tag}/point-selection-smear/traj-{traj}.txt",
            f"{job_tag}/eig/traj-{traj}",
            f"{job_tag}/eig/traj-{traj}/metadata.txt",
            f"{job_tag}/eig/traj-{traj}/eigen-values.txt",
            f"{job_tag}/eig-strange/traj-{traj}",
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
    get_gt = run_gt(job_tag, traj_gf, get_gf)
    #
    get_psel = run_psel(job_tag, traj)
    get_fsel = run_fsel(job_tag, traj, get_psel)
    assert get_psel is not None
    assert get_fsel is not None
    #
    get_gf_ape = run_gf_ape(job_tag, get_gf)
    #
    get_psel_smear = run_psel_smear(job_tag, traj)
    assert get_psel_smear is not None
    #
    def run_with_eig():
        get_eig = run_eig(job_tag, traj_gf, get_gf)
        run_get_inverter(job_tag, traj, inv_type = 0, get_gf = get_gf, get_eig = get_eig)
        run_prop_smear(job_tag, traj, inv_type = 0, get_gf = get_gf, get_gf_ape = get_gf_ape, get_eig = get_eig, get_gt = get_gt, get_psel = get_psel, get_fsel = get_fsel, get_psel_smear = get_psel_smear)
    #
    def run_with_eig_strange():
        get_eig = run_eig_strange(job_tag, traj_gf, get_gf)
        # get_eig = None
        run_get_inverter(job_tag, traj, inv_type = 1, get_gf = get_gf, get_eig = get_eig)
        run_prop_smear(job_tag, traj, inv_type = 1, get_gf = get_gf, get_gf_ape = get_gf_ape, get_eig = get_eig, get_gt = get_gt, get_psel = get_psel, get_fsel = get_fsel, get_psel_smear = get_psel_smear)
    #
    run_with_eig()
    run_with_eig_strange()
    #
    q.clean_cache()
    q.timer_display()

qg.begin_with_gpt()

# ADJUST ME
job_tags = [
        "test-4nt8", "test-4nt16",
        # "32IH1",
        # "32IH2",
        # "24IH1",
        # "24IH2",
        # "24IH3",
        # "24D",
        # "32Dfine",
        # "24DH",
        # "16IH2",
        # "32IfineH",
        # "32IcoarseH1",
        # "48I",
        ]

q.check_time_limit()

for job_tag in job_tags:
    q.displayln_info(pprint.pformat(rup.dict_params[job_tag]))
    for traj in rup.dict_params[job_tag]["trajs"]:
        run_job(job_tag, traj)

qg.end_with_gpt()
