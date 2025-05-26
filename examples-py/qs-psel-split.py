#!/usr/bin/env python3

import sys, os
import numpy as np
import qlat as q

import qlat_scripts.v1 as qs

from qlat_scripts.v1 import (
    set_param,
    get_param,
)

@q.timer(is_timer_fork=True)
def run_check_psel(get_psel):
    q.json_results_append(q.get_fname())
    psel = get_psel()
    q.json_results_append(f"len(psel) = {len(psel)}")
    psel_str = f"{psel.total_site} {psel[:].tolist()}"
    q.json_results_append(f"hash(psel)={q.hash_sha256(psel_str)}")

@q.timer(is_timer_fork=True)
def run_check_psel_list(get_psel_list):
    q.json_results_append(q.get_fname())
    psel_list = get_psel_list()
    q.json_results_append(f"len(psel_list) = {len(psel_list)}")
    for idx, psel in enumerate(psel_list):
        q.json_results_append(f"idx={idx} ; len(psel) = {len(psel)}")
        psel_str = f"{psel.total_site} {psel[:].tolist()}"
        q.json_results_append(f"idx={idx} ; hash(psel)={q.hash_sha256(psel_str)}")

@q.timer(is_timer_fork=True)
def run_job(job_tag, traj):
    get_f_weight = qs.run_f_weight_uniform(job_tag, traj)
    get_f_rand_01 = qs.run_f_rand_01(job_tag, traj)
    get_fsel_prob = qs.run_fsel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
    get_psel_prob = qs.run_psel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
    get_fsel = qs.run_fsel_from_fsel_prob(get_fsel_prob)
    get_psel = qs.run_psel_from_psel_prob(get_psel_prob)
    num_piece = 32
    get_psel_list = qs.run_psel_split(job_tag, traj, get_psel=get_psel, num_piece=num_piece)
    get_fsel_psel_list = qs.run_fsel_split(job_tag, traj, get_fsel=get_fsel, num_piece=num_piece)
    q.json_results_append("run_check_psel(get_psel)")
    run_check_psel(get_psel)
    q.json_results_append("run_check_psel_list(get_psel_list)")
    run_check_psel_list(get_psel_list)
    q.json_results_append("run_check_psel_list(get_fsel_psel_list)")
    run_check_psel_list(get_fsel_psel_list)

# --------------------------------------------

job_tag = "test-4nt8-checker"

set_param(job_tag, "traj_list")([ 1000, 1100, ])
set_param(job_tag, "total_site")([ 4, 4, 4, 8, ])
set_param(job_tag, "field_selection_fsel_rate")(0.1)
set_param(job_tag, "field_selection_psel_rate")(0.01)

job_tag = "test-8nt16-checker"

set_param(job_tag, "traj_list")([ 1000, 1100, ])
set_param(job_tag, "total_site")([ 8, 8, 8, 16, ])
set_param(job_tag, "field_selection_fsel_rate")(0.01)
set_param(job_tag, "field_selection_psel_rate")(0.001)

# --------------------------------------------

job_tag_list = q.get_arg("--job_tag_list", default="").split(",")

q.begin_with_mpi()

job_tag_list_default = [
        "test-4nt8-checker",
        "test-8nt16-checker",
        ]

if job_tag_list == [ "", ]:
    job_tag_list = job_tag_list_default

for job_tag in job_tag_list:
    for traj in get_param(job_tag, "traj_list"):
        run_job(job_tag, traj)

q.check_log_json(__file__, check_eps=1e-10)

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
