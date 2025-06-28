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
    q.json_results_append(f"psel.xg_arr {psel.xg_arr}")

@q.timer(is_timer_fork=True)
def run_check_fsel(get_fsel):
    q.json_results_append(q.get_fname())
    fsel = get_fsel()
    q.json_results_append(f"len(fsel) = {len(fsel)}")
    sig = q.get_data_sig(fsel.to_psel().xg_arr, q.RngState("seed-sig"))
    q.json_results_append(f"fsel.to_psel().xg_arr sig", sig, 1e-14)

@q.timer(is_timer_fork=True)
def run_check_fsel_prob(get_fsel_prob):
    q.json_results_append(q.get_fname())
    fsel_prob = get_fsel_prob()
    psel = fsel_prob.fsel.to_psel()
    psel_prob = q.SelectedPointsRealD(psel, 1)
    psel_prob @= fsel_prob
    sig = q.get_data_sig(psel_prob[:], q.RngState("seed-sig"))
    q.json_results_append(f"psel_prob from fsel_prob sig", sig, 1e-14)

@q.timer(is_timer_fork=True)
def run_check_psel_prob(get_psel_prob):
    q.json_results_append(q.get_fname())
    psel_prob = get_psel_prob()
    sig = q.get_data_sig(psel_prob[:], q.RngState("seed-sig"))
    q.json_results_append(f"psel_prob sig", sig, 1e-14)

@q.timer(is_timer_fork=True)
def run_job(job_tag, traj):
    get_f_weight = qs.run_f_weight_uniform(job_tag, traj)
    get_f_rand_01 = qs.run_f_rand_01(job_tag, traj)
    get_fsel_prob = qs.run_fsel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
    get_psel_prob = qs.run_psel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
    get_fsel = qs.run_fsel_from_fsel_prob(get_fsel_prob)
    get_psel = qs.run_psel_from_psel_prob(get_psel_prob)
    get_fsel_prob_sub = qs.run_fsel_prob_sub_sampling(
        job_tag, traj,
        sub_sampling_rate=0.5,
        get_fsel_prob=get_fsel_prob,
        get_f_rand_01=get_f_rand_01,
        get_f_weight=get_f_weight,
    )
    get_psel_prob_sub = qs.run_psel_prob_sub_sampling(
        job_tag, traj,
        sub_sampling_rate=0.5,
        get_psel_prob=get_psel_prob,
        get_f_rand_01=get_f_rand_01,
        get_f_weight=get_f_weight,
    )
    get_fsel_sub = qs.run_fsel_from_fsel_prob(get_fsel_prob_sub)
    get_psel_sub = qs.run_psel_from_psel_prob(get_psel_prob_sub)
    q.json_results_append("run_check_psel(get_psel)")
    run_check_psel(get_psel)
    q.json_results_append("run_check_psel(get_psel_sub)")
    run_check_psel(get_psel_sub)
    q.json_results_append("run_check_fsel(get_fsel)")
    run_check_fsel(get_fsel)
    q.json_results_append("run_check_fsel(get_fsel_sub)")
    run_check_fsel(get_fsel_sub)
    q.json_results_append("run_check_psel_prob(get_psel_prob)")
    run_check_psel_prob(get_psel_prob)
    q.json_results_append("run_check_psel_prob(get_psel_prob_sub)")
    run_check_psel_prob(get_psel_prob_sub)
    q.json_results_append("run_check_fsel_prob(get_fsel_prob)")
    run_check_fsel_prob(get_fsel_prob)
    q.json_results_append("run_check_fsel_prob(get_fsel_prob_sub)")
    run_check_fsel_prob(get_fsel_prob_sub)

# --------------------------------------------

job_tag = "test-4nt8-checker"

set_param(job_tag, "traj_list")([ 1000, 1100, ])
set_param(job_tag, "total_site")([ 4, 4, 4, 8, ])
set_param(job_tag, "field_selection_fsel_rate")(0.1)
set_param(job_tag, "field_selection_psel_rate")(0.01)

# --------------------------------------------

q.begin_with_mpi()

for traj in get_param(job_tag, "traj_list"):
    run_job(job_tag, traj)

q.check_log_json(__file__, check_eps=1e-10)

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
