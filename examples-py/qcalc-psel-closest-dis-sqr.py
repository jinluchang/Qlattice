#!/usr/bin/env python3

import sys, os
import numpy as np
import glob
import pickle
import qlat as q

from qlat_scripts.v1 import (
    set_param,
    get_param,
    is_test,
    run_params,
    run_f_weight_uniform,
    run_f_rand_01,
    run_fsel_prob,
    run_psel_prob,
    run_fsel_from_fsel_prob,
    run_psel_from_psel_prob,
    run_psel_split,
    run_fsel_split,
    get_load_path,
)

usage = f"""
{__file__} --test
# Generate some test data and then perform the conversion.
{__file__} --usage
# Show this message.
{""}
{__file__} --src PATH_SRC_1 PATH_SRC_2 ...
# Show information for a list of `psel` provided by `PATH_SRC_1`, `PATH_SRC_2`, ...
# E.g.: {__file__} --src results/test-4nt8/points-selection/traj-1000.lati results/test-4nt8/points-selection/traj-2000.lati
"""

@q.timer(is_timer_fork=True)
def gen_test_data():
    job_tag_list = [
            "test-4nt8-checker",
            "test-8nt16-checker",
            ]
    job_tag_traj_list = []
    for job_tag in job_tag_list:
        run_params(job_tag)
        traj_list = get_param(job_tag, "traj_list")
        for traj in traj_list:
            job_tag_traj_list.append((job_tag, traj,))
    fn_list = []
    for job_tag, traj in job_tag_traj_list:
        get_f_weight = run_f_weight_uniform(job_tag, traj)
        get_f_rand_01 = run_f_rand_01(job_tag, traj)
        get_fsel_prob = run_fsel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
        get_psel_prob = run_psel_prob(job_tag, traj, get_f_rand_01=get_f_rand_01, get_f_weight=get_f_weight)
        get_fsel = run_fsel_from_fsel_prob(get_fsel_prob)
        get_psel = run_psel_from_psel_prob(get_psel_prob)
        num_piece = 16
        get_psel_list = run_psel_split(job_tag, traj, get_psel=get_psel, num_piece=num_piece)
        get_fsel_psel_list = run_fsel_split(job_tag, traj, get_fsel=get_fsel, num_piece=num_piece)
        path_psel_list = f"{job_tag}/points-selection-split/traj-{traj}/num-piece-{num_piece}"
        for idx in range(num_piece):
            fn = get_load_path(f"{path_psel_list}/idx-piece-{idx}.lati")
            fn_list.append(fn)
        path_psel_list = f"{job_tag}/field-selection-split/traj-{traj}/num-piece-{num_piece}"
        for idx in range(num_piece):
            fn = get_load_path(f"{path_psel_list}/idx-piece-{idx}.lati")
            fn_list.append(fn)
    return fn_list

@q.timer(is_timer_fork=True)
def run_psel_closest_dis_sqr(fn_list):
    fn_list_init = fn_list
    fn_list = []
    psel_list = []
    for fn in fn_list_init:
        if not fn.endswith(".lati"):
            q.displayln_info(-1, f"WARNING: '{fn}' does not endswith '.lati'. Skip this file.")
            continue
        if not q.does_file_exist_qar_sync_node(fn):
            q.displayln_info(-1, f"WARNING: '{fn}' does not exist. Skip this file.")
            continue
        psel = q.PointsSelection()
        psel.load(fn)
        fn_list.append(fn)
        psel_list.append(psel)
    closest_dis_sqr_list = q.find_closest_dis_sqr_for_psel_list(psel_list)
    assert len(fn_list) == len(psel_list)
    assert len(fn_list) == len(closest_dis_sqr_list)
    ret = (fn_list, psel_list, closest_dis_sqr_list,)
    q.displayln_info(f"INFO: idx closest_dis_sqr_list[idx] len(psel_list[idx]) psel_list[idx].total_site fn_list[idx]")
    for idx in range(len(fn_list)):
        info_str = f"RESULT: {idx:10} {closest_dis_sqr_list[idx]:10} {len(psel_list[idx]):10} {str(psel_list[idx].total_site.to_list()):>20} {fn_list[idx]}"
        q.displayln_info(info_str)
        if is_test():
            check_str = f"closest_dis_sqr_list[{idx}]={closest_dis_sqr_list[idx]} len(psel_list[{idx}])={len(psel_list[idx])} psel_list[{idx}].total_site={psel_list[idx].total_site} fn_list[idx]='{fn_list[idx]}'"
            q.json_results_append(check_str)
    return ret

@q.timer(is_timer_fork=True)
def run():
    if is_test():
        q.displayln_info(f"Usage:{usage}")
        q.displayln_info(f"Will now generate test data and run conversion.")
        fn_list = gen_test_data()
    else:
        fn_list = q.get_all_arg_list("--src")
        if fn_list is None:
            q.displayln_info(f"Usage:{usage}")
            return
    run_psel_closest_dis_sqr(fn_list)

# --------------------------------------------

job_tag = "test-4nt8-checker"
set_param(job_tag, "traj_list")([ 1000, 1100, ])
set_param(job_tag, "total_site")([ 4, 4, 4, 8, ])
set_param(job_tag, "field_selection_fsel_rate")(1)
set_param(job_tag, "field_selection_psel_rate")(1 / 4)

job_tag = "test-8nt16-checker"
set_param(job_tag, "traj_list")([ 1000, 1100, ])
set_param(job_tag, "total_site")([ 8, 8, 8, 16, ])
set_param(job_tag, "field_selection_fsel_rate")(1 / 32)
set_param(job_tag, "field_selection_psel_rate")(1 / 128)

job_tag = "test-16nt32-checker"
set_param(job_tag, "traj_list")([ 1000, 1100, ])
set_param(job_tag, "total_site")([ 16, 16, 16, 32, ])
set_param(job_tag, "field_selection_fsel_rate")(1 / 32)
set_param(job_tag, "field_selection_psel_rate")(1 / 128)

job_tag = "test-48nt96-checker"
set_param(job_tag, "traj_list")([ 1000, ])
set_param(job_tag, "total_site")([ 48, 48, 48, 96, ])
set_param(job_tag, "field_selection_fsel_rate")(4096 / (48**3 * 96))
set_param(job_tag, "field_selection_psel_rate")(2048 / (48**3 * 96))

job_tag = "test-64nt128-checker"
set_param(job_tag, "traj_list")([ 1000, ])
set_param(job_tag, "total_site")([ 64, 64, 64, 64, ])
set_param(job_tag, "field_selection_fsel_rate")(4096 / (64**3 * 128))
set_param(job_tag, "field_selection_psel_rate")(2048 / (64**3 * 128))

# --------------------------------------------

if __name__ == "__main__":
    is_show_usage = q.get_option("--usage")
    if is_show_usage:
        q.displayln_info(f"Usage:{usage}")
        exit()
    q.begin_with_mpi()
    run()
    q.timer_display()
    if is_test():
        q.check_log_json(__file__, check_eps=1e-10)
    q.end_with_mpi()
    q.displayln_info(f"CHECK: finished successfully.")
