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
    run_gf,
    get_load_path,
    get_save_path,
)

usage = f"""
{__file__} --usage
# Show this message and exit.
{__file__} --test
# Generate some test data and then perform the topo flow.
{""}
{__file__} --gf PATH_GAUGE_FIELD --out INST_MAP.pickle
{__file__} --gf PATH_GAUGE_FIELD_1 --gf PATH_GAUGE_FIELD_2 --out INST_MAP_1.pickle --out INST_MAP_2.pickle
...
# E.g.: {__file__} --gf qcddata-1/16IH2/configs/ckpoint_lat.IEEE64BIG.1000 --out results/16IH2/instanton-map/traj-1000.pickle
# E.g.: {__file__} --gf qcddata-1/16IH2/configs/ckpoint_lat.IEEE64BIG.1000 --gf qcddata-1/16IH2/configs/ckpoint_lat.IEEE64BIG.1010 --out results/16IH2/instanton-map/traj-1000.pickle --out results/16IH2/instanton-map/traj-1010.pickle
# The program does not overwrite output files. If the output file already exist, the program will simply display that output file content and skip that input and output file pair and continue.
{""}
{__file__} --info INST_MAP_1.pickle INST_MAP_2.pickle ...
# E.g.: {__file__} --info results/16IH2/instanton-map/traj-1000.pickle results/16IH2/instanton-map/traj-1010.pickle
# Display these INST_MAP file contents.
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
    fn_gf_list = []
    fn_out_list = []
    for job_tag, traj in job_tag_traj_list:
        traj_gf = traj
        get_gf = run_gf(job_tag, traj_gf)
        fn_gf = get_load_path(f"{job_tag}/configs/ckpoint_lat.{traj_gf}", f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj_gf}")
        assert fn_gf is not None
        fn_out = get_save_path(f"{job_tag}/instanton-map/traj-{traj_gf}.pickle")
        fn_gf_list.append(fn_gf)
        fn_out_list.append(fn_out)
    argv = []
    for fn in fn_gf_list:
        argv += [ "--gf", fn, ]
    for fn in fn_out_list:
        argv += [ "--out", fn, ]
    argv += [ "--info", ]
    for fn in fn_out_list:
        argv += [ fn, ]
    return argv

@q.timer(is_timer_fork=True)
def run_inst_map(fn_out, fn_gf=None):
    fname = q.get_fname()
    if not fn_out.endswith(".pickle"):
        q.displayln_info(-1, f"{fname}: WARNING: '{fn_out}' does not endswith '.pickle'. Skip this file.")
        return
    if fn_gf is None:
        if not q.does_file_exist_qar_sync_node(fn_out):
            return
        q.json_results_append(f"{fname}: Start show inst_map_obj fn='{fn_out}'")
        inst_map_obj = q.load_pickle_obj(fn_out)
        q.displayln_info_inst_map_obj(inst_map_obj)
        q.json_results_append(f"{fname}: End show inst_map_obj fn='{fn_out}'")
        return
    if q.does_file_exist_qar_sync_node(fn_out):
        q.displayln_info(-1, f"{fname}: WARNING: '{fn_out}' for '{fn_gf}' already exist. Simply display its contents.")
        q.json_results_append(f"{fname}: Start show inst_map_obj fn='{fn_out}'")
        inst_map_obj = q.load_pickle_obj(fn_out)
        q.displayln_info_inst_map_obj(inst_map_obj)
        q.json_results_append(f"{fname}: End show inst_map_obj fn='{fn_out}'")
        return
    if not q.does_file_exist_qar_sync_node(fn_gf):
        q.displayln_info(-1, f"{fname}: WARNING: '{fn_gf}' does not exist. Skip this file.")
        return
    q.json_results_append(f"{fname}: Start compute inst_map_obj fn='{fn_out}' for '{fn_gf}'")
    gf = q.GaugeField()
    gf.load(fn_gf)
    inst_map_obj = q.compute_inst_map(gf)
    q.save_pickle_obj(inst_map_obj, fn_out)
    q.json_results_append(f"{fname}: End compute inst_map_obj fn='{fn_out}' for '{fn_gf}'")
    q.json_results_append(f"{fname}: Start show inst_map_obj fn='{fn_out}'")
    inst_map_obj = q.load_pickle_obj(fn_out)
    q.displayln_info_inst_map_obj(inst_map_obj)
    q.json_results_append(f"{fname}: End show inst_map_obj fn='{fn_out}'")

def show_usage():
    q.displayln_info(f"Usage:{usage}")

@q.timer(is_timer_fork=True)
def run():
    if is_test():
        q.displayln_info(f"Will now generate test data and run instanton map.")
        argv = gen_test_data()
    else:
        argv = sys.argv
    fn_gf_list = q.get_arg_list("--gf", argv=argv)
    fn_out_list = q.get_arg_list("--out", argv=argv)
    fn_info_list = q.get_all_arg_list("--info", argv=argv)
    assert len(fn_gf_list) == len(fn_out_list)
    for fn_gf, fn_out in zip(fn_gf_list, fn_out_list):
        run_inst_map(fn_out, fn_gf)
    if fn_info_list is not None:
        for fn_info in fn_info_list:
            run_inst_map(fn_info)

# --------------------------------------------

job_tag = "test-4nt8-checker"
set_param(job_tag, "traj_list")([ 1000, 1100, ])
set_param(job_tag, "total_site")([ 4, 4, 4, 8, ])
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(20)
set_param(job_tag, "mk_sample_gauge_field", "rand_sigma")(0.25)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(4)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
set_param(job_tag, "mk_sample_gauge_field", "hmc_beta")(3.0)

job_tag = "test-8nt16-checker"
set_param(job_tag, "traj_list")([ 1000, ])
set_param(job_tag, "total_site")([ 8, 8, 8, 16, ])
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(10)
set_param(job_tag, "mk_sample_gauge_field", "rand_sigma")(0.25)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(8)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
set_param(job_tag, "mk_sample_gauge_field", "hmc_beta")(5.0)

job_tag = "test-16nt32-checker"
set_param(job_tag, "traj_list")([ 1000, ])
set_param(job_tag, "total_site")([ 16, 16, 16, 32, ])
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(5)
set_param(job_tag, "mk_sample_gauge_field", "rand_sigma")(0.25)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(8)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
set_param(job_tag, "mk_sample_gauge_field", "hmc_beta")(5.0)

job_tag = "test-48nt96-checker"
set_param(job_tag, "traj_list")([ 1000, ])
set_param(job_tag, "total_site")([ 48, 48, 48, 96, ])
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(5)
set_param(job_tag, "mk_sample_gauge_field", "rand_sigma")(0.25)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(4)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
set_param(job_tag, "mk_sample_gauge_field", "hmc_beta")(5.0)

job_tag = "test-64nt128-checker"
set_param(job_tag, "traj_list")([ 1000, ])
set_param(job_tag, "total_site")([ 64, 64, 64, 64, ])
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(5)
set_param(job_tag, "mk_sample_gauge_field", "rand_sigma")(0.25)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(4)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
set_param(job_tag, "mk_sample_gauge_field", "hmc_beta")(5.0)

# --------------------------------------------

if __name__ == "__main__":
    is_show_usage = q.get_option("--usage")
    if is_show_usage:
        show_usage()
        exit()
    q.begin_with_mpi()
    run()
    q.timer_display()
    if is_test():
        q.check_log_json(__file__, check_eps=1e-10)
    q.end_with_mpi()
    q.displayln_info(f"CHECK: finished successfully.")
