#!/usr/bin/env python3

import argparse
import sys
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

default_params = q.default_run_flow_scale_params

usage = f"""
{__file__} --usage
# Show this message and exit.
{__file__} --test
# Generate some test data and then perform the scale flow, "Wilson flow" is the default. "Spatial flow" is performed if "--is-spatial True" option is specified.
{""}
{__file__} [ --step_size {default_params["step_size"]} ] [ --num_step {default_params["num_step"]} ] [ --is_spatial {default_params["is_spatial"]} ] [ --t_dir {default_params["t_dir"]} ] [ --integrator_type {default_params["integrator_type"]} ] --gf PATH_GAUGE_FIELD --out SCALE_FLOW_RECORD.pickle ...
...
# E.g.: {__file__} --gf qcddata-1/16IH2/configs/ckpoint_lat.IEEE64BIG.1000 --out results/16IH2/gf-flow-record/traj-1000.pickle
# E.g.: {__file__} --is_spatial True --gf qcddata-1/16IH2/configs/ckpoint_lat.IEEE64BIG.1000 --gf qcddata-1/16IH2/configs/ckpoint_lat.IEEE64BIG.1010 --out results/16IH2/gf-flow-record-spatial/traj-1000.pickle --out results/16IH2/gf-flow-record-spatial/traj-1010.pickle
# The program does not overwrite output files. If the output file already exist, the program will simply display that output file content and skip that input and output file pair and continue.
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
            job_tag_traj_list.append(
                (
                    job_tag,
                    traj,
                )
            )
    fn_gf_list = []
    fn_out_list = []
    for job_tag, traj in job_tag_traj_list:
        traj_gf = traj
        run_gf(job_tag, traj_gf)
        fn_gf = get_load_path(
            f"{job_tag}/configs/ckpoint_lat.{traj_gf}",
            f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj_gf}",
        )
        assert fn_gf is not None
        fn_out = get_save_path(f"{job_tag}/gf-flow-record/traj-{traj_gf}.pickle")
        fn_gf_list.append(fn_gf)
        fn_out_list.append(fn_out)
    num_job = len(job_tag_traj_list)
    argv_op_list = [
        [
            "--integrator_type",
            "euler",
        ],
        [
            "--num_step",
            "10",
            "--is_spatial",
            "True",
        ],
        [
            "--num_step",
            "10",
            "--is_spatial",
            "True",
            "--integrator_type",
            "euler",
        ],
        [
            "--num_step",
            "10",
            "--is_spatial",
            "True",
            "--integrator_type",
            "euler",
            "--t_dir",
            "0",
        ],
        [
            "--num_step",
            "10",
            "--step_size0.1--integrator_type",
            "euler",
        ],
    ]
    argv_list = []
    for id_job in range(num_job):
        if id_job >= len(argv_op_list):
            argv_op = []
        else:
            argv_op = argv_op_list[id_job]
        fn_gf = fn_gf_list[id_job]
        fn_out = fn_out_list[id_job]
        argv = (
            argv_op
            + [
                "--gf",
                fn_gf,
            ]
            + [
                "--out",
                fn_out,
            ]
        )
        argv_list.append(argv)
    return argv_list

def parse_params(argv):
    sys_args = parse_args(argv)
    params = dict()
    params["step_size"] = sys_args.step_size
    params["num_step"] = sys_args.num_step
    params["is_spatial"] = sys_args.is_spatial.lower() in ["true", "t", "yes", "y"]
    params["t_dir"] = sys_args.t_dir
    params["integrator_type"] = sys_args.integrator_type
    return params

@q.timer(is_timer_fork=True)
def run_job(argv):
    sys_args = parse_args(argv)
    fn_gf_list = sys_args.gf
    fn_out_list = sys_args.out
    assert len(fn_gf_list) == len(fn_out_list)
    params = parse_params(argv)
    for fn_gf, fn_out in zip(fn_gf_list, fn_out_list):
        q.run_flow_scale(fn_out, fn_gf=fn_gf, params=params)

@q.timer(is_timer_fork=True)
def run():
    if is_test():
        q.displayln_info("Will now generate test data and run instanton map.")
        argv_list = gen_test_data()
    else:
        argv_list = [
            sys.argv[1:],
        ]
    for argv in argv_list:
        run_job(argv)

def show_usage():
    q.displayln_info(f"Usage:{usage}")

class ExtendAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        items = getattr(namespace, self.dest, None) or []
        items.extend(values)
        setattr(namespace, self.dest, items)

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Scale flow (Wilson flow or Spatial flow)."
    )
    parser.add_argument(
        "--usage",
        action="store_true",
        default=False,
        help="Show usage",
    )
    parser.add_argument(
        "--step_size",
        type=float,
        default=default_params["step_size"],
        help=f"Step size (default: {default_params['step_size']})",
    )
    parser.add_argument(
        "--num_step",
        type=int,
        default=default_params["num_step"],
        help=f"Number of steps (default: {default_params['num_step']})",
    )
    parser.add_argument(
        "--is_spatial",
        type=str,
        default=f"{default_params['is_spatial']}",
        help=f"Whether to perform spatial flow (default: {default_params['is_spatial']})",
    )
    parser.add_argument(
        "--t_dir",
        type=int,
        default=default_params["t_dir"],
        help=f"Time direction (default: {default_params['t_dir']})",
    )
    parser.add_argument(
        "--integrator_type",
        type=str,
        default=default_params["integrator_type"],
        help=f"Integrator type (default: {default_params['integrator_type']})",
    )
    parser.add_argument(
        "--gf",
        action=ExtendAction,
        nargs="+",
        default=[],
        help="Gauge field paths",
    )
    parser.add_argument(
        "--out",
        action=ExtendAction,
        nargs="+",
        default=[],
        help="Output paths",
    )
    args, _ = parser.parse_known_args(args=argv)
    return args

# --------------------------------------------

job_tag = "test-4nt8-checker"
set_param(job_tag, "traj_list")(
    [
        1000,
        1100,
        1200,
        1300,
    ]
)
set_param(job_tag, "total_site")(
    [
        4,
        4,
        4,
        8,
    ]
)
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(20)
set_param(job_tag, "mk_sample_gauge_field", "rand_sigma")(0.25)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(4)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(1)
set_param(job_tag, "mk_sample_gauge_field", "hmc_beta")(3.0)

job_tag = "test-8nt16-checker"
set_param(job_tag, "traj_list")(
    [
        1000,
    ]
)
set_param(job_tag, "total_site")(
    [
        8,
        8,
        8,
        16,
    ]
)
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(10)
set_param(job_tag, "mk_sample_gauge_field", "rand_sigma")(0.25)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(4)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(1)
set_param(job_tag, "mk_sample_gauge_field", "hmc_beta")(5.0)

job_tag = "test-16nt32-checker"
set_param(job_tag, "traj_list")(
    [
        1000,
    ]
)
set_param(job_tag, "total_site")(
    [
        16,
        16,
        16,
        32,
    ]
)
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(5)
set_param(job_tag, "mk_sample_gauge_field", "rand_sigma")(0.25)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(8)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
set_param(job_tag, "mk_sample_gauge_field", "hmc_beta")(5.0)

job_tag = "test-48nt96-checker"
set_param(job_tag, "traj_list")(
    [
        1000,
    ]
)
set_param(job_tag, "total_site")(
    [
        48,
        48,
        48,
        96,
    ]
)
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(5)
set_param(job_tag, "mk_sample_gauge_field", "rand_sigma")(0.25)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(4)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
set_param(job_tag, "mk_sample_gauge_field", "hmc_beta")(5.0)

job_tag = "test-64nt128-checker"
set_param(job_tag, "traj_list")(
    [
        1000,
    ]
)
set_param(job_tag, "total_site")(
    [
        64,
        64,
        64,
        64,
    ]
)
set_param(job_tag, "mk_sample_gauge_field", "rand_n_step")(5)
set_param(job_tag, "mk_sample_gauge_field", "rand_sigma")(0.25)
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(4)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
set_param(job_tag, "mk_sample_gauge_field", "hmc_beta")(5.0)

# --------------------------------------------

if __name__ == "__main__":
    sys_args = parse_args()
    if sys_args.usage:
        show_usage()
        exit()
    q.begin_with_mpi()
    run()
    q.timer_display()
    if is_test():
        q.check_log_json(__file__, check_eps=1e-10)
    q.end_with_mpi()
    q.displayln_info("CHECK: finished successfully.")
