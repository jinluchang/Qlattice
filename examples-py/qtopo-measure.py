#!/usr/bin/env python3

import sys
import qlat as q
import numpy as np

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
Topological charge measurement based on Wilson flow with Qlattice
by Luchang Jin
2026/03/04
{""}
{__file__} --usage
# Show this message and exit.
{__file__} --test
# Generate some test data and then perform the topo flow.
{""}
{__file__} --gf PATH_GAUGE_FIELD --out PATH_OUTPUT --step_size 0.05 --n_step 80 --flow_type Wilson
{__file__} \
    --step_size 0.05 --n_step 80 --flow_type Wilson --integrator_type euler \
    --step_size 0.1 --n_step 100 --flow_type Localize --integrator_type euler \
    --step_size 0.1 --n_step 100 --flow_type Localize --integrator_type euler \
    --step_size 0.1 --n_step 100 --flow_type Localize --integrator_type euler \
    --step_size 0.1 --n_step 100 --flow_type Localize --integrator_type euler \
    --gf PATH_GAUGE_FIELD \
    --out PATH_OUTPUT \
    --out_df PATH_DENSITY_FIELD_OUTPUT
# Multiple input and output paths can be specified (they will be processed the same way).
# Default PATH_DENSITY_FIELD_OUTPUT is PATH_OUTPUT
# The program does not overwrite output files. If the output file already exist, the program will simply display that output file content and skip that input and output file pair and continue.
"""

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 3],
    [1, 1, 1, 4],
    [1, 1, 1, 6],
    [1, 1, 1, 8],
    [1, 2, 2, 4],
    [2, 2, 2, 4],
    [2, 2, 2, 4],
]


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
    fn_out_df_list = []
    for job_tag, traj in job_tag_traj_list:
        traj_gf = traj
        run_gf(job_tag, traj_gf)
        fn_gf = get_load_path(
            f"{job_tag}/configs/ckpoint_lat.{traj_gf}",
            f"{job_tag}/configs/ckpoint_lat.IEEE64BIG.{traj_gf}",
        )
        assert fn_gf is not None
        fn_out = get_save_path(f"{job_tag}/flow-topo-info/traj-{traj_gf}")
        fn_out_df = get_save_path(f"{job_tag}/flow-density-field/traj-{traj_gf}")
        fn_gf_list.append(fn_gf)
        fn_out_list.append(fn_out)
        fn_out_df_list.append(fn_out_df)
    argv = []
    for fn in fn_gf_list:
        argv += [
            "--gf",
            fn,
        ]
    for fn in fn_out_list:
        argv += [
            "--out",
            fn,
        ]
    for fn in fn_out_df_list:
        argv += [
            "--out_df",
            fn,
        ]
    argv += [
        "--step_size",
        "0.05",
        "--n_step",
        "80",
        "--flow_type",
        "Wilson",
        "--integrator_type",
        "euler",
        "--step_size",
        "0.1",
        "--n_step",
        "100",
        "--flow_type",
        "Localize",
        "--integrator_type",
        "euler",
        "--step_size",
        "0.1",
        "--n_step",
        "100",
        "--flow_type",
        "Localize",
        "--integrator_type",
        "euler",
    ]
    return argv


@q.timer(is_timer_fork=True)
def run_topo_measure(fn_gf, fn_out, fn_out_df, smear_info_list):
    fname = q.get_fname()
    if q.does_file_exist_qar_sync_node(fn_out):
        q.displayln_info(
            -1, f"{fname}: WARNING: '{fn_out}' for '{fn_gf}' already exist. Skip."
        )
        return
    if q.does_file_exist_qar_sync_node(fn_out_df):
        q.displayln_info(
            -1, f"{fname}: WARNING: '{fn_out_df}' for '{fn_gf}' already exist. Skip."
        )
        return
    if not q.does_file_exist_qar_sync_node(fn_gf):
        q.displayln_info(
            -1, f"{fname}: WARNING: '{fn_gf}' does not exist. Skip this file."
        )
        return
    q.json_results_append(
        f"{fname}: Start compute topo_measure out='{fn_out}' (out_df='{fn_out_df}') for '{fn_gf}'"
    )
    gf = q.GaugeField()
    gf.load(fn_gf)
    (
        topo_list,
        energy_list,
    ) = q.smear_measure_topo(
        gf,
        smear_info_list=smear_info_list,
        info_path=fn_out,
        density_field_path=fn_out_df,
        energy_derivative_info=None,
        is_show_topo_terms=False,
    )
    for name in [
        "flow_time",
        "plaq",
        "plaq_min",
        "plaq_max",
        "abs_topo",
        "topo",
        "topo_tslice",
        "topo_clf",
        "plaq_action_density",
        "energy_density",
        "energy_deriv",
    ]:
        q.json_results_append(
            f"{fname}: {name} in topo_list",
            np.array([d[name] for d in topo_list]),
            1e-8,
        )
    for name in [
        "flow_time",
        "plaq",
        "plaq_min",
        "plaq_max",
        "energy_density",
        "energy_density_tslice",
    ]:
        q.json_results_append(
            f"{fname}: {name} in energy_list",
            q.get_data_sig_arr(
                np.array([d[name] for d in energy_list]),
                q.RngState(),
                3,
            ),
            1e-8,
        )
    q.json_results_append(
        f"{fname}: End compute topo_measure out='{fn_out}' (out_df='{fn_out_df}') for '{fn_gf}'"
    )


def show_usage():
    q.displayln_info(f"Usage:{usage}")


@q.timer(is_timer_fork=True)
def run():
    if is_test():
        q.displayln_info("Will now generate test data and run topo measure.")
        argv = gen_test_data()
    else:
        argv = sys.argv
    fn_gf_list = q.get_arg_list("--gf", argv=argv)
    fn_out_list = q.get_arg_list("--out", argv=argv)
    fn_out_df_list = q.get_arg_list("--out_df", argv=argv)
    assert len(fn_gf_list) == len(fn_out_list)
    if len(fn_out_df_list) == 0:
        fn_out_df_list = fn_out_list
    assert len(fn_out_df_list) == len(fn_out_list)
    p_step_size_list = q.get_arg_list("--step_size", argv=argv)
    p_n_step_list = q.get_arg_list("--n_step", argv=argv)
    flow_type_list = q.get_arg_list("--flow_type", argv=argv)
    integrator_type_list = q.get_arg_list("--integrator_type", argv=argv)
    assert len(p_step_size_list) == len(p_n_step_list)
    assert len(p_step_size_list) == len(flow_type_list)
    assert len(p_step_size_list) == len(integrator_type_list)
    smear_info_list = []
    for step_size, n_step, flow_type, integrator_type in zip(
        p_step_size_list, p_n_step_list, flow_type_list, integrator_type_list
    ):
        smear_info = [
            float(step_size),
            int(n_step),
            flow_type,
            integrator_type,
        ]
        smear_info_list.append(smear_info)
    for fn_gf, fn_out, fn_out_df in zip(fn_gf_list, fn_out_list, fn_out_df_list):
        run_topo_measure(
            fn_gf,
            fn_out,
            fn_out_df,
            smear_info_list=smear_info_list,
        )


# --------------------------------------------

job_tag = "test-4nt8-checker"
set_param(job_tag, "traj_list")(
    [
        1000,
        1100,
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
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
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
set_param(job_tag, "mk_sample_gauge_field", "flow_n_step")(8)
set_param(job_tag, "mk_sample_gauge_field", "hmc_n_traj")(5)
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
    is_show_usage = q.get_option("--usage")
    if is_show_usage:
        show_usage()
        exit()
    q.begin_with_mpi(size_node_list)
    run()
    q.timer_display()
    if is_test():
        q.check_log_json(__file__, check_eps=1e-10)
    q.end_with_mpi()
    q.displayln_info("CHECK: finished successfully.")
