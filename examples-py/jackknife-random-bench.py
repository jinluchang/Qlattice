#!/usr/bin/env python3

import qlat as q
import numpy as np


@q.timer
def bench_rjk(n_rand_sample, n_data_sample_1, n_data_sample_2):
    fname = q.get_fname()
    q.default_g_jk_kwargs["jk_type"] = "rjk"
    q.default_g_jk_kwargs["n_rand_sample"] = n_rand_sample
    q.default_g_jk_kwargs["rng_state"] = q.RngState("rejk")
    q.default_g_jk_kwargs["block_size"] = 1
    q.default_g_jk_kwargs["block_size_dict"] = {
        "job_tag_1": 1,
        "job_tag_2": 4,
    }
    #
    q.json_results_append(f"n_rand_sample={n_rand_sample}")
    q.json_results_append(f"n_data_sample_1={n_data_sample_1}")
    q.json_results_append(f"n_data_sample_2={n_data_sample_2}")
    #
    rs = q.RngState("seed1")
    job_tag = "job_tag_1"
    traj_list = list(range(n_data_sample_1))
    #
    data_arr = rs.g_rand_arr((len(traj_list), 5,))  # can be list or np.array
    jk_arr_1 = q.g_mk_jk(data_arr, [(job_tag, traj) for traj in traj_list])
    avg, err = q.g_jk_avg_err(jk_arr_1)
    #
    for i in range(len(avg)):
        q.json_results_append(f"avg[{i}]", avg[i])
        q.json_results_append(f"err[{i}]", err[i])
    #
    rs = q.RngState("seed2")
    job_tag = "job_tag_2"
    traj_list = list(range(n_data_sample_2))
    #
    data_arr = rs.g_rand_arr((len(traj_list), 5,))  # can be list or np.array
    jk_arr_2 = q.g_mk_jk(data_arr, [(job_tag, traj) for traj in traj_list])
    avg, err = q.g_jk_avg_err(jk_arr_2)
    #
    for i in range(len(avg)):
        q.json_results_append(f"avg[{i}]", avg[i])
        q.json_results_append(f"err[{i}]", err[i])
    #
    jk_arr = jk_arr_1 + jk_arr_2
    avg, err = q.g_jk_avg_err(jk_arr)
    #
    for i in range(len(avg)):
        q.json_results_append(f"avg[{i}]", avg[i])
        q.json_results_append(f"err[{i}]", err[i])


q.begin_with_mpi()

for n_rand_sample in [1024, 128, 32, 8, 2, ]:
    for n_data_sample_1, n_data_sample_2 in [
        (1024, 1024,),
        (256, 128,),
        (32, 64,),
        (8, 57,),
        (2, 2,),
    ]:
        bench_rjk(n_rand_sample, n_data_sample_1, n_data_sample_2)

q.check_log_json(__file__, check_eps=1e-10)
q.timer_display()
q.end_with_mpi()
q.displayln_info(f"CHECK: finished successfully.")
