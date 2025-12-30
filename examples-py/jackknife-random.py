#!/usr/bin/env python3

import qlat as q
import numpy as np
import functools

q.begin_with_mpi()

q.default_g_jk_kwargs["jk_type"] = "rjk"
q.default_g_jk_kwargs["n_rand_sample"] = 1024
q.default_g_jk_kwargs["rng_state"] = q.RngState("rejk")
q.default_g_jk_kwargs["block_size"] = 1
q.default_g_jk_kwargs["block_size_dict"] = {
    "job_tag_1": 1,
    "job_tag_2": 4,
}

rs = q.RngState("seed1")
job_tag = "job_tag_1"
traj_list = list(range(20))

data_arr = rs.g_rand_arr((len(traj_list), 5,))  # can be list or np.array
jk_arr_1 = q.g_mk_jk(data_arr, [(job_tag, traj) for traj in traj_list])
avg, err = q.g_jk_avg_err(jk_arr_1)

for i in range(len(avg)):
    q.json_results_append(f"avg[{i}]", avg[i])
    q.json_results_append(f"err[{i}]", err[i])

rs = q.RngState("seed2")
job_tag = "job_tag_2"
traj_list = list(range(30))

data_arr = rs.g_rand_arr((len(traj_list), 5,))  # can be list or np.array
jk_arr_2 = q.g_mk_jk(data_arr, [(job_tag, traj) for traj in traj_list])
avg, err = q.g_jk_avg_err(jk_arr_2)

for i in range(len(avg)):
    q.json_results_append(f"avg[{i}]", avg[i])
    q.json_results_append(f"err[{i}]", err[i])

jk_arr = jk_arr_1 + jk_arr_2
avg, err = q.g_jk_avg_err(jk_arr)

for i in range(len(avg)):
    q.json_results_append(f"avg[{i}]", avg[i])
    q.json_results_append(f"err[{i}]", err[i])

q.check_log_json(__file__, check_eps=1e-10)
q.end_with_mpi()
q.displayln_info(f"CHECK: finished successfully.")
