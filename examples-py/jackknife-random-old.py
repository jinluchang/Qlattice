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
        "job_tag": 1,
        "test1": 1,
        }
q.default_g_jk_kwargs["is_normalizing_rand_sample"] = True
q.default_g_jk_kwargs["is_apply_rand_sample_jk_idx_blocking_shift"] = True

@functools.lru_cache
def get_traj_list(job_tag):
    return list(range(25))

rs = q.RngState("seed")
job_tag = "test1"
traj_list = get_traj_list(job_tag)

data_arr = rs.g_rand_arr((len(traj_list), 5,)) # can be list or np.array
jk_arr = q.g_jk(data_arr)
jk_idx_list = [ "avg", ] + [ (job_tag, traj) for traj in traj_list ]
jk_arr = q.g_rejk(jk_arr, jk_idx_list)
avg, err = q.g_jk_avg_err(jk_arr)

q.displayln_info(f"CHECK: {avg}")
q.displayln_info(f"CHECK: {err}")

for i in range(len(avg)):
    q.json_results_append(f"avg[{i}]", avg[i])
for i in range(len(avg)):
    q.json_results_append(f"err[{i}]", err[i])

q.check_log_json(__file__, check_eps=1e-10)

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
