#!/usr/bin/env python3

import qlat as q
import numpy as np
import functools

q.begin_with_mpi()

q.default_g_jk_kwargs["jk_type"] = "rjk"
q.default_g_jk_kwargs["n_rand_sample"] = 1024
q.default_g_jk_kwargs["rng_state"] = q.RngState("rejk")
q.default_g_jk_kwargs["jk_blocking_func"] = None
q.default_g_jk_kwargs["is_normalizing_rand_sample"] = True

@functools.lru_cache
def get_trajs(job_tag):
    return list(range(25))

rs = q.RngState("seed")
job_tag = "test1"
trajs = get_trajs(job_tag)

data_arr = rs.g_rand_arr((len(trajs), 5,)) # can be list or np.array
jk_arr = q.g_jk(data_arr)
jk_idx_list = [ "avg", ] + [ (job_tag, traj) for traj in trajs ]
jk_arr = q.g_rejk(jk_arr, jk_idx_list)
avg, err = q.g_jk_avg_err(jk_arr)

q.displayln_info(f"CHECK: {avg}")
q.displayln_info(f"CHECK: {err}")

json_results = []
check_eps = 1e-10

for i in range(len(avg)):
    json_results.append((f"avg[{i}]", avg[i],))
for i in range(len(avg)):
    json_results.append((f"err[{i}]", err[i],))

q.check_log_json(__file__, json_results)

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
