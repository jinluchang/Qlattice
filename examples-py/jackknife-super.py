#!/usr/bin/env python3

import qlat as q
import numpy as np
import functools

q.begin_with_mpi()

job_tag_list = [ 'test1', 'test2', ]

q.default_g_jk_kwargs["jk_type"] = "super"
@functools.lru_cache
def get_all_jk_idx():
    jk_idx_list = [ 'avg', ]
    for job_tag in job_tag_list:
        trajs = get_trajs(job_tag)
        for traj in trajs:
            jk_idx_list.append((job_tag, traj,))
    return jk_idx_list
q.default_g_jk_kwargs["get_all_jk_idx"] = get_all_jk_idx

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

for i in range(len(avg)):
    q.json_results_append(f"avg[{i}]", avg[i])
for i in range(len(avg)):
    q.json_results_append(f"err[{i}]", err[i])

q.check_log_json(__file__, check_eps=1e-10)

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
