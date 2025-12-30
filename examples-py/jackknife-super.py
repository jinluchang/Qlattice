#!/usr/bin/env python3

import qlat as q
import numpy as np
import functools

q.begin_with_mpi()

job_tag_list = ['job_tag_1', 'job_tag_2', ]

@functools.lru_cache
def get_traj_list(job_tag):
    fname = q.get_fname()
    if job_tag == "job_tag_1":
        return list(range(20))
    elif job_tag == "job_tag_2":
        return list(range(30))
    else:
        raise Exception(f"{fname}: job_tag='{job_tag}'")
    return None

@functools.lru_cache
def get_all_jk_idx():
    jk_idx_list = ['avg', ]
    for job_tag in job_tag_list:
        traj_list = get_traj_list(job_tag)
        for traj in traj_list:
            jk_idx_list.append((job_tag, traj,))
    return jk_idx_list


q.default_g_jk_kwargs["jk_type"] = "super"
q.default_g_jk_kwargs["get_all_jk_idx"] = get_all_jk_idx

rs = q.RngState("seed1")
job_tag = "job_tag_1"
traj_list = get_traj_list(job_tag)

data_arr = rs.g_rand_arr((len(traj_list), 5,))  # can be list or np.array
jk_arr_1 = q.g_mk_jk(data_arr, [(job_tag, traj) for traj in traj_list])
avg, err = q.g_jk_avg_err(jk_arr_1)

for i in range(len(avg)):
    q.json_results_append(f"avg[{i}]", avg[i])
    q.json_results_append(f"err[{i}]", err[i])

rs = q.RngState("seed2")
job_tag = "job_tag_2"
traj_list = get_traj_list(job_tag)

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
