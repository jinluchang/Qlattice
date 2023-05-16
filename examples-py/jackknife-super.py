#!/usr/bin/env python3

import qlat as q
import numpy as np
import functools

q.begin_with_mpi()

job_tags = [ 'test1', 'test2', ]

q.default_g_jk_kwargs["jk_type"] = "super"
@functools.cache
def get_all_jk_idx():
    jk_idx_list = [ 'avg', ]
    for job_tag in job_tags:
        trajs = get_trajs(job_tag)
        for traj in trajs:
            jk_idx_list.append((job_tag, traj,))
    return jk_idx_list
q.default_g_jk_kwargs["get_all_jk_idx"] = get_all_jk_idx

@functools.cache
def get_trajs(job_tag):
    return list(range(25))

rs = q.RngState("seed")
job_tag = "test1"
trajs = get_trajs(job_tag)

data_list = np.zeros((len(trajs), 5,)) # can be list or np.array
rs.g_rand_fill(data_list)
jk_list = q.g_jk(data_list)
jk_idx_list = [ "avg", ] + [ (job_tag, traj) for traj in trajs ]
jk_list = q.g_rejk(jk_list, jk_idx_list)
avg, err = q.g_jk_avg_err(jk_list)

q.displayln_info(f"CHECK: {avg}")
q.displayln_info(f"CHECK: {err}")

q.displayln_info(f"CHECK: finished successfully.")

q.end_with_mpi()
