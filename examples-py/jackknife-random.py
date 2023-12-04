#!/usr/bin/env python3

import qlat as q
import numpy as np
import functools

q.begin_with_mpi()

q.default_g_jk_kwargs["jk_type"] = "rjk"
q.default_g_jk_kwargs["n_rand_sample"] = 1024
q.default_g_jk_kwargs["rng_state"] = q.RngState("rejk")

@functools.lru_cache
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

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
