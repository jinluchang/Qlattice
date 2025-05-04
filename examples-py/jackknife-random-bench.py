#!/usr/bin/env python3

import qlat as q
import numpy as np
import functools

q.begin_with_mpi()

@q.timer
def bench_rjk(n_data_sample, n_rand_sample, is_normalizing_rand_sample):
    fname = q.get_fname()
    q.default_g_jk_kwargs["jk_type"] = "rjk"
    q.default_g_jk_kwargs["n_rand_sample"] = n_rand_sample
    q.default_g_jk_kwargs["rng_state"] = q.RngState("rejk")
    q.default_g_jk_kwargs["jk_blocking_func"] = None
    q.default_g_jk_kwargs["is_normalizing_rand_sample"] = is_normalizing_rand_sample
    @functools.lru_cache
    def get_trajs(job_tag):
        return list(range(n_data_sample))
    rs = q.RngState("seed")
    job_tag = "test1"
    trajs = get_trajs(job_tag)
    data_arr = rs.g_rand_arr((len(trajs), 16,)) # can be list or np.array
    jk_arr = q.g_jk(data_arr)
    jk_idx_list = [ "avg", ] + [ (job_tag, traj) for traj in trajs ]
    jk_arr = q.g_rejk(jk_arr, jk_idx_list)
    avg, err = q.g_jk_avg_err(jk_arr)
    q.displayln_info(f"CHECK: n_data_sample={n_data_sample}")
    q.displayln_info(f"CHECK: n_rand_sample={n_rand_sample}")
    q.displayln_info(f"CHECK: is_normalizing_rand_sample={is_normalizing_rand_sample}")
    q.displayln_info(f"CHECK: ", " ".join([ f"{x:.4f}" for x in avg]))
    q.displayln_info(f"CHECK: ", " ".join([ f"{x:.4f}" for x in err]))
    q.json_results_append(f"{fname}: {n_data_sample} {n_rand_sample} {is_normalizing_rand_sample} avg", q.get_data_sig(avg, q.RngState()))
    for i in range(len(avg)):
        q.json_results_append(f"avg[{i}]", avg[i])
    q.json_results_append(f"{fname}: {n_data_sample} {n_rand_sample} {is_normalizing_rand_sample} err", q.get_data_sig(err, q.RngState()))
    for i in range(len(avg)):
        q.json_results_append(f"err[{i}]", err[i])

for n_data_sample in [ 1024, 128, 32, 8, 2, ]:
    for n_rand_sample in [ 1024, 128, 32, 8, 2, ]:
        for is_normalizing_rand_sample in [ True, False, ]:
            bench_rjk(n_data_sample, n_rand_sample, is_normalizing_rand_sample)

q.check_log_json(__file__, check_eps=1e-10)

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
