#!/usr/bin/env python3

import qlat as q
import numpy as np


@q.timer
def compute_avg_err_direct(data_list, block_size, rng_state):
    q.reset_default_g_jk_kwargs()
    avg, err = q.avg_err(data_list, block_size=block_size)
    return avg, err


@q.timer
def compute_avg_err_jk(data_list, block_size, rng_state):
    q.reset_default_g_jk_kwargs()
    jk_data_list = q.jackknife(data_list.tolist())
    avg, err = q.jk_avg_err(jk_data_list, block_size=block_size)
    return avg, err


@q.timer
def compute_avg_err_sjk_hash(data_list, block_size, rng_state):
    q.reset_default_g_jk_kwargs()
    q.default_g_jk_kwargs["jk_type"] = "super"
    q.default_g_jk_kwargs["eps"] = 1
    q.default_g_jk_kwargs["is_hash_jk_idx"] = True
    q.default_g_jk_kwargs["jk_idx_hash_size"] = 32
    q.default_g_jk_kwargs["block_size"] = 1
    q.default_g_jk_kwargs["block_size_dict"] = {
        "job_tag": block_size,
    }
    q.default_g_jk_kwargs["rng_state"] = rng_state
    q.default_g_jk_kwargs["all_jk_idx"] = None
    q.default_g_jk_kwargs["get_all_jk_idx"] = None
    q.default_g_jk_kwargs["all_jk_idx_set"] = set()
    n = len(data_list)
    jk_idx_list = [("job_tag", traj,) for traj in range(n)]
    jk_data_arr = q.g_mk_jk(data_list, jk_idx_list)
    avg, err = q.g_jk_avg_err(jk_data_arr)
    return avg, err


@q.timer
def compute_avg_err_sjk(data_list, block_size, rng_state):
    q.reset_default_g_jk_kwargs()
    q.default_g_jk_kwargs["jk_type"] = "super"
    q.default_g_jk_kwargs["eps"] = 1
    q.default_g_jk_kwargs["is_hash_jk_idx"] = True
    q.default_g_jk_kwargs["jk_idx_hash_size"] = 32
    q.default_g_jk_kwargs["block_size"] = 1
    q.default_g_jk_kwargs["block_size_dict"] = {
        "job_tag": block_size,
    }
    q.default_g_jk_kwargs["rng_state"] = rng_state
    q.default_g_jk_kwargs["all_jk_idx"] = None
    q.default_g_jk_kwargs["get_all_jk_idx"] = None
    q.default_g_jk_kwargs["all_jk_idx_set"] = set()
    n = len(data_list)
    jk_idx_list = [("job_tag", traj,) for traj in range(n)]
    q.default_g_jk_kwargs["all_jk_idx"] = (
        ["avg", ]
        + list(set(
            q.g_jk_blocking_func(0, jk_idx)
            for jk_idx in jk_idx_list
        ))
    )
    jk_data_arr = q.g_mk_jk(data_list, jk_idx_list)
    avg, err = q.g_jk_avg_err(jk_data_arr)
    return avg, err


@q.timer
def compute_avg_err_rjk(data_list, block_size, rng_state):
    q.reset_default_g_jk_kwargs()
    q.default_g_jk_kwargs["jk_type"] = "rjk"
    q.default_g_jk_kwargs["eps"] = 1
    q.default_g_jk_kwargs["n_rand_sample"] = 16
    q.default_g_jk_kwargs["is_normalizing_rand_sample"] = False
    q.default_g_jk_kwargs["is_apply_rand_sample_jk_idx_blocking_shift"] = True
    q.default_g_jk_kwargs["block_size"] = 1
    q.default_g_jk_kwargs["block_size_dict"] = {
        "job_tag": block_size,
    }
    q.default_g_jk_kwargs["rng_state"] = rng_state
    q.default_g_jk_kwargs["all_jk_idx_set"] = set()
    n = len(data_list)
    jk_idx_list = [("job_tag", traj,) for traj in range(n)]
    jk_data_arr = q.g_mk_jk(data_list, jk_idx_list)
    avg, err = q.g_jk_avg_err(jk_data_arr)
    return avg, err


@q.timer
def jk_test(data_list_size, block_size, compute_avg_err):
    #
    size = 128
    #
    avg_list = []
    err_list = []
    #
    for i in range(size):
        #
        data_list = q.RngState(f"seed-{i}").g_rand_arr(data_list_size)
        #
        avg, err = compute_avg_err(
            data_list, block_size, q.RngState(f"jk-seed-{i}"))
        #
        avg_list.append(avg)
        err_list.append(err)
    #
    avg_arr = np.array(avg_list, dtype=np.float64)
    err_arr = np.array(err_list, dtype=np.float64)
    #
    # Assuming true average be zero.
    actual_err_arr = avg_arr
    jk_actual_err_arr = np.sqrt(q.jackknife(actual_err_arr**2))
    jk_err_arr = np.sqrt(q.jackknife(err_arr**2))
    #
    q.json_results_append("Test results")
    q.json_results_append(f"data_list_size={data_list_size}")
    q.json_results_append(f"block_size={block_size}")
    q.json_results_append(f"compute_avg_err={compute_avg_err.__qualname__}")
    #
    q.displayln_info(
        0,
        "\n",
        q.show_val_err(q.avg_err(avg_arr)),
        "\n",
        q.show_val_err(q.jk_avg_err(jk_err_arr)),
        "\n",
        q.show_val_err(q.jk_avg_err(jk_actual_err_arr)),
        "\nDifference: ",
        q.show_val_err(q.jk_avg_err(jk_err_arr - jk_actual_err_arr)),
        f"({data_list_size} {block_size} {compute_avg_err.__qualname__})",
    )
    #
    q.json_results_append(
        f"q.avg_err(avg_arr)",
        np.array(q.avg_err(avg_arr), dtype=np.float64),
    )
    q.json_results_append(
        f"q.jk_avg_err(jk_err_arr)",
        np.array(q.jk_avg_err(jk_err_arr), dtype=np.float64),
    )
    q.json_results_append(
        f"q.jk_avg_err(jk_actual_err_arr)",
        np.array(q.jk_avg_err(jk_actual_err_arr), dtype=np.float64),
    )
    q.json_results_append(
        f"q.jk_avg_err(jk_err_arr - jk_actual_err_arr)",
        np.array(q.jk_avg_err(jk_err_arr - jk_actual_err_arr), dtype=np.float64),
    )


q.begin_with_mpi()

for data_list_size in [1, 2, 4, 7, 15, 32, ]:
    for block_size in [1, 2, 4, 7, ]:
        for compute_avg_err in [
            compute_avg_err_direct,
            compute_avg_err_jk,
            compute_avg_err_sjk_hash,
            compute_avg_err_sjk,
            compute_avg_err_rjk,
        ]:
            if data_list_size < 2 * block_size:
                continue
            jk_test(data_list_size, block_size, compute_avg_err)

q.check_log_json(__file__, check_eps=1e-10)
q.timer_display()
q.end_with_mpi()
q.displayln_info(f"CHECK: finished successfully.")
