#!/usr/bin/env python3

import inspect
import numpy as np
import qlat as q

def setup_g_jk_kwargs(
    *,
    n_rand_sample,
    block_size,
    block_size_dict,
    jk_type="rjk",
    all_jk_idx=None,
    is_hash_jk_idx=True,
    jk_idx_hash_size=1024,
    is_normalizing_rand_sample=False,
    is_use_old_rand_alg=False,
):
    q.default_g_jk_kwargs["jk_type"] = jk_type
    q.default_g_jk_kwargs["eps"] = 1
    q.default_g_jk_kwargs["n_rand_sample"] = n_rand_sample
    q.default_g_jk_kwargs["is_normalizing_rand_sample"] = is_normalizing_rand_sample
    q.default_g_jk_kwargs["is_apply_rand_sample_jk_idx_blocking_shift"] = True
    q.default_g_jk_kwargs["is_use_old_rand_alg"] = is_use_old_rand_alg
    q.default_g_jk_kwargs["block_size"] = block_size
    q.default_g_jk_kwargs["block_size_dict"] = block_size_dict
    q.default_g_jk_kwargs["rng_state"] = q.RngState("rejk")
    q.default_g_jk_kwargs["all_jk_idx"] = all_jk_idx
    q.default_g_jk_kwargs["get_all_jk_idx"] = None
    q.default_g_jk_kwargs["is_hash_jk_idx"] = is_hash_jk_idx
    q.default_g_jk_kwargs["jk_idx_hash_size"] = jk_idx_hash_size
    q.default_g_jk_kwargs["all_jk_idx_set"] = set()
    q.default_g_jk_kwargs["is_sync_node"] = False

def get_max_rel_diff(a, b):
    d = np.abs(a - b)
    scale = np.maximum(np.abs(b), 1e-100)
    return float(np.max(d / scale)) if d.size > 0 else 0.0

def check_signatures():
    params = inspect.signature(q.g_mk_jk).parameters
    ok = "is_sync_node" in params and params["is_sync_node"].default is False
    q.json_results_append(f"g_mk_jk has is_sync_node (default False) = {ok}")
    assert ok
    params = inspect.signature(q.rjackknife).parameters
    ok = "is_sync_node" in params and params["is_sync_node"].default is False
    q.json_results_append(f"rjackknife has is_sync_node (default False) = {ok}")
    assert ok

def check_comm():
    comm = q.get_comm()
    q.json_results_append(f"get_comm is not None = {comm is not None}")
    q.json_results_append(
        f"get_comm().rank == get_id_node() = {comm.rank == q.get_id_node()}"
    )
    q.json_results_append(
        f"get_comm().size == get_num_node() = {comm.size == q.get_num_node()}"
    )
    assert comm is not None
    assert comm.rank == q.get_id_node()
    assert comm.size == q.get_num_node()

def check_sync(data_arr, jk_idx_list, tag):
    jk_seq = q.g_mk_jk(data_arr, jk_idx_list, is_sync_node=False)
    jk_sync = q.g_mk_jk(data_arr, jk_idx_list, is_sync_node=True)
    rel = get_max_rel_diff(jk_sync, jk_seq)
    ok = bool(jk_sync.shape == jk_seq.shape and rel < 1e-9)
    q.json_results_append(f"{tag}: sync == sequential = {ok}")
    q.json_results_append(f"{tag}: sync shape = {tuple(jk_sync.shape)}")
    assert ok
    avg, err = q.g_jk_avg_err(jk_sync)
    q.json_results_append(f"{tag} avg", np.array(avg, dtype=np.float64))
    q.json_results_append(f"{tag} err", np.array(err, dtype=np.float64))

q.begin_with_mpi()

check_signatures()
check_comm()

# ---- job_tag_1, block_size = 1, 64 trajectories, uneven 2-rank split ----
setup_g_jk_kwargs(
    n_rand_sample=1023,
    block_size=1,
    block_size_dict={"job_tag_1": 1},
)
data_arr_1 = q.RngState("sync-node-seed-1").g_rand_arr((64, 3))
jk_idx_list_1 = [("job_tag_1", traj) for traj in range(64)]

check_sync(data_arr_1, jk_idx_list_1, "bs1")

# ---- is_sync_node via the global default kwargs ----
q.default_g_jk_kwargs["is_sync_node"] = True
jk_default_sync = q.g_mk_jk(data_arr_1, jk_idx_list_1)
q.default_g_jk_kwargs["is_sync_node"] = False
jk_default_seq = q.g_mk_jk(data_arr_1, jk_idx_list_1)
ok = bool(
    jk_default_seq.shape == jk_default_sync.shape
    and get_max_rel_diff(jk_default_seq, jk_default_sync) < 1e-9
)
q.json_results_append(f"global default is_sync_node: sync == sequential = {ok}")
assert ok

# ---- job_tag_2, block_size = 4, 100 trajectories ----
setup_g_jk_kwargs(
    n_rand_sample=1023,
    block_size=4,
    block_size_dict={"job_tag_2": 4},
)
data_arr_2 = q.RngState("sync-node-seed-2").g_rand_arr((100, 3))
jk_idx_list_2 = [("job_tag_2", traj) for traj in range(100)]

check_sync(data_arr_2, jk_idx_list_2, "bs4")

# ---- mixed job tags in one call ----
setup_g_jk_kwargs(
    n_rand_sample=1023,
    block_size=1,
    block_size_dict={"job_tag_1": 1, "job_tag_2": 4},
)
data_arr_3 = np.concatenate([data_arr_1, data_arr_2], axis=0)
jk_idx_list_3 = jk_idx_list_1 + jk_idx_list_2

check_sync(data_arr_3, jk_idx_list_3, "mixed")

# ---- is_normalizing_rand_sample = True ----
setup_g_jk_kwargs(
    n_rand_sample=1023,
    block_size=1,
    block_size_dict={"job_tag_1": 1},
    is_normalizing_rand_sample=True,
)
check_sync(data_arr_1, jk_idx_list_1, "normalizing")

# ---- n_rand_sample = 0 ----
setup_g_jk_kwargs(
    n_rand_sample=0,
    block_size=1,
    block_size_dict={"job_tag_1": 1},
)
check_sync(data_arr_1, jk_idx_list_1, "n_rand_sample=0")

# ---- jk_type = "super" with is_sync_node ----
all_jk_idx = ["avg"] + [("job_tag_1", b) for b in range(16)]
setup_g_jk_kwargs(
    n_rand_sample=1023,
    block_size=1,
    block_size_dict={"job_tag_1": 1},
    jk_type="super",
    all_jk_idx=all_jk_idx,
)
check_sync(data_arr_1, jk_idx_list_1, "super")

# ---- jk_type = "super" with hash based samples ----
setup_g_jk_kwargs(
    n_rand_sample=1023,
    block_size=1,
    block_size_dict={"job_tag_1": 1},
    jk_type="super",
    is_hash_jk_idx=True,
    jk_idx_hash_size=32,
)
check_sync(data_arr_1, jk_idx_list_1, "super-hash")

# ---- is_sync_node without a qlat communicator must raise ----
setup_g_jk_kwargs(
    n_rand_sample=16,
    block_size=1,
    block_size_dict={"job_tag_1": 1},
)
saved_comm = q.get_comm()
q.set_comm(None)
comm_raised = False
try:
    q.g_mk_jk(data_arr_1, jk_idx_list_1, is_sync_node=True)
except Exception:
    comm_raised = True
q.set_comm(saved_comm)
q.json_results_append(f"is_sync_node without get_comm raises = {comm_raised}")
assert comm_raised
assert q.get_comm() is saved_comm

q.timer_display()
if q.is_test():
    q.check_log_json(__file__, check_eps=1e-10)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
