#!/usr/bin/env python3

import numpy as np
import qlat as q

def setup_g_jk_kwargs(
    *,
    n_rand_sample,
    block_size=1,
    block_size_dict=None,
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
    q.default_g_jk_kwargs["block_size_dict"] = (
        dict() if block_size_dict is None else block_size_dict
    )
    q.default_g_jk_kwargs["rng_state"] = q.RngState("rejk")
    q.default_g_jk_kwargs["all_jk_idx"] = all_jk_idx
    q.default_g_jk_kwargs["get_all_jk_idx"] = None
    q.default_g_jk_kwargs["is_hash_jk_idx"] = is_hash_jk_idx
    q.default_g_jk_kwargs["jk_idx_hash_size"] = jk_idx_hash_size
    q.default_g_jk_kwargs["all_jk_idx_set"] = set()
    q.default_g_jk_kwargs["is_sync_node"] = False

def get_local_part(x, id_node, num_node):
    return x[slice(*q.get_distributed_range(len(x), id_node, num_node))]

def get_max_rel_diff(a, b):
    d = np.abs(a - b)
    scale = np.maximum(np.abs(b), 1e-100)
    return float(np.max(d / scale)) if d.size > 0 else 0.0

def check_distributed(data_arr, jk_idx_list, tag):
    """
    Check that the distributed result of ``g_mk_jk_distributed`` agrees with
    the sequential ``g_mk_jk`` and that concatenating the local parts in the
    order of the nodes reproduces the sequential result.
    """
    comm = q.get_comm()
    rank, num_node = comm.rank, comm.size
    jk_seq = q.g_mk_jk(data_arr, jk_idx_list)
    jk_local = q.g_mk_jk_distributed(
        get_local_part(data_arr, rank, num_node),
        get_local_part(jk_idx_list, rank, num_node),
    )
    jk_local_list = comm.allgather(jk_local)
    jk_glb = np.concatenate(jk_local_list)
    total_size = q.g_jk_size()
    ok_split = all(
        len(jk_local_list[r])
        == q.get_distributed_range(total_size, r, num_node)[1]
        - q.get_distributed_range(total_size, r, num_node)[0]
        for r in range(num_node)
    )
    ok_shape = jk_glb.shape == jk_seq.shape
    rel = get_max_rel_diff(jk_glb, jk_seq)
    ok = bool(ok_split and ok_shape and rel < 1e-9)
    q.json_results_append(f"{tag}: nodes own all the samples = {ok_split}")
    q.json_results_append(f"{tag}: distributed shape = {tuple(jk_glb.shape)}")
    q.json_results_append(f"{tag}: distributed matches g_mk_jk = {ok}")
    assert ok_split
    assert ok_shape
    assert rel < 1e-9, (tag, rel)

q.begin_with_mpi()

comm = q.get_comm()
q.json_results_append(f"get_comm().size = {comm.size}")
q.json_results_append(
    f"get_comm().size == get_num_node() = {comm.size == q.get_num_node()}"
)

# ---- rjk, block_size = 1, 64 trajectories (equal input split) ----
setup_g_jk_kwargs(
    n_rand_sample=1023,
    block_size=1,
    block_size_dict={"job_tag_1": 1},
)
data_arr_1 = q.RngState("distributed-seed-1").g_rand_arr((64, 3))
jk_idx_list_1 = [("job_tag_1", traj) for traj in range(64)]

check_distributed(data_arr_1, jk_idx_list_1, "bs1")

avg, err = q.g_jk_avg_err(
    np.concatenate(
        comm.allgather(
            q.g_mk_jk_distributed(
                get_local_part(data_arr_1, comm.rank, comm.size),
                get_local_part(jk_idx_list_1, comm.rank, comm.size),
            )
        )
    )
)
q.json_results_append("bs1 avg", np.array(avg, dtype=np.float64))
q.json_results_append("bs1 err", np.array(err, dtype=np.float64))

# ---- the distributed result also matches the is_sync_node result ----
jk_sync = q.g_mk_jk(data_arr_1, jk_idx_list_1, is_sync_node=True)
jk_dist_local = q.g_mk_jk_distributed(
    get_local_part(data_arr_1, comm.rank, comm.size),
    get_local_part(jk_idx_list_1, comm.rank, comm.size),
)
jk_dist = np.concatenate(comm.allgather(jk_dist_local))
rel = get_max_rel_diff(jk_dist, jk_sync)
q.json_results_append(f"bs1: distributed matches is_sync_node = {rel < 1e-9}")
assert rel < 1e-9, rel

# ---- rjk, block_size = 4, uneven input split (65 trajectories) ----
setup_g_jk_kwargs(
    n_rand_sample=1023,
    block_size=4,
    block_size_dict={"job_tag_2": 4},
)
data_arr_2 = q.RngState("distributed-seed-2").g_rand_arr((65, 3))
jk_idx_list_2 = [("job_tag_2", traj) for traj in range(65)]

check_distributed(data_arr_2, jk_idx_list_2, "bs4-uneven")

# ---- rjk, 100 trajectories, block_size = 4 ----
data_arr_3 = q.RngState("distributed-seed-3").g_rand_arr((100, 3))
jk_idx_list_3 = [("job_tag_2", traj) for traj in range(100)]

check_distributed(data_arr_3, jk_idx_list_3, "bs4")

# ---- rjk, mixed job tags in one call ----
setup_g_jk_kwargs(
    n_rand_sample=1023,
    block_size=1,
    block_size_dict={"job_tag_1": 1, "job_tag_2": 4},
)
data_arr_4 = np.concatenate([data_arr_1, data_arr_3], axis=0)
jk_idx_list_4 = jk_idx_list_1 + jk_idx_list_3

check_distributed(data_arr_4, jk_idx_list_4, "mixed")

# ---- rjk, is_normalizing_rand_sample = True ----
setup_g_jk_kwargs(
    n_rand_sample=1023,
    block_size=1,
    block_size_dict={"job_tag_1": 1},
    is_normalizing_rand_sample=True,
)
check_distributed(data_arr_1, jk_idx_list_1, "normalizing")

# ---- rjk, is_use_old_rand_alg = "v1" ----
setup_g_jk_kwargs(
    n_rand_sample=1023,
    block_size=1,
    block_size_dict={"job_tag_1": 1},
    is_use_old_rand_alg="v1",
)
check_distributed(data_arr_1, jk_idx_list_1, "old-rand-alg")

# ---- rjk, n_rand_sample = 0 (only the last node owns the average) ----
setup_g_jk_kwargs(
    n_rand_sample=0,
    block_size=1,
    block_size_dict={"job_tag_1": 1},
)
check_distributed(data_arr_1, jk_idx_list_1, "n_rand_sample=0")

# ---- rjk, n_rand_sample = 1 (one sample per node) ----
setup_g_jk_kwargs(
    n_rand_sample=1,
    block_size=1,
    block_size_dict={"job_tag_1": 1},
)
check_distributed(data_arr_1, jk_idx_list_1, "n_rand_sample=1")

# ---- rjk, a single trajectory (some nodes have no data) ----
setup_g_jk_kwargs(
    n_rand_sample=1023,
    block_size=1,
    block_size_dict={"job_tag_1": 1},
)
data_arr_5 = q.RngState("distributed-seed-5").g_rand_arr((1, 2))
jk_idx_list_5 = [("job_tag_1", 0)]

check_distributed(data_arr_5, jk_idx_list_5, "single-traj")

# ---- super, all_jk_idx provided ----
all_jk_idx = ["avg"] + [("job_tag_1", b) for b in range(16)]
setup_g_jk_kwargs(
    n_rand_sample=1023,
    block_size=1,
    block_size_dict={"job_tag_1": 1},
    jk_type="super",
    all_jk_idx=all_jk_idx,
)
check_distributed(data_arr_1, jk_idx_list_1, "super")

# ---- super, hash based samples ----
setup_g_jk_kwargs(
    n_rand_sample=1023,
    block_size=1,
    block_size_dict={"job_tag_1": 1},
    jk_type="super",
    is_hash_jk_idx=True,
    jk_idx_hash_size=32,
)
check_distributed(data_arr_1, jk_idx_list_1, "super-hash")

# ---- the jk_type specific distributed functions ----
setup_g_jk_kwargs(
    n_rand_sample=1023,
    block_size=1,
    block_size_dict={"job_tag_1": 1},
)
data_local_1 = get_local_part(data_arr_1, comm.rank, comm.size)
jk_idx_local_1 = get_local_part(jk_idx_list_1, comm.rank, comm.size)
jk_r_dist = np.concatenate(
    comm.allgather(
        q.rjackknife_distributed(
            data_local_1,
            jk_idx_local_1,
            n_rand_sample=1023,
            rng_state=q.RngState("rejk"),
            jk_blocking_func=q.jk_blocking_func_default,
            is_normalizing_rand_sample=False,
            is_apply_rand_sample_jk_idx_blocking_shift=True,
            is_use_old_rand_alg=False,
            eps=1,
        )
    )
)
jk_r_seq = q.rjackknife(
    data_arr_1,
    jk_idx_list_1,
    n_rand_sample=1023,
    rng_state=q.RngState("rejk"),
    jk_blocking_func=q.jk_blocking_func_default,
    is_normalizing_rand_sample=False,
    is_apply_rand_sample_jk_idx_blocking_shift=True,
    is_use_old_rand_alg=False,
    eps=1,
)
rel = get_max_rel_diff(jk_r_dist, jk_r_seq)
q.json_results_append(f"rjackknife_distributed matches rjackknife = {rel < 1e-9}")
assert rel < 1e-9, rel

all_jk_idx = ["avg"] + [("job_tag_1", b) for b in range(16)]
setup_g_jk_kwargs(
    n_rand_sample=1023,
    block_size=1,
    block_size_dict={"job_tag_1": 1},
    jk_type="super",
    all_jk_idx=all_jk_idx,
)
jk_s_dist = np.concatenate(
    comm.allgather(
        q.sjackknife_distributed(
            data_local_1,
            jk_idx_local_1,
            all_jk_idx=all_jk_idx,
            rng_state=q.RngState("rejk"),
            jk_blocking_func=q.jk_blocking_func_default,
            eps=1,
        )
    )
)
jk_s_seq = q.sjackknife(
    data_arr_1,
    jk_idx_list_1,
    all_jk_idx=all_jk_idx,
    rng_state=q.RngState("rejk"),
    jk_blocking_func=q.jk_blocking_func_default,
    eps=1,
)
rel = get_max_rel_diff(jk_s_dist, jk_s_seq)
q.json_results_append(f"sjackknife_distributed matches sjackknife = {rel < 1e-9}")
assert rel < 1e-9, rel

# ---- the gather / reduce-scatter helpers accept non-contiguous views ----
# The send buffers of the collectives must be contiguous; the 1-D strided
# views below used to fail inside mpi4py with "ndarray is not contiguous",
# because reshape(-1) is not a copy for an array which is already 1-D.
total_g = 5 * comm.size
rows_g = np.zeros((5, 4))
rows_g[:, 1] = np.arange(5) + 10.0 * comm.rank
jk_col = q.get_gathered_jk_arr(rows_g[:, 1], comm, comm.size)
jk_col_cpy = q.get_gathered_jk_arr(np.ascontiguousarray(rows_g[:, 1]), comm, comm.size)
ok = bool(jk_col.shape == (total_g,) and np.array_equal(jk_col, jk_col_cpy))
q.json_results_append(f"gather of a strided 1-D view = {ok}")
assert ok

buf_g = np.zeros(10)
jk_stride = q.get_gathered_jk_arr(buf_g[::2], comm, comm.size)
jk_stride_cpy = q.get_gathered_jk_arr(np.ascontiguousarray(buf_g[::2]), comm, comm.size)
ok = bool(jk_stride.shape == (total_g,) and np.array_equal(jk_stride, jk_stride_cpy))
q.json_results_append(f"gather of a stride 2 1-D view = {ok}")
assert ok

total_r = 5 * comm.size
i_start_r, i_end_r = q.get_distributed_range(total_r, comm.rank, comm.size)
rows_r = np.zeros((total_r, 3))
rows_r[:, 0] = np.arange(total_r) + 100.0 * comm.rank
rs_a = q.get_reduce_scattered_jk_arr(rows_r[:, 0], comm, comm.rank, comm.size, 0.0)
rs_b = q.get_reduce_scattered_jk_arr(
    np.ascontiguousarray(rows_r[:, 0]), comm, comm.rank, comm.size, 0.0
)
rs_expected = np.sum(
    [np.arange(total_r) + 100.0 * r for r in range(comm.size)], axis=0
)[i_start_r:i_end_r]
ok = bool(
    rs_a.shape == (i_end_r - i_start_r,)
    and np.array_equal(rs_a, rs_b)
    and np.array_equal(rs_a, rs_expected)
)
q.json_results_append(f"reduce scatter of a strided 1-D view = {ok}")
assert ok

# ---- g_mk_jk_distributed without a qlat communicator must raise ----
setup_g_jk_kwargs(
    n_rand_sample=16,
    block_size=1,
    block_size_dict={"job_tag_1": 1},
)
saved_comm = q.get_comm()
q.set_comm(None)
comm_raised = False
try:
    q.g_mk_jk_distributed(data_arr_1, jk_idx_list_1)
except Exception:
    comm_raised = True
q.set_comm(saved_comm)
q.json_results_append(f"distributed without get_comm raises = {comm_raised}")
assert comm_raised
assert q.get_comm() is saved_comm

q.timer_display()
if q.is_test():
    q.check_log_json(__file__, check_eps=1e-10)
q.end_with_mpi()
q.displayln_info("CHECK: finished successfully.")
