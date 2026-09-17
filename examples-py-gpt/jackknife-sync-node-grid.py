#!/usr/bin/env python3

import numpy as np
import qlat_grid as q

q.begin_with_grid()

# ``begin_with_grid`` must set the qlat communicator (like ``begin_with_gpt``)
# so that ``q.get_comm().rank == q.get_id_node()`` even when the Grid processor
# coordinates do not follow the ``MPI_COMM_WORLD`` rank order.
q.json_results_append(f"get_comm is not None = {q.get_comm() is not None}")
q.json_results_append(
    f"get_comm().rank == get_id_node() = {q.get_comm().rank == q.get_id_node()}"
)
q.json_results_append(
    f"get_comm().size == get_num_node() = {q.get_comm().size == q.get_num_node()}"
)
assert q.get_comm() is not None
assert q.get_comm().rank == q.get_id_node()
assert q.get_comm().size == q.get_num_node()

gathered = sorted(q.get_comm().allgather(q.get_id_node()))
q.json_results_append(
    f"get_comm().allgather == range(num_node) = {gathered == list(range(q.get_num_node()))}"
)
assert gathered == list(range(q.get_num_node()))

q.default_g_jk_kwargs["jk_type"] = "rjk"
q.default_g_jk_kwargs["eps"] = 1
q.default_g_jk_kwargs["n_rand_sample"] = 1023
q.default_g_jk_kwargs["is_normalizing_rand_sample"] = False
q.default_g_jk_kwargs["is_apply_rand_sample_jk_idx_blocking_shift"] = True
q.default_g_jk_kwargs["block_size"] = 1
q.default_g_jk_kwargs["block_size_dict"] = {"job_tag": 1}
q.default_g_jk_kwargs["rng_state"] = q.RngState("rejk")
q.default_g_jk_kwargs["is_sync_node"] = False

data_arr = q.RngState("grid-sync-node-seed").g_rand_arr((32, 3))
jk_idx_list = [("job_tag", traj) for traj in range(32)]

jk_seq = q.g_mk_jk(data_arr, jk_idx_list, is_sync_node=False)
jk_sync = q.g_mk_jk(data_arr, jk_idx_list, is_sync_node=True)

ok = bool(np.array_equal(jk_seq, jk_sync))
q.json_results_append(f"grid: sync == sequential = {ok}")
q.json_results_append(f"grid: sync shape = {tuple(jk_sync.shape)}")
assert ok

q.timer_display()
if q.is_test():
    q.check_log_json(__file__, check_eps=1e-10)
q.end_with_grid()
q.displayln_info("CHECK: finished successfully.")
