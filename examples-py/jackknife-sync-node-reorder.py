#!/usr/bin/env python3

import numpy as np
import qlat as q
from mpi4py import MPI

def check_sync(data_arr, jk_idx_list, tag):
    jk_seq = q.g_mk_jk(data_arr, jk_idx_list, is_sync_node=False)
    jk_sync = q.g_mk_jk(data_arr, jk_idx_list, is_sync_node=True)
    ok = bool(np.array_equal(jk_seq, jk_sync))
    q.json_results_append(f"{tag}: sync == sequential = {ok}")
    q.json_results_append(f"{tag}: sync shape = {tuple(jk_sync.shape)}")
    assert ok

# Start qlat with a deliberately reordered id_node so that
# q.get_id_node() != MPI.COMM_WORLD.rank (as happens with begin_with_gpt and
# begin_with_grid for some Grid processor layouts), and set the qlat
# communicator exactly like the fixed begin_with_grid / begin_with_gpt do.
comm = MPI.COMM_WORLD
id_node = comm.size - 1 - comm.rank
q.begin(id_node, q.Coordinate([1, 1, 1, comm.size]))
q.set_comm(comm.Split(color=0, key=q.get_id_node()))

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
q.json_results_append(f"id_node != world.rank = {q.get_id_node() != comm.rank}")
if comm.size > 1:
    assert q.get_id_node() != comm.rank

q.default_g_jk_kwargs["jk_type"] = "rjk"
q.default_g_jk_kwargs["eps"] = 1
q.default_g_jk_kwargs["n_rand_sample"] = 1023
q.default_g_jk_kwargs["is_normalizing_rand_sample"] = False
q.default_g_jk_kwargs["is_apply_rand_sample_jk_idx_blocking_shift"] = True
q.default_g_jk_kwargs["block_size"] = 1
q.default_g_jk_kwargs["block_size_dict"] = {"job_tag": 1}
q.default_g_jk_kwargs["rng_state"] = q.RngState("rejk")
q.default_g_jk_kwargs["is_sync_node"] = False

data_arr = q.RngState("reorder-seed").g_rand_arr((32, 2))
jk_idx_list = [("job_tag", traj) for traj in range(32)]

check_sync(data_arr, jk_idx_list, "reorder")

q.timer_display()
if q.is_test():
    q.check_log_json(__file__, check_eps=1e-10)
q.end()
q.displayln_info("CHECK: finished successfully.")
