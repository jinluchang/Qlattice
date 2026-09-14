#!/usr/bin/env python3

# Tests for the cqlat MPI bootstrap functions (qlat/cqlat/init.cpp):
#     cbegin, cend
#
# ``cbegin`` is the low-level entry point that initializes qlat (and the MPI
# geometry node) and ``cend`` is its counterpart.  They are not used by
# ``q.begin_with_mpi`` / ``q.end_with_mpi`` (which go through ``qlat.mpi``), so
# this is the only place where they are exercised: the test intentionally
# bootstraps qlat with ``q.c.cbegin(id_node, [2, 1, 1, 1])`` (the branch that
# takes an existing MPI rank/size instead of initializing MPI itself) instead of
# calling ``q.begin_with_mpi``, and shuts it down with ``q.c.cend``.

import numpy as np

import qlat as q

from mpi4py import MPI

size_node = [2, 1, 1, 1]

def get_latt_size(geo):
    return tuple(int(geo.total_site[i]) for i in range(4))

def mk_global_field_arr(rs, latt_size):
    # rank independent global pattern; the local part is picked by global
    # coordinate so that the two MPI ranks hold different values
    return rs.split("f").u_rand_arr(latt_size)

def mk_field_from_global(geo, g):
    f = q.FieldRealD(geo, 1)
    buf = np.asarray(f)
    for index in range(geo.local_volume):
        xg = geo.coordinate_g_from_l(geo.coordinate_from_index(index))
        buf[index, 0] = g[xg[0], xg[1], xg[2], xg[3]]
    return f

comm = MPI.COMM_WORLD
assert comm.size == 2, comm.size
assert q.is_test()

# --- cbegin with an already initialized MPI (id_node, size_node) ---
q.c.cbegin(comm.rank, size_node)

geo = q.Geometry(q.Coordinate([4, 4, 4, 8]))
q.json_results_append(f"cqlat-begin-end: geo.show()={geo.show()}")
q.json_results_append("cqlat-begin-end: num_node", float(q.get_num_node()), 0.0)
q.json_results_append("cqlat-begin-end: id_node", float(q.get_id_node()), 0.0)
q.json_results_append("cqlat-begin-end: local_volume", float(geo.local_volume), 0.0)
assert q.get_num_node() == 2
assert q.get_size_node().volume() == 2
assert geo.local_volume * q.get_num_node() == 4 * 4 * 4 * 8

# --- qlat works normally after the cbegin bootstrap ---
latt_size = get_latt_size(geo)
g = mk_global_field_arr(q.RngState("cqlat-begin-end"), latt_size)
f = mk_field_from_global(geo, g)
f_sum = q.glb_sum(float(np.asarray(f).sum()))
q.json_results_append("cqlat-begin-end: field sum", f_sum, 1e-12)
assert abs(f_sum - float(g.sum())) < 1e-10 * max(1.0, abs(float(g.sum())))
f_sum_sq = q.glb_sum(float((np.asarray(f) ** 2).sum()))
q.json_results_append("cqlat-begin-end: field sum of squares", f_sum_sq, 1e-10)
assert abs(f_sum_sq - float((g**2).sum())) < 1e-10 * max(1.0, float((g**2).sum()))

# --- global reductions still work ---
q.json_results_append("cqlat-begin-end: total_site", float(np.prod(latt_size)), 0.0)
assert q.glb_sum(float(q.get_num_node())) == 4.0

q.timer_display()
if q.is_test():
    q.check_log_json(__file__)

# --- cend, then finish MPI by hand (begin_with_mpi was never used) ---
q.c.cend(False)
MPI.Finalize()
q.displayln_info("CHECK: finished successfully.")
