#!/usr/bin/env python3

import qlat as q
import numpy as np

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 1, 2],
        [1, 1, 1, 4],
        [1, 1, 1, 8],
        [2, 2, 2, 2],
        [2, 2, 2, 4],
        ]

q.begin_with_mpi(size_node_list)

q.qremove_all_info("results")
q.qmkdir_info("results")

num_node = q.get_num_node()

q.displayln_info(f"CHECK: num_node={num_node}")

a = 12 + q.get_id_node()
gs_a = q.glb_sum(a)

q.displayln_info(f"CHECK: {a} {gs_a}")

b = 12.4 + q.get_id_node()
gs_b = q.glb_sum(b)

q.displayln_info(f"CHECK: {b} {gs_b}")

c = np.arange(3.0) + 1.1 + q.get_id_node()
gs_c = q.glb_sum(c)

q.displayln_info(f"CHECK: {c.tolist()} {gs_c.tolist()}")

d = np.arange(3.0) + 1.1 + q.get_id_node() * 1.0j
gs_d = q.glb_sum(d)

q.displayln_info(f"CHECK: {d.tolist()} {gs_d.tolist()}")

e = (np.arange(6.0) + 1.1 + q.get_id_node()).reshape(3,2)
gs_e = q.glb_sum(e)

q.displayln_info(f"CHECK: {e.tolist()} {gs_e.tolist()}")

f = (np.arange(6.0) + 1.1 + q.get_id_node() * 1.0j).reshape(3,2)
gs_f = q.glb_sum(f)

q.displayln_info(f"CHECK: {f.tolist()} {gs_f.tolist()}")

q.timer_display()
q.end_with_mpi()
q.displayln_info(f"CHECK: finished successfully.")
