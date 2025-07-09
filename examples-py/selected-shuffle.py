#!/usr/bin/env python3

import qlat as q
import numpy as np
import sys

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 1, 2],
        [1, 1, 1, 4],
        [1, 1, 1, 8],
        [2, 2, 2, 2],
        [2, 2, 2, 4],
        ]

q.begin_with_mpi(size_node_list)

total_site_list = [
        q.Coordinate([ 4, 4, 4, 4, ]),
        q.Coordinate([ 6, 6, 6, 6, ]),
        q.Coordinate([ 8, 8, 8, 8, ]),
        ]

multiplicity_list = [
        1, 2, 3,
        ]

@q.timer
def selected_shuffle_random(total_site, multiplicity, seed):
    fname = q.get_fname()
    q.json_results_append(f"{fname}: {total_site} {multiplicity} {seed}")
    rs = q.RngState(f"seed {fname} {seed}")
    n_points = total_site.volume() // 16
    q.json_results_append(f"n_points={n_points}")
    psel = q.PointsSelection()
    psel.set_rand(total_site, n_points, rs.split("psel"))
    q.json_results_append(f"hash(psel)={q.hash_sha256(psel)}")
    fsel = q.FieldSelection(psel)
    psel_l = q.PointsSelection(fsel)
    ssp = q.SelectedShufflePlan(psel_l, rs.split("ssp"))
    psel_s = q.PointsSelection(psel_l, ssp)
    q.displayln_info(f"len(psel)={len(psel)} ; psel={psel}");
    psel_str = f"len(psel_s)={len(psel_s)} ; psel_s={psel_s}"
    psel_str_list = q.get_comm().allgather(psel_str)
    for id_node, psel_str in enumerate(psel_str_list):
        q.displayln_info(f"id_node={id_node} ; {psel_str}");
    q.json_results_append(f"hash(psel_s)={q.hash_sha256(psel_s)}")

for total_site in total_site_list:
    for multiplicity in multiplicity_list:
        for seed in range(1):
            selected_shuffle_random(total_site, multiplicity, seed)

q.timer_display()
q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info(f"CHECK: finished successfully.")
