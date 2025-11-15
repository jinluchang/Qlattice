#!/usr/bin/env python3

import sys
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

rs = q.RngState("seed")

total_site = q.Coordinate([ 4, 4, 4, 8, ])
geo = q.Geometry(total_site)
gf = q.GaugeField(geo)
gf.set_rand(rs.split("gf-init"), 0.5, 10)
q.json_results_append(f"gf-init plaq", gf.plaq(), 1e-8)
q.json_results_append(f"gf-init", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

q.gf_block_stout_smear(gf, q.Coordinate(), 0.1)
q.json_results_append(f"gf_block_stout_smear after plaq", gf.plaq(), 1e-8)
q.json_results_append(f"gf_block_stout_smear after", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)
q.gf_block_stout_smear(gf, q.Coordinate([ 2, 2, 2, 2, ]), 0.1)
q.json_results_append(f"gf_block_stout_smear after block plaq", gf.plaq(), 1e-8)
q.json_results_append(f"gf_block_stout_smear after block", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)
q.gf_block_stout_smear(gf, q.Coordinate([ 4, 4, 4, 4, ]), 0.1)
q.json_results_append(f"gf_block_stout_smear after block 4 plaq", gf.plaq(), 1e-8)
q.json_results_append(f"gf_block_stout_smear after block 4", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

q.timer_display()
q.check_log_json(__file__, check_eps=1e-5)
q.end_with_mpi()
q.displayln_info(f"CHECK: finished successfully.")
