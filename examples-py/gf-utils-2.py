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
gf0 = gf.copy()

gf @= gf0
q.gf_block_stout_smear(gf, q.Coordinate(), 0.1)
q.json_results_append(f"gf_block_stout_smear after plaq", gf.plaq(), 1e-8)
q.json_results_append(f"gf_block_stout_smear after", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

gf @= gf0
q.gf_block_stout_smear(gf, q.Coordinate([ 2, 2, 2, 2, ]), 0.1)
q.json_results_append(f"gf_block_stout_smear after block plaq", gf.plaq(), 1e-8)
q.json_results_append(f"gf_block_stout_smear after block", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

gf @= gf0
q.gf_block_stout_smear(gf, q.Coordinate([ 4, 4, 4, 4, ]), 0.1)
q.json_results_append(f"gf_block_stout_smear after block 4 plaq", gf.plaq(), 1e-8)
q.json_results_append(f"gf_block_stout_smear after block 4", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

gf @= gf0
gf = q.field_shift(gf, q.Coordinate([ 1, 2, 3, 4, ]))
q.json_results_append(f"field_shift plaq", gf.plaq(), 1e-8)
q.json_results_append(f"field_shift", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

gf @= gf0
q.gf_local_stout_smear(gf, 0.1)
q.json_results_append(f"gf_local_stout_smear after block 4 plaq", gf.plaq(), 1e-8)
q.json_results_append(f"gf_local_stout_smear after block 4", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

fc = q.FieldChar(geo, 1)
fc[:] = np.arange(len(fc[:]))[:, None]
fc_init = fc[:].copy()
fc_list = q.shuffle_field_char(fc, q.Coordinate([ 2, 2, 2, 4, ]))
q.json_results_append(f"shuffle_field_char {len(fc_list)}")
fc[:] = 0
q.shuffle_field_char_back(fc, fc_list, q.Coordinate([ 2, 2, 2, 4, ]))
assert np.all(fc[:] == fc_init)

gf @= gf0
gf_list = q.shuffle_field(gf, q.Coordinate([ 2, 2, 2, 4, ]))
q.json_results_append(f"shuffle_field {len(gf_list)}")
for gf_local in gf_list:
    q.gf_local_stout_smear(gf_local, 0.1)
q.shuffle_field_back(gf, gf_list, q.Coordinate([ 2, 2, 2, 4, ]))
q.json_results_append(f"gf_local_stout_smear plaq", gf.plaq(), 1e-8)
q.json_results_append(f"gf_local_stout_smear", q.get_data_sig_arr(gf, q.RngState(), 3), 1e-8)

q.timer_display()
q.check_log_json(__file__, check_eps=1e-5)
q.end_with_mpi()
q.displayln_info(f"CHECK: finished successfully.")
