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
rs_sig = q.RngState()

total_site = q.Coordinate([ 4, 4, 4, 8, ])
geo = q.Geometry(total_site)
gf = q.GaugeField(geo)
gf.set_rand(rs.split("gf-init"), 0.5, 10)

q.json_results_append(f"total_site={total_site}")

f_plaq = q.gf_plaq_field(gf)
plaq_min = np.min(f_plaq.glb_min()[:])
plaq = f_plaq.glb_sum()[:].sum().item() / geo.total_volume / 6
gf_sig = q.get_data_sig_arr(gf, rs_sig, 3)
plaq_sig = q.get_data_sig_arr(f_plaq, rs_sig, 3)
q.json_results_append(f"gf plaq", plaq, 1e-10)
q.json_results_append(f"gf plaq_min", plaq_min, 1e-10)
q.json_results_append(f"gf_sig", gf_sig, 1e-10)
q.json_results_append(f"plaq_sig", plaq_sig, 1e-10)

assert abs(plaq - gf.plaq()) < 1e-10

num_step = 5
step_size = 0.1
for i in range(num_step):
    q.gf_flow_topo(gf, step_size)
q.json_results_append(f"gf_flow_topo ; num_step={num_step} ; step_size={step_size:.3f}")

f_plaq = q.gf_plaq_field(gf)
plaq_min = np.min(f_plaq.glb_min()[:])
plaq = f_plaq.glb_sum()[:].sum().item() / geo.total_volume / 6
gf_sig = q.get_data_sig_arr(gf, rs_sig, 3)
plaq_sig = q.get_data_sig_arr(f_plaq, rs_sig, 3)
q.json_results_append(f"gf plaq", plaq, 1e-10)
q.json_results_append(f"gf plaq_min", plaq_min, 1e-10)
q.json_results_append(f"gf_sig", gf_sig, 1e-10)
q.json_results_append(f"plaq_sig", plaq_sig, 1e-10)

num_step = 5
step_size = 0.1
for i in range(num_step):
    q.gf_flow_topo(gf, step_size, "Freeze")
q.json_results_append(f"gf_flow_topo ; Freeze ; num_step={num_step} ; step_size={step_size:.3f}")

f_plaq = q.gf_plaq_field(gf)
plaq_min = np.min(f_plaq.glb_min()[:])
plaq = f_plaq.glb_sum()[:].sum().item() / geo.total_volume / 6
gf_sig = q.get_data_sig_arr(gf, rs_sig, 3)
plaq_sig = q.get_data_sig_arr(f_plaq, rs_sig, 3)
q.json_results_append(f"gf plaq", plaq, 1e-10)
q.json_results_append(f"gf plaq_min", plaq_min, 1e-10)
q.json_results_append(f"gf_sig", gf_sig, 1e-10)
q.json_results_append(f"plaq_sig", plaq_sig, 1e-10)

num_step = 5
step_size = 0.1
for i in range(num_step):
    q.gf_flow_topo(gf, step_size, "Shrink")
q.json_results_append(f"gf_flow_topo ; Shrink ; num_step={num_step} ; step_size={step_size:.3f}")

f_plaq = q.gf_plaq_field(gf)
plaq_min = np.min(f_plaq.glb_min()[:])
plaq = f_plaq.glb_sum()[:].sum().item() / geo.total_volume / 6
gf_sig = q.get_data_sig_arr(gf, rs_sig, 3)
plaq_sig = q.get_data_sig_arr(f_plaq, rs_sig, 3)
q.json_results_append(f"gf plaq", plaq, 1e-10)
q.json_results_append(f"gf plaq_min", plaq_min, 1e-10)
q.json_results_append(f"gf_sig", gf_sig, 1e-10)
q.json_results_append(f"plaq_sig", plaq_sig, 1e-10)

q.timer_display()
q.check_log_json(__file__, check_eps=1e-0)
q.end_with_mpi()
q.displayln_info(f"CHECK: finished successfully.")
