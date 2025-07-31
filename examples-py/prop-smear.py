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
        # q.Coordinate([ 8, 8, 8, 8, ]),
        # q.Coordinate([ 12, 12, 12, 12, ]),
        # q.Coordinate([ 16, 16, 16, 16, ]),
        # q.Coordinate([ 20, 20, 20, 20, ]),
        ]

def get_f_list_sig(f_list, rs, n):
    sig = np.zeros(n, dtype=np.complex128)
    for idx, sp in enumerate(f_list):
        sig += q.get_data_sig_arr(sp[:], rs.split(f"sig {idx}"), n)
    sig = q.glb_sum(sig)
    return sig

@q.timer(is_timer_fork=True)
def benchmark_prop_spatial_smear(total_site, seed):
    fname = q.get_fname()
    q.json_results_append(f"{fname}: {total_site} {seed}")
    rs = q.RngState(f"seed {fname} {seed}")
    geo = q.Geometry(total_site)
    #
    coef = 0.9375
    step = 20
    #
    gf = q.GaugeField(geo)
    gf.set_rand(rs.split(f"gf"), 0.5, 2)
    gf_sig = q.get_data_sig_arr(gf, rs, 3)
    q.json_results_append(f"gf sig", gf_sig, 1e-10)
    #
    mom = q.CoordinateD([ 0.0, 0.1, -0.2, 0.2, ])
    #
    prop = q.Prop(geo)
    prop.set_rand(rs.split(f"prop"))
    prop_sig = q.get_data_sig_arr(prop, rs, 3)
    q.json_results_append(f"prop sig", prop_sig, 1e-10)
    ss_prop = q.prop_spatial_smear(prop, gf, coef, step, mom)
    ss_prop_sig = q.get_data_sig_arr(ss_prop, rs, 3)
    q.json_results_append(f"ss_prop sig", ss_prop_sig, 1e-10)
    ss_prop = q.prop_smear(prop, gf, coef, step, mom)
    ss_prop_sig = q.get_data_sig_arr(ss_prop, rs, 3)
    q.json_results_append(f"ss_prop sig", ss_prop_sig, 1e-10)
    ss_prop = q.prop_smear(prop, gf, coef, step, mom, mode_smear=0)
    ss_prop_sig = q.get_data_sig_arr(ss_prop, rs, 3)
    q.json_results_append(f"ss_prop sig", ss_prop_sig, 1e-10)

for total_site in total_site_list:
    for seed in range(2):
        benchmark_prop_spatial_smear(total_site, seed)

q.timer_display()
q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info(f"CHECK: finished successfully.")
