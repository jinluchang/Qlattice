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

total_site_list = [
        q.Coordinate([ 4, 4, 4, 4, ]),
        q.Coordinate([ 6, 6, 6, 6, ]),
        q.Coordinate([ 8, 8, 8, 8, ]),
        ]

multiplicity_list = [
        1, 2, 3, 4, 5, 7, 9, 11, 13, 16, 20, 40,
        ]

for total_site in total_site_list:
    for seed in range(2):
        q.json_results_append(f"smear {total_site} {seed}")
        rs = q.RngState(f"{total_site} {seed}")
        geo = q.Geometry(total_site)
        prop = q.Prop(geo)
        prop.set_rand(rs)
        gf = q.GaugeField(geo)
        gf.set_rand(rs)
        gf_ape = q.gf_spatial_ape_smear(gf, 0.5, 3)
        gf1 = q.mk_left_expanded_gauge_field(gf_ape)
        coef = 0.9375
        step = 3
        prop0 = q.prop_smear(prop, gf1, coef, step, mode_smear=0)
        prop1 = q.prop_smear(prop, gf1, coef, step)
        q.json_results_append(f"smear {total_site} {seed} prop", q.get_data_sig_arr(prop, rs, 2), 1e-12)
        q.json_results_append(f"smear {total_site} {seed} prop0", q.get_data_sig_arr(prop0, rs, 2), 1e-12)
        q.json_results_append(f"smear {total_site} {seed} prop1", q.get_data_sig_arr(prop1, rs, 2), 1e-12)
        assert q.qnorm(prop0[:] - prop1[:]) < 1e-15

for is_forward in [ True, False, ]:
    for is_only_spatial in [ False, True, ]:
        fft0 = q.mk_fft(is_forward=is_forward, is_only_spatial=is_only_spatial, mode_fft=0)
        fft1 = q.mk_fft(is_forward=is_forward, is_only_spatial=is_only_spatial, mode_fft=1)
        for total_site in total_site_list:
            for multiplicity in multiplicity_list:
                for seed in range(1):
                    q.json_results_append(f"fft {total_site} {multiplicity} {seed}")
                    rs = q.RngState(f"{total_site} {multiplicity} {seed}")
                    geo = q.Geometry(total_site)
                    f = q.FieldComplexD(geo, multiplicity)
                    f.set_rand(q.RngState())
                    f0 = fft0 * f
                    f1 = fft1 * f
                    q.json_results_append(f"fft {total_site} {multiplicity} {seed} f", q.get_data_sig(f, rs), 1e-12)
                    q.json_results_append(f"fft {total_site} {multiplicity} {seed} f0", q.get_data_sig(f0, rs), 1e-12)
                    q.json_results_append(f"fft {total_site} {multiplicity} {seed} f1", q.get_data_sig(f1, rs), 1e-12)
                    assert q.qnorm(f1[:] - f0[:]) < 1e-15

q.timer_display()
q.check_log_json(__file__)
q.end_with_mpi()
q.displayln_info(f"CHECK: finished successfully.")
