#!/usr/bin/env python3

check_eps = 1e-10

import sys
import qlat as q
import numpy as np

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
    [2, 2, 2, 2],
    [2, 2, 2, 4]]

q.begin_with_mpi(size_node_list)

q.qremove_all_info("results")
q.qmkdir_info("results")

total_site = q.Coordinate([ 4, 4, 4, 8, ])
geo = q.Geometry(total_site)
q.displayln_info("CHECK: geo.show() =", geo.show())
rs = q.RngState("seed")

f1 = q.FieldComplexD(geo, 3)
f2 = q.FieldComplexD(geo, 5)
f1.set_rand(rs.split("f1"))
f2.set_rand(rs.split("f2"))

q.json_results_append("f1 data sig", q.get_data_sig(f1, rs), check_eps)
q.json_results_append("f2 data sig", q.get_data_sig(f2, rs), check_eps)

radius = 2.0
sf1 = q.smear_field(f1, radius)
sf2 = q.smear_field(f2, radius, is_only_spatial=True)

q.json_results_append("sf1 data sig", q.get_data_sig(sf1, rs), check_eps)
q.json_results_append("sf2 data sig", q.get_data_sig(sf2, rs), check_eps)

idx1 = np.array([ 0, 1, ], dtype=np.int32)
idx2 = np.array([ 0, 2, ], dtype=np.int32)

ff = q.field_convolution(f1, f2, idx1, idx2)

q.json_results_append("ff data sig", q.get_data_sig(ff, rs), check_eps)

ff3d = q.field_convolution(f1, f2, idx1, idx2, is_only_spatial=True)

q.json_results_append("ff3d data sig", q.get_data_sig(ff3d, rs), check_eps)

def check_convolution_4d(f1, f2, idx1, idx2, ff, xg_rel):
    fsel = q.FieldSelection(geo, 0)
    psel = fsel.to_psel()
    s = 0
    for xg in psel:
        xg1 = xg
        xg2 = (xg + xg_rel) % total_site
        v1 = f1.get_elems_xg(xg1)[0, idx1]
        v2 = f2.get_elems_xg(xg2)[0, idx2]
        s += v1 * v2
    s_ref = ff.get_elems_xg(xg_rel)[0]
    diff = s - s_ref
    assert q.qnorm(diff) <= 1e-12 * q.qnorm(s)
    assert q.qnorm(diff) <= 1e-12 * q.qnorm(s_ref)

def check_convolution_3d(f1, f2, idx1, idx2, ff, xg_rel):
    t = xg_rel[3]
    xg_rel_3d = xg_rel.copy()
    xg_rel_3d[3] = 0
    fsel = q.FieldSelection(geo, 0)
    psel = fsel.to_psel()
    s = 0
    for xg in psel:
        if xg[3] != t:
            continue
        xg1 = xg
        xg2 = (xg + xg_rel_3d) % total_site
        v1 = f1.get_elems_xg(xg1)[0, idx1]
        v2 = f2.get_elems_xg(xg2)[0, idx2]
        s += v1 * v2
    s_ref = ff.get_elems_xg(xg_rel)[0]
    diff = s - s_ref
    assert q.qnorm(diff) <= 1e-12 * q.qnorm(s)
    assert q.qnorm(diff) <= 1e-12 * q.qnorm(s_ref)

rsi = rs.split("c-rand-gen")
for i in range(4):
    xg_rel = rsi.c_rand_gen(total_site)
    check_convolution_4d(f1, f2, idx1, idx2, ff, xg_rel)
    check_convolution_3d(f1, f2, idx1, idx2, ff3d, xg_rel)

q.check_log_json(__file__)

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
