#!/usr/bin/env python3

import sys
import qlat as q

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

total_site = q.Coordinate([ 4, 4, 4, 8, ])
geo = q.Geometry(total_site)
q.displayln_info("CHECK: geo.show() =", geo.show())
rs = q.RngState("seed")

gf = q.GaugeField(geo)

q.displayln_info("CHECK: ", gf.geo.show())

q.set_unit(gf)

q.gf_show_info(gf)

gf.set_rand(rs.split("gf-init"), 0.3, 1)

gf.show_info()

plaq = gf.plaq()

gf.save("results/ckpoint_lat.0")

gf = q.GaugeField()

gf.load("results/ckpoint_lat.0")

gf.show_info()

plaq1 = gf.plaq()

q.displayln_info(f"CHECK: plaq: {plaq:.12E} {plaq1:.12E}")

assert abs(plaq - plaq1) < 1e-12

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
