#!/usr/bin/env python3

import sys
import qlat as q

size_node_list = [
    [1, 1, 1, 1],
    [1, 1, 1, 2],
    [1, 1, 1, 4],
    [1, 1, 1, 8],
    [2, 2, 2, 2],
    [2, 2, 2, 4]]

q.begin(sys.argv, size_node_list)

q.qremove_all_info("results")
q.qmkdir_info("results")

total_site = [4, 4, 4, 8]
geo = q.Geometry(total_site, 1)
q.displayln_info("geo.show() =", geo.show())
rs = q.RngState("seed")

gf = q.GaugeField(geo)

q.displayln_info(gf.geo().show_all())

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

q.displayln_info(f"CHECK: plaq: {plaq:.14E} {plaq1:.14E}")

assert abs(plaq - plaq1) < 1e-12

q.timer_display()

q.displayln_info(f"CHECK: finished successfully.")

q.end()
