#!/usr/bin/env python3

import sys
import qlat as q
import qlat_gpt as qg

qg.begin_with_gpt()

q.qremove_all_info("results")
q.qmkdir_info("results")

total_site = q.Coordinate([4, 4, 4, 8])
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

qg.save_gauge_field(gf, "results/ckpoint_lat.0")

gf.save("results/ckpoint_lat.1")

gf = q.GaugeField()

gf.load("results/ckpoint_lat.0")

gf.show_info()

plaq1 = gf.plaq()

assert abs(plaq - plaq1) < 1e-12

gf = qg.load_gauge_field("results/ckpoint_lat.1")

gf.show_info()

plaq1 = gf.plaq()

assert abs(plaq - plaq1) < 1e-12

q.timer_display()

qg.end_with_gpt()

q.displayln_info(f"CHECK: finished successfully.")
