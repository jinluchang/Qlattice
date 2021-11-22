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
    [2, 2, 2, 4]]

q.begin(sys.argv, size_node_list)

q.qremove_all_info("results")
q.qmkdir_info("results")

total_site = [ 4, 4, 4, 8, ]
geo = q.Geometry(total_site, 1)
q.displayln_info("geo.show() =", geo.show())
rs = q.RngState("seed")

gf = q.GaugeField(geo)

gf.set_rand(rs.split("gf-init"), 0.5, 10)

gf.show_info()

gf_f = gf.copy()

flow_time = 1.0
t = 0
for i in range(4):
    q.gf_wilson_flow(gf_f, flow_time, 50, existing_flow_time = t, c1 = -0.331)
    gf_f.show_info()
    t += flow_time
    topo = q.gf_topology(gf_f)
    q.displayln_info(f"t={t} topo={topo}")

q.timer_display()

q.end()
