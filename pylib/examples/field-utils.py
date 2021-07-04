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

total_site = [4, 4, 4, 8]
geo = q.Geometry(total_site, 1)
q.displayln_info("geo.show() =", geo.show())
rs = q.RngState("seed")

gf = q.GaugeField(geo)

gf.set_rand(rs.split("gf-init"), 0.3, 1)

gf.show_info()

f_factor = q.mk_phase_field(gf.geo(), [1, 0, 0, 0,])

gf *= f_factor

gf.show_info()

q.displayln_info(np.array(gf.get_elems([0, 0, 0, 0,])))

gf_sum = np.array(gf.glb_sum())

q.displayln_info(gf_sum)

gf_sum_tslice = np.array(gf.glb_sum_tslice())

for t in range(total_site[3]):
    if t == 0:
        gf_sum -= gf_sum_tslice[t]
    else:
        gf_sum -= gf_sum_tslice[t]

q.displayln_info(np.linalg.norm(gf_sum))

f = gf.as_complex_field()

q.displayln_info(np.array(f.get_elems([0, 0, 0, 0,])))

gf1 = q.GaugeField()

gf1.from_complex_field(f)

gf1 -= gf

q.displayln_info("diff norm", gf1.qnorm())

q.timer_display()

q.end()
