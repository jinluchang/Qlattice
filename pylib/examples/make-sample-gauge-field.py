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

total_site = [4, 4, 4, 8]

geo = q.Geometry(total_site)

q.displayln_info(geo.show_all())

gf = q.GaugeField(geo)

q.displayln_info(gf.geo().show_all())

q.set_unit(gf)

q.gf_show_info(gf)

rs = q.RngState("gf-init")

gf.set_rand(rs, 0.3, 1)

q.gf_show_info(gf)

q.timer_display()

q.end()
