#!/usr/bin/env python3

import qlat as q
import os

q.begin()

q.qremove_all_info("results")
q.qmkdir_info("results")

ld = q.LatData()

ld.set_dim_sizes([ 2, 4, 8 ])

ld.set_zero()

q.displayln_info("ld:")
q.displayln_info(ld.show())

q.displayln_info(len(ld[0,]))
q.displayln_info(len(ld[1,]))

ld.save("results/test.lat")

if q.get_id_node() == 0:
    q.displayln(os.listdir("results"))

q.timer_display()

q.end()
