#!/usr/bin/env python3

import qlat_cps as q

total_site = q.Coordinate([ 8, 8, 8, 16, ])

q.begin_with_cps(total_site)

rs = q.RngState()

geo = q.Geometry(total_site)

gf = q.GaugeField(geo)

gf.set_unit()

prop = q.Prop(geo)

prop.set_zero()

prop1 = prop.copy()

prop.set_rand(rs)

q.qremove_all_info(f"results")

q.save_cps_prop_double(prop, f"results/prop-d.cps.field")

q.load_cps_prop_double(prop1, f"results/prop-d.cps.field")

prop1 -= prop

q.displayln_info(f"CHECK: diff ratio with double {q.qnorm(prop1) / q.qnorm(prop)}")

assert q.qnorm(prop1) == 0

q.timer_display()

q.end_with_cps()

q.displayln_info(f"CHECK: finished successfully.")
