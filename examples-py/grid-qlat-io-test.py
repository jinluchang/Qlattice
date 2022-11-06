#!/usr/bin/env python3

import qlat_grid as qgi
import qlat as q

qgi.begin_with_grid()

rs = q.RngState()

total_site = [ 8, 8, 8, 16, ]

geo = q.Geometry(total_site, 1)

gf = q.GaugeField(geo)

gf.set_unit()

prop = q.Prop(geo)

prop.set_zero()

prop1 = prop.copy()

prop.set_rand(rs)

q.qremove_all_info(f"results")

qgi.save_prop_float(prop, f"results/prop.grid.field")

qgi.load_prop_float(prop1, f"results/prop.grid.field")

prop1 -= prop

q.displayln_info(f"diff ratio {q.qnorm(prop1) / q.qnorm(prop)}")

assert q.qnorm(prop1) / q.qnorm(prop) < 1e-15

q.timer_display()

q.displayln_info(f"CHECK: finished successfully.")

qgi.end_with_grid()
