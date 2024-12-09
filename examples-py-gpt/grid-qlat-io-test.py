#!/usr/bin/env python3

import qlat_grid as q

q.begin_with_grid()

rs = q.RngState()

total_site = q.Coordinate([ 8, 8, 8, 16, ])

geo = q.Geometry(total_site)

gf = q.GaugeField(geo)

gf.set_unit()

prop = q.Prop(geo)

prop.set_zero()

prop1 = prop.copy()

prop.set_rand(rs)

q.qremove_all_info(f"results")

q.save_grid_prop_float(prop, f"results/prop.grid.field")

q.load_grid_prop_float(prop1, f"results/prop.grid.field")

prop1 -= prop

q.displayln_info(f"diff ratio {q.qnorm(prop1) / q.qnorm(prop)}")

assert q.qnorm(prop1) / q.qnorm(prop) < 1e-15

q.save_grid_prop_double(prop, f"results/prop-d.grid.field")

q.load_grid_prop_double(prop1, f"results/prop-d.grid.field")

prop1 -= prop

q.displayln_info(f"CHECK: diff ratio with double {q.qnorm(prop1) / q.qnorm(prop)}")

assert q.qnorm(prop1) == 0

q.timer_display()

q.end_with_grid()

q.displayln_info(f"CHECK: finished successfully.")
