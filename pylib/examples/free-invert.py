#!/usr/bin/env python3

import qlat as q

q.begin()

q.qremove_all_info("results")
q.qmkdir_info("results")
rs = q.RngState("seed")
geo = q.Geometry([4, 4, 4, 8], 1)
q.displayln_info("geo.show() =", geo.show())

gf = q.GaugeField(geo)
gf.set_unit()
gf.show_info()

fa = q.FermionAction(mass = 0.05, ls = 16, m5 = 1.0)

qinv_free = q.InverterDwfFreeField(mass = fa.mass(), m5 = fa.m5(), timer = q.Timer("py:InverterDwfFreeField"))

qinv_dwf = q.InverterDomainWall(gf = gf, fa = fa, timer = q.Timer("py:InverterDomainWall"))

src_p = q.mk_point_src(geo, [0, 0, 0, 0])

src_r = q.Prop(geo)
src_r.set_rand(rs.split("src_r"))

for src in [src_p, src_r]:
    sol = qinv_free * src

    sol1 = qinv_dwf * src

    sol_diff = sol1.copy()
    sol_diff -= sol

    q.displayln_info(sol.qnorm(), sol1.qnorm(), sol_diff.qnorm())

q.timer_display()

q.end()
