#!/usr/bin/env python3

import qlat as q

q.begin_with_mpi()

total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
q.displayln_info("CHECK: geo.show() =", geo.show())
rs = q.RngState("seed")

gf = q.GaugeField(geo)
gf.set_unit()
gf.show_info()

fa = q.FermionAction(mass=0.05, ls=16, m5=1.0)

qinv_free = q.InverterDwfFreeField(mass=fa.mass(), m5=fa.m5(), qtimer=q.Timer("py:InverterDwfFreeField"))

qinv_dwf = q.InverterDomainWall(gf=gf, fa=fa, qtimer=q.Timer("py:InverterDomainWall"))

src_p = q.mk_point_src(geo, q.Coordinate([0, 0, 0, 0]))

src_r = q.Prop(geo)
src_r.set_rand(rs.split("src_r"))

for src in [src_p, src_r]:
    sol = qinv_free * src

    sol1 = qinv_dwf * src

    sol_diff = sol1.copy()
    sol_diff -= sol

    q.displayln_info(f"CHECK: {sol.qnorm():.13E}", f"{sol1.qnorm():.9E}", f"{sol_diff.qnorm():.3E}")
    assert sol_diff.qnorm() < 1e-7

    sol2 = q.free_invert(src, mass=fa.mass(), m5=fa.m5(), momtwist=q.CoordinateD([0, 0, 0, 0]))
    sol2_diff = sol2.copy()
    sol2_diff -= sol
    assert sol2_diff.qnorm() < 1e-9

q.timer_display()
q.end_with_mpi()
q.displayln_info(f"CHECK: finished successfully.")
