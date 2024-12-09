#!/usr/bin/env python3

import qlat_gpt as qg
import qlat as q
import gpt as g

qg.begin_with_gpt()

q.qremove_all_info("results")
q.qmkdir_info("results")
total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
q.displayln_info("CHECK: geo.show() =", geo.show())
rs = q.RngState("seed")

gf = q.GaugeField(geo)
gf.set_unit()
gf.show_info()

mass = 0.05

m5 = 1.0

qinv = q.InverterDwfFreeField(mass = mass, m5 = m5, qtimer = q.Timer("py:InverterDwfFreeField"))

mobius_params = {
        "mass": mass,
        "M5": m5,
        "b": 1.0,
        "c": 0.0,
        "Ls": 32,
        "boundary_phases": [1.0, 1.0, 1.0, 1.0],
        }

gpt_gf = qg.gpt_from_qlat(gf)
qm = g.qcd.fermion.mobius(gpt_gf, mobius_params)
pc = g.qcd.fermion.preconditioner
inv = g.algorithms.inverter
cg = inv.cg({"eps": 1e-11, "maxiter": 10000})
slv_5d = inv.preconditioned(pc.eo2_ne(), cg)
slv_qm = qm.propagator(slv_5d)
ginv = qg.InverterGPT(inverter = slv_qm, qtimer = q.Timer("py:InverterGPT"))

src_p = q.mk_point_src(geo, q.Coordinate([0, 0, 0, 0]))

src_r = q.Prop(geo)
src_r.set_rand(rs.split("src_r"))

for src in [src_p, src_r]:
    sol = qinv * src

    sol1 = ginv * src

    sol_diff = sol1.copy()
    sol_diff -= sol

    q.displayln_info(f"CHECK: {sol.qnorm():.10E} {sol1.qnorm():.10E}")
    assert sol_diff.qnorm() < 1e-15

q.timer_display()

qg.end_with_gpt()

q.displayln_info(f"CHECK: finished successfully.")
