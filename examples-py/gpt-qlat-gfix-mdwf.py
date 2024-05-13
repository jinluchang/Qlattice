#!/usr/bin/env python3

import qlat_gpt as qg
import qlat as q
import gpt as g

qg.begin_with_gpt()

q.qremove_all_info("results")
q.qmkdir_info("results")
rs = q.RngState("seed")
total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
q.displayln_info("CHECK: geo.show() =", geo.show())

gf = q.GaugeField(geo)

gf.set_rand(rs.split("gf-init"), 0.05, 2)

gf.show_info()

gt = q.GaugeTransform(geo)

gt.set_rand(rs.split("gt-init"), 0.1, 1)

gf_gfix = gt * gf

gf_gfix.show_info()

def mk_inverter_gpt(gf, qtimer):
    mobius_params = {
            "mass": 0.1,
            "M5": 1.0,
            "b": 1.0,
            "c": 0.0,
            "Ls": 8,
            "boundary_phases": [1.0, 1.0, 1.0, 1.0],
            }
    gpt_gf = qg.gpt_from_qlat(gf)
    qm = g.qcd.fermion.mobius(gpt_gf, mobius_params)
    pc = g.qcd.fermion.preconditioner
    inv = g.algorithms.inverter
    cg = inv.cg({"eps": 1e-8, "maxiter": 10000})
    slv_5d = inv.preconditioned(pc.eo2_ne(), cg)
    slv_qm = qm.propagator(slv_5d)
    inv_qm = qg.InverterGPT(inverter = slv_qm, qtimer = qtimer)
    return inv_qm

inv = mk_inverter_gpt(gf_gfix, q.Timer("py:slv_qm", True))

inv_gt = q.InverterGaugeTransform(
        inverter = mk_inverter_gpt(gf, q.Timer("py:slv_qm", True)),
        gt = gt,
        qtimer = q.Timer("py:inv_gfix", True))

src = q.Prop(geo)
src.set_rand(rs.split("src_r"))

sol = inv * src

sol1 = inv_gt * src

sol_diff = sol1.copy()
sol_diff -= sol1

q.displayln_info(f"CHECK: {sol.qnorm():.14E} {sol1.qnorm():.14E} {sol_diff.qnorm()}")

q.timer_display()

qg.end_with_gpt()

q.displayln_info(f"CHECK: finished successfully.")
