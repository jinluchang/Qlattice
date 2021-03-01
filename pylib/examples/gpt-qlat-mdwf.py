#!/usr/bin/env python3

import qlat_gpt as qg
import qlat as q
import gpt as g

qg.begin_with_gpt()

q.qremove_all_info("results")
q.qmkdir_info("results")
rs = q.RngState("seed")
geo = q.Geometry((4, 4, 4, 8), 1)
q.displayln_info("geo.show() =", geo.show())

gf = q.GaugeField(geo)

gf.set_rand(rs.split("gf-init"), 0.05, 2)

gf.show_info()

gpt_gf = qg.gpt_from_qlat(gf)

q.displayln_info("g.qcd.gauge.plaquette = {g.qcd.gauge.plaquette(gpt_gf):.17f}")

gf1 = qg.qlat_from_gpt(gpt_gf)

gf1.show_info()

mobius_params = {
        "mass": 0.1,
        "M5": 1.0,
        "b": 1.0,
        "c": 0.0,
        "Ls": 8,
        "boundary_phases": [1.0, 1.0, 1.0, 1.0],
        }

qm = g.qcd.fermion.mobius(gpt_gf, mobius_params)

pc = g.qcd.fermion.preconditioner
inv = g.algorithms.inverter
cg = inv.cg({"eps": 1e-10, "maxiter": 100})
slv_5d = inv.preconditioned(pc.eo2_ne(), cg)
slv_qm = qm.propagator(slv_5d)

grid = qg.mk_grid(geo)
src = g.mspincolor(grid)
g.create.point(src, [0, 0, 0, 0])
dst_qm = g.mspincolor(grid)
dst_qm @= slv_qm * src

prop = qg.qlat_from_gpt(dst_qm)

ld = q.contract_pion_field(prop, 0)

q.displayln_info(q.show(ld))

q.timer_display()

qg.end_with_gpt()
