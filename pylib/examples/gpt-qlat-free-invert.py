#!/usr/bin/env python3

import qlat_gpt as qg
import qlat as q
import gpt as g

qg.begin_with_gpt()

q.qremove_all_info("results")
q.qmkdir_info("results")
rs = q.RngState("seed")
geo = q.Geometry([4, 4, 4, 8], 1)
q.displayln_info("geo.show() =", geo.show())

gf = q.GaugeField(geo)
gf.set_unit()
gf.show_info()

mobius_params = {
        "mass": 0.05,
        "M5": 1.0,
        "b": 1.0,
        "c": 0.0,
        "Ls": 16,
        "boundary_phases": [1.0, 1.0, 1.0, 1.0],
        }

src = q.mk_point_src(geo, [0, 0, 0, 0])

sol = q.free_invert(src, mobius_params["mass"], mobius_params["M5"])

gpt_gf = qg.gpt_from_qlat(gf)
qm = g.qcd.fermion.mobius(gpt_gf, mobius_params)
pc = g.qcd.fermion.preconditioner
inv = g.algorithms.inverter
cg = inv.cg({"eps": 1e-8, "maxiter": 10000})
slv_5d = inv.preconditioned(pc.eo2_ne(), cg)
slv_qm = qm.propagator(slv_5d)

sol1 = qg.qlat_invert(src, slv_qm)
sol1 -= sol

q.displayln_info(sol.qnorm(), sol1.qnorm())
