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

slv_qm_timer = q.Timer("py:slv_qm", True)

def mk_src(geo):
    src = q.mk_point_src(geo, [0, 0, 0, 0])
    grid = qg.mk_grid(geo)
    g_src = g.mspincolor(grid)
    g.create.point(g_src, [0, 0, 0, 0])
    src1 = qg.qlat_from_gpt(g_src)
    src1 -= src
    assert src1.qnorm() == 0.0
    return src

def test_inv(geo, slv, slv_timer = None):
    src = mk_src(geo)
    q.displayln_info(f"src info {src.qnorm()} {src.crc32()}")
    sol = qg.qlat_invert(src, slv, slv_timer)
    q.displayln_info(f"sol info {sol.qnorm()} {sol.crc32()}")
    sol1 = qg.qlat_invert(sol, slv, slv_timer)
    q.displayln_info(f"sol1 info {sol1.qnorm()} {sol1.crc32()}")
    return src, sol, sol1

src, sol, sol1 = test_inv(geo, slv_qm, slv_qm_timer)

ld = q.contract_pion_field(sol, 0)

q.displayln_info("q.contract_pion_field(sol)", q.show(ld))

ld1 = q.contract_pion_field(sol1, 0)

q.displayln_info("q.contract_pion_field(sol1)", q.show(ld1))

q.timer_display()

qg.end_with_gpt()
