#!/usr/bin/env python3

import qlat_gpt as qg
import qlat as q
import gpt as g

qg.begin_with_gpt()

q.qremove_all_info("results")
q.qmkdir_info("results")
rs = q.RngState("seed")
geo = q.Geometry((4, 4, 4, 4), 1)
q.displayln_info("geo.show() =", geo.show())

grid = qg.mk_grid(geo)
rng = g.random("test")
gpt_gf = g.qcd.gauge.random(grid, rng, scale=0.5)

q.displayln_info(f"g.qcd.gauge.plaquette = {g.qcd.gauge.plaquette(gpt_gf):.17f}")

gpt_gf_f = g.convert(gpt_gf, g.single)

q.displayln_info(f"g.qcd.gauge.plaquette = {g.qcd.gauge.plaquette(gpt_gf_f):.17f} single precision")

gf = qg.qlat_from_gpt(gpt_gf)

gf.show_info()

mobius_params = {
        "mass": 0.08,
        "M5": 1.8,
        "b": 1.5,
        "c": 0.5,
        "Ls": 12,
        "boundary_phases": [1.0, 1.0, 1.0, 1.0],
        }

zmobius_params = {
        "mass": 0.08,
        "M5": 1.8,
        "b": 1.0,
        "c": 0.0,
        "omega": [
            0.17661651536320583 + 1j * (0.14907774771612217),
            0.23027432016909377 + 1j * (-0.03530801572584271),
            0.3368765581549033 + 1j * (0),
            0.7305711010541054 + 1j * (0),
            1.1686138337986505 + 1j * (0.3506492418109086),
            1.1686138337986505 + 1j * (-0.3506492418109086),
            0.994175013717952 + 1j * (0),
            0.5029903152251229 + 1j * (0),
            0.23027432016909377 + 1j * (0.03530801572584271),
            0.17661651536320583 + 1j * (-0.14907774771612217),
            ],
        "boundary_phases": [1.0, 1.0, 1.0, 1.0],
        }

qm = g.qcd.fermion.mobius(gpt_gf, mobius_params)

qz = g.qcd.fermion.zmobius(gpt_gf, zmobius_params)

qz_f = g.qcd.fermion.zmobius(gpt_gf_f, zmobius_params)

pc = g.qcd.fermion.preconditioner
inv = g.algorithms.inverter
cg = inv.cg({"eps": 1e-8, "maxiter": 10000})
cg_mp = inv.cg({"eps": 3e-5, "maxiter": 300})
cg_f = inv.cg({"eps": 5e-4, "maxiter": 1000})
cg_pv_f = inv.cg({"eps": 5e-4, "maxiter": 1000})
slv_5d = inv.preconditioned(pc.eo2_ne(), cg)
slv_5d_mp = inv.preconditioned(pc.eo2_ne(), cg_mp)
slv_5d_f = inv.preconditioned(pc.eo2_ne(), cg_f)
slv_5d_pv_f = inv.preconditioned(pc.eo2_ne(), cg_pv_f)

slv_qm = qm.propagator(slv_5d)
slv_qm_mp = qm.propagator(
        inv.defect_correcting(
            inv.mixed_precision(
                slv_5d_mp, g.single, g.double),
            eps=1e-8, maxiter=100))
slv_qz_f = qz.propagator(
        inv.mixed_precision(
            slv_5d_f, g.single, g.double))
slv_qm_madwf = qm.propagator(
        inv.defect_correcting(
            inv.mixed_precision(
                pc.mixed_dwf(slv_5d_f, slv_5d_pv_f, qz_f),
                g.single, g.double),
            eps=1e-8, maxiter=100))

slv_qm_timer = q.Timer("py:slv_qm", True)
slv_qm_mp_timer = q.Timer("py:slv_qm_mp", True)
slv_qz_f_timer = q.Timer("py:slv_qz_f", True)
slv_qm_madwf_timer = q.Timer("py:slv_qm_madwf", True)

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

src_mp, sol_mp, sol1_mp = test_inv(geo, slv_qm_mp, slv_qm_mp_timer)

src_f, sol_f, sol1_f = test_inv(geo, slv_qz_f, slv_qz_f_timer)

src_madwf, sol_madwf, sol1_madwf = test_inv(geo, slv_qm_madwf, slv_qm_madwf_timer)

src_f -= src
sol_f -= sol
sol1_f -= sol1

src_mp -= src
sol_mp -= sol
sol1_mp -= sol1

src_madwf -= src
sol_madwf -= sol
sol1_madwf -= sol1

q.displayln_info(f"sol_f {sol_f.qnorm()} {sol1_f.qnorm()}")
q.displayln_info(f"sol_mp {sol_mp.qnorm()} {sol1_mp.qnorm()}")
q.displayln_info(f"sol_madwf {sol_madwf.qnorm()} {sol1_madwf.qnorm()}")

q.timer_display()

qg.end_with_gpt()
