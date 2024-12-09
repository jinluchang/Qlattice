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

grid = qg.mk_grid(geo)
rng = g.random("test")
gpt_gf = g.qcd.gauge.random(grid, rng, scale=0.5)

q.displayln_info(f"CHECK: g.qcd.gauge.plaquette = {g.qcd.gauge.plaquette(gpt_gf):.10f}")

gpt_gf_f = g.convert(gpt_gf, g.single)

q.displayln_info(f"CHECK: g.qcd.gauge.plaquette = {g.qcd.gauge.plaquette(gpt_gf_f):.4f} single precision")

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
cg_split = inv.split(cg_mp, mpi_split = g.default.get_ivec("--mpi_split", None, 4))
cg_f = inv.cg({"eps": 5e-4, "maxiter": 1000})
cg_pv_f = inv.cg({"eps": 5e-4, "maxiter": 1000})

slv_5d = inv.preconditioned(pc.eo2_ne(), cg)
slv_5d_mp = inv.preconditioned(pc.eo2_ne(), cg_mp)
slv_5d_split = inv.preconditioned(pc.eo2_ne(), cg_split)
slv_5d_f = inv.preconditioned(pc.eo2_ne(), cg_f)
slv_5d_pv_f = inv.preconditioned(pc.eo2_ne(), cg_pv_f)

slv_qm = qm.propagator(slv_5d)

slv_qz_f = qz.propagator(
        inv.mixed_precision(
            slv_5d_f, g.single, g.double)).grouped(4)

slv_qm_mp = qm.propagator(
        inv.defect_correcting(
            inv.mixed_precision(
                slv_5d_mp, g.single, g.double),
            eps=1e-8, maxiter=100)).grouped(4)

slv_qm_split_sloppy = qm.propagator(
        inv.defect_correcting(
            inv.mixed_precision(
                slv_5d_split, g.single, g.double),
            eps=1e-8, maxiter=1)).grouped(4)

slv_qm_split = qm.propagator(
        inv.defect_correcting(
            inv.mixed_precision(
                slv_5d_split, g.single, g.double),
            eps=1e-8, maxiter=100)).grouped(4)

slv_qm_madwf = qm.propagator(
        inv.defect_correcting(
            inv.mixed_precision(
                pc.mixed_dwf(slv_5d_f, slv_5d_pv_f, qz_f),
                g.single, g.double),
            eps=1e-8, maxiter=100)).grouped(4)

inv_qm = qg.InverterGPT(inverter = slv_qm, qtimer = q.Timer("py:slv_qm", True))
inv_qz_f = qg.InverterGPT(inverter = slv_qz_f, qtimer = q.Timer("py:slv_qz_f", True))
inv_qm_mp = qg.InverterGPT(inverter = slv_qm_mp, qtimer = q.Timer("py:slv_qm_mp", True))
inv_qm_split = qg.InverterGPT(inverter = slv_qm_split, qtimer = q.Timer("py:slv_qm_split", True))
inv_qm_split_sloppy = qg.InverterGPT(inverter = slv_qm_split_sloppy, qtimer = q.Timer("py:slv_qm_split_sloppy", True))
inv_qm_madwf = qg.InverterGPT(inverter = slv_qm_madwf, qtimer = q.Timer("py:slv_qm_madwf", True))

def mk_src(geo):
    src = q.mk_point_src(geo, q.Coordinate([0, 0, 0, 0]))
    grid = qg.mk_grid(geo)
    g_src = g.mspincolor(grid)
    g.create.point(g_src, [0, 0, 0, 0])
    src1 = qg.qlat_from_gpt(g_src)
    src1 -= src
    assert src1.qnorm() == 0.0
    return src

def test_inv(geo, inverter):
    src = mk_src(geo)
    q.displayln_info(f"CHECK: src info {src.qnorm()}")
    q.displayln_info(f"CHECK: src info {src.crc32():08X}")
    sol = inverter * src
    q.displayln_info(f"CHECK: sol info {sol.qnorm():.7E}")
    q.displayln_info(f"sol info {sol.crc32():08X}")
    sol1 = inverter * sol
    q.displayln_info(f"CHECK: sol1 info {sol1.qnorm():.4E}")
    q.displayln_info(f"sol1 info {sol1.crc32():08X}")
    return src, sol, sol1

tags = [ "qm", "qz_f", "qm_mp", "qm_split", "qm_split_sloppy", "inv_qm_madwf" ]
invs = [ inv_qm, inv_qz_f, inv_qm_mp, inv_qm_split, inv_qm_split_sloppy, inv_qm_madwf]

q.displayln_info(f"CHECK: tag={tags[0]} start")
src, sol, sol1 = test_inv(geo, invs[0])

for tag, inv in zip(tags[1:], invs[1:]) :
    q.displayln_info(f"CHECK: tag={tag} start")
    src_n, sol_n, sol1_n = test_inv(geo, inv)
    src_n -= src
    sol_n -= sol
    sol1_n -= sol1
    q.displayln_info(f"CHECK: tag={tag} diff src {src_n.qnorm()} sol {sol_n.qnorm():.1E} sol1 {sol1_n.qnorm():.1E}")

q.timer_display()

qg.end_with_gpt()

q.displayln_info(f"CHECK: finished successfully.")
