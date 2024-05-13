#!/usr/bin/env python3

import qlat_gpt as qg
import qlat as q
import gpt as g

@q.timer
def test_gf(gf):
    q.displayln_info(f"CHECK: test_gf")
    #
    plaq = gf.plaq()
    #
    gpt_gf = qg.gpt_from_qlat(gf)
    g_plaq = g.qcd.gauge.plaquette(gpt_gf)
    #
    q.displayln_info(f"CHECK: gf.plaq()={plaq:.13E} ; g.qcd.gauge.plaquette(gpt_gf) = {g_plaq:.13F}")
    assert abs(plaq - g_plaq) < 1.0e-10

@q.timer
def test_src(xg):
    q.displayln_info(f"CHECK: test_src: xg={xg.to_list()}")
    #
    src_q = q.mk_point_src(geo, xg)
    #
    grid = qg.mk_grid(geo)
    g_src_gpt = g.mspincolor(grid)
    g.create.point(g_src_gpt, xg.to_list())
    src_g = qg.qlat_from_gpt(g_src_gpt)
    #
    src_diff = src_q.copy()
    src_diff -= src_g
    #
    q.displayln_info(f"CHECK: {src_q.qnorm():.12E} {src_g.qnorm():.12E} {src_diff.qnorm():.12E}")
    assert src_diff.qnorm() == 0

qg.begin_with_gpt()

rs = q.RngState("seed")
total_site = q.Coordinate([4, 4, 4, 8])
geo = q.Geometry(total_site)
q.displayln_info("CHECK: geo.show() =", geo.show())

gf = q.GaugeField(geo)
gf.set_rand(rs, 0.2, 2)
gf.show_info()

test_gf(gf)

for i in range(16):
    xg = q.Coordinate([ rs.rand_gen() % total_site[0],
                       rs.rand_gen() % total_site[1],
                       rs.rand_gen() % total_site[2],
                       rs.rand_gen() % total_site[3] ])
    test_src(xg)

q.timer_display()

qg.end_with_gpt()

q.displayln_info(f"CHECK: finished successfully.")
