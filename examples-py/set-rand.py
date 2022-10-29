#!/usr/bin/env python3

import sys
import qlat as q

size_node_list = [
        [1, 1, 1, 1],
        [1, 1, 1, 2],
        [1, 1, 1, 4],
        [1, 1, 1, 8],
        [2, 2, 2, 2],
        [2, 2, 2, 4],
        ]

q.begin(sys.argv, size_node_list)

q.qremove_all_info("results")
q.qmkdir_info("results")

rs = q.RngState("seed")

total_site = [ 4, 4, 4, 8, ]
geo = q.Geometry(total_site, 1)

q.displayln_info(f"CHECK: total_site = {total_site}")
q.displayln_info(f"CHECK: geo = {geo}")

rs_prop = rs.split("prop")

prop = q.Prop(geo)
prop.set_rand(rs_prop, 1.0, 0.0)

q.displayln_info(f"CHECK: prop.crc32() = {prop.crc32()} ; prop.qnorm() = {prop.qnorm()}")

n_points = 16
psel = q.PointSelection()
psel.set_rand(rs.split("psel"), total_site, n_points)

n_per_tslice = 16
fsel = q.FieldSelection(geo.total_site(), n_per_tslice, rs.split("fsel"))

fselc = fsel.copy()
fselc.add_psel(psel)

sc_prop = q.SelProp(fselc)
sc_prop @= prop
q.displayln_info(f"CHECK: sc_prop.qnorm() = {sc_prop.qnorm()}")

sc_prop1 = q.SelProp(fselc)
sc_prop1.set_rand(rs_prop, 1.0, 0.0)
sc_prop1 -= sc_prop
q.displayln_info(f"CHECK: sc_prop1.qnorm() = {sc_prop1.qnorm()}")

s_prop = q.SelProp(fsel)
s_prop @= prop
q.displayln_info(f"CHECK: s_prop.qnorm() = {s_prop.qnorm()}")

s_prop1 = q.SelProp(fsel)
s_prop1.set_rand(rs_prop, 1.0, 0.0)
s_prop1 -= s_prop
q.displayln_info(f"CHECK: s_prop1.qnorm() = {s_prop1.qnorm()}")

s_prop2 = q.SelProp(fsel)
s_prop2 @= sc_prop
s_prop2 -= s_prop
q.displayln_info(f"CHECK: s_prop2.qnorm() = {s_prop2.qnorm()}")

q.timer_display()

q.displayln_info(f"CHECK: finished successfully.")

q.end()
