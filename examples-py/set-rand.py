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

q.begin_with_mpi(size_node_list)

rs = q.RngState("seed")
rs_prop = rs.split("prop")

total_site = q.Coordinate([ 4, 4, 4, 8, ])
geo = q.Geometry(total_site)

q.json_results_append(f"total_site = {total_site}")
q.json_results_append(f"geo = {geo}")

prop = q.Prop(geo)
prop.set_rand(rs_prop, 1.0, 0.0)
q.json_results_append(f"prop.crc32() = {prop.crc32()} ; prop.qnorm()", prop.qnorm())

n_points = 16
psel = q.PointsSelection()
psel.set_rand(total_site, n_points, rs.split("psel"))

n_per_tslice = 16
fsel = q.FieldSelection()
fsel.set_rand(geo.total_site, n_per_tslice, rs.split("fsel"))

fselc = fsel.copy()
fselc.add_psel(psel)

sc_prop = q.SelProp(fselc)
sc_prop @= prop
q.json_results_append(f"sc_prop.qnorm()", sc_prop.qnorm())

sc_prop1 = q.SelProp(fselc)
sc_prop1.set_rand(rs_prop, 1.0, 0.0)
sc_prop1 -= sc_prop
q.json_results_append(f"sc_prop1.qnorm()", sc_prop1.qnorm())

s_prop = q.SelProp(fsel)
s_prop @= prop
q.json_results_append(f"s_prop.qnorm()", s_prop.qnorm())

s_prop1 = q.SelProp(fsel)
s_prop1.set_rand(rs_prop, 1.0, 0.0)
s_prop1 -= s_prop
q.json_results_append(f"s_prop1.qnorm()", s_prop1.qnorm())

s_prop2 = q.SelProp(fsel)
s_prop2 @= sc_prop
s_prop2 -= s_prop
q.json_results_append(f"s_prop2.qnorm()", s_prop2.qnorm())

sp_prop = q.PselProp(psel)
sp_prop @= prop
q.json_results_append(f"sp_prop.qnorm()", sp_prop.qnorm())

sp_prop1 = q.PselProp(psel)
sp_prop1 @= sc_prop
sp_prop1 -= sp_prop
q.json_results_append(f"sp_prop1.qnorm()", sp_prop1.qnorm())

prop_norm = q.sqrt_field(q.qnorm_field(prop))
q.json_results_append(f"prop_norm.qnorm()", prop_norm.qnorm())

s_prop_norm = q.sqrt_field(q.qnorm_field(s_prop))
q.json_results_append(f"s_prop_norm.qnorm()", s_prop_norm.qnorm())

s_prop_norm1 = q.SelectedField(q.ElemTypeRealD, fsel)
s_prop_norm1 @= prop_norm
s_prop_norm1 -= s_prop_norm
q.json_results_append(f"s_prop_norm1.qnorm()", s_prop_norm1.qnorm())

sp_prop_norm = q.sqrt_field(q.qnorm_field(sp_prop))
q.json_results_append(f"sp_prop_norm.qnorm()", sp_prop_norm.qnorm())

sp_prop_norm1 = q.SelectedPoints(q.ElemTypeRealD, psel)
sp_prop_norm1 @= prop_norm
sp_prop_norm1 -= sp_prop_norm
q.json_results_append(f"sp_prop_norm1.qnorm()", sp_prop_norm1.qnorm())

q.check_log_json(__file__, check_eps=1e-14)

q.timer_display()

q.end_with_mpi()

q.displayln_info(f"CHECK: finished successfully.")
